#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#define N 6

double complex A[N*N];
double complex copyA[N*N];

double complex Gij[N*N];
double complex oldGij[N*N];

double complex parn;
double complex parm;

double complex oldparn;
double complex oldparm;

double w[N];
double complex work[3*N];
double rwork[3*N];
double doublekmax=2.0*M_PI;
int Niter=500;
double mixing=0.9;
double orderpara_eps=1.0e-12;

double Qx=2.0*M_PI/3.0;
double Qy=2.0*M_PI/3.0;
int RBZ=3;// reduced BZ
int intkmax=6;// L=intkmax, Ns=L^2
double filling=0.5;// 1/2 filling

double t=1.0;
double U=0.0;
double V=0.0;

double complex tmp;
double tmp_err;
double shift;

int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);

struct myData {
  double data;
  double complex vec[N];
  double px;
  double py;
  int orig_pos; // will hold the original position of data on your array
};

int myData_compare (const void* a, const void* b) {
  double xx = ((struct myData*)a)->data;
  double yy = ((struct myData*)b)->data;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

int main(int argc, char *argv[])
{
  int c;
  int show_band = 0;

  if (argv[optind] == NULL || argv[optind + 1] == NULL) {
    fprintf(stderr,"arguments missing\n");
    fprintf(stderr,"usage: %s -U [U] -V [V] -k [kmax] -i [Niter] -m [mix] -b (show band)\n",argv[0]);
    exit(1);
  }
  while((c = getopt(argc, argv, "U:V:k:i:m:b")) != -1){
    switch (c){
      case 'U':
        U = strtod(optarg,NULL);
        break;
      case 'V':
        V = strtod(optarg,NULL);
        break;
      case 'k':
        intkmax = strtol(optarg,NULL,10);
        break;
      case 'i':
        Niter = strtol(optarg,NULL,10);
        break;
      case 'm':
        mixing = strtod(optarg,NULL);
        break;
      case 'b':
        show_band = 1;
        break;
      default:
        fprintf(stderr,"usage: %s -U [U] -V [V] -k [kmax] -i [Niter] -m [mix] -b (show band)\n",argv[0]);
        exit(1);
    }
  }

  int i,j,k;
  char jobz = 'V';
  char uplo ='U';
  int n=N,lda=N,lwork=3*N,info;

  int intkx,intky;
  int numk;
  int Norbital=N*intkmax*intkmax/RBZ;// 2*3*L^2/3 = 2*L^2 = (# eigenvalue)
  int Ne=(int)Norbital*filling;
  int Nfourier=intkmax*intkmax;
  double invintkmax=1.0/intkmax;
  double kx,ky;
  double Ene=0.0;
  double Ene_0=0.0;
  double Ene_U=0.0;

  int Ndiv=intkmax;
  int intdos[Ndiv];
  double Efermi=0.0;
  double Emin;
  double Emax;
  double Edelta;
  double dos[Ndiv];

  struct myData* array = (struct myData*) malloc(Norbital*sizeof(struct myData));

  double complex **Epsilon;
  Epsilon = (double complex**)malloc((intkmax)*sizeof(double complex*));
  Epsilon[0] = (double complex*)malloc((intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(intkmax); i++){
    Epsilon[i] = Epsilon[0] + i*(intkmax);
  }
  double complex **EpsilonQ;
  EpsilonQ = (double complex**)malloc((intkmax)*sizeof(double complex*));
  EpsilonQ[0] = (double complex*)malloc((intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(intkmax); i++){
    EpsilonQ[i] = EpsilonQ[0] + i*(intkmax);
  }
  double complex **Epsilon2Q;
  Epsilon2Q = (double complex**)malloc((intkmax)*sizeof(double complex*));
  Epsilon2Q[0] = (double complex*)malloc((intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(intkmax); i++){
    Epsilon2Q[i] = Epsilon2Q[0] + i*(intkmax);
  }

  int step;

  for(i=0; i<N*N; i++){
    oldGij[i]=0.0;
    oldGij[i]=0.0;
  }

// half fill (U=0)
  parn = 0.5;
  parm = 0.5;

  oldparn = 0.0;
  oldparm = 0.0;

  printf("# filling: %f\n",filling);
  printf("# t: %f\n",t);
  printf("# U: %f\n",U);
  printf("# V: %f\n",V);
  printf("# intkmax: %d\n",intkmax);
  printf("# mixing: %f\n",mixing);
  printf("# orderpara_eps: %e\n",orderpara_eps);
  printf("\n\n");
  printf("# ene ene_0 ene_U\n");

  for(intkx=0; intkx<intkmax; intkx++){
    for(intky=0; intky<intkmax; intky++){
      kx = doublekmax*((intkx+0.0)*invintkmax-0.5);
//      kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
      ky = doublekmax*((intky+0.0)*invintkmax-0.5);
      Epsilon[intkx][intky] = cos(kx)+cos(ky)+cos(kx+ky);
      EpsilonQ[intkx][intky] = cos(kx+Qx)+cos(ky+Qy)+cos(kx+Qx+ky+Qy);
      Epsilon2Q[intkx][intky] = cos(kx+2*Qx)+cos(ky+2*Qy)+cos(kx+2*Qx+ky+2*Qy);
    }
  }

for(step=0; step<Niter; step++){

  numk=0;
  shift = doublekmax*0.25*invintkmax;
    for(intkx=0; intkx<intkmax; intkx++){
      for(intky=0; intky<intkmax; intky++){
        kx = doublekmax*((intkx+0.0)*invintkmax-0.5);
//        kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
        ky = doublekmax*((intky+0.0)*invintkmax-0.5);

        if(
             kx+ky + shift < doublekmax/3.0
          && kx+ky + shift > -doublekmax/3.0
          && kx-2*ky + shift > -doublekmax*5.0/6.0
          && kx-2*ky + shift < doublekmax*5.0/6.0
          && 2*kx-ky + shift > -doublekmax*5.0/6.0
          && 2*kx-ky + shift < doublekmax*5.0/6.0
        ){// reduced BZ

          for(i=0; i<N*N; i++){
            A[i]=0.0;
            copyA[i]=0.0;
          }

//  0x  1o  2o  3o  4o  5o
//  6x  7x  8o  9o 10o 11o
// 12x 13x 14x 15o 16o 17o
// 18x 19x 20x 21x 22o 23o
// 24x 25x 26x 27x 28x 29o
// 30x 31x 32x 33x 34x 35x

          A[0] = U*parn + (-2.0*t) * Epsilon[intkx][intky];
          A[7] = U*parn + (-2.0*t) * EpsilonQ[intkx][intky];
          A[14] = U*parn + (-2.0*t) * Epsilon2Q[intkx][intky];
          A[21] = A[0];
          A[28] = A[7];
          A[35] = A[14];

          A[20] = -U*parm;
          A[24] = A[20];
          A[31] = A[20];

          A[1]=conj(A[6]);
          A[2]=conj(A[12]);
          A[3]=conj(A[18]);
          A[4]=conj(A[24]);
          A[5]=conj(A[30]);
          A[8]=conj(A[13]);
          A[9]=conj(A[19]);
          A[10]=conj(A[25]);
          A[11]=conj(A[31]);
          A[15]=conj(A[20]);
          A[16]=conj(A[26]);
          A[17]=conj(A[32]);
          A[22]=conj(A[27]);
          A[23]=conj(A[33]);
          A[29]=conj(A[34]);

          for(i=0; i<N; i++){
            for(j=0; j<N; j++){
              copyA[i+j*N]=A[i*N+j];
            }
          }

          zheev_(&jobz, &uplo, &n, copyA, &lda, w, work, &lwork, rwork, &info);

//printf("%d %f %f %d\n",numk,kx/doublekmax*intkmax,ky/doublekmax*intkmax,Norbital);
          for(i=0; i<N; i++){
            array[i + numk*N].data = w[i];
            array[i + numk*N].orig_pos = i + numk*N;
            array[i + numk*N].px = kx;
            array[i + numk*N].py = ky;
            for(j=0; j<N; j++){
              array[i + numk*N].vec[j] = copyA[i*N+j];
            }
          }

          numk++;
if(numk == intkmax*intkmax/RBZ + 1){
  printf("# numk=%d > intkmax*intkmax/RBZ=%d\n",numk,intkmax*intkmax/RBZ);
  printf("# failed to take correct RBZ\n");
  exit(2);
}

        }
      }
    }

  qsort(array, Norbital, sizeof(struct myData), &myData_compare);

  Ene_0=0.0;
  for(i=0; i<Ne; i++){
    Ene_0 += array[i].data;
  }

/*
  for(i=0; i<Norbital; i++){
//printf("%d %d %f\n",i,array[i].sigma,array[i].data);
printf("%d %d %f %f %f\n",i,array[i].sigma,array[i].px,array[i].py,array[i].data);
  }
*/

  Ene_U = 0.0;
  Ene_U += -U*Nfourier*(conj(parn)*parn - conj(parm)*parm);

  Ene = Ene_0 + Ene_U;

  printf("%6d %13.10f ",step,Ene/Nfourier);
  printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
  printf("%13.10f ",creal(parn));
  printf("%13.10f ",creal(parm));
  printf("\n");

  for(i=0; i<N*N; i++){
    Gij[i]=0.0;
    Gij[i]=0.0;
  }
  for(k=0; k<Ne; k++){
    for(i=0; i<N; i++){
      tmp = conj(array[k].vec[i]);
      for(j=0; j<N; j++){
        Gij[i*N+j] += tmp * array[k].vec[j];
      }
    }
  }

  parn = (Gij[0*N+0]+Gij[1*N+1]+Gij[2*N+2]+Gij[3*N+3]+Gij[4*N+4]+Gij[5*N+5])/(2.0*Nfourier);
  parm = (Gij[2*N+3]+Gij[0*N+4]+Gij[1*N+5]+Gij[3*N+2]+Gij[4*N+0]+Gij[5*N+1])/(2.0*Nfourier);

//  parn = (Gij[0*N+0]+Gij[1*N+1]+Gij[2*N+2])/(Nfourier);
//  parm = (Gij[2*N+3]+Gij[0*N+4]+Gij[1*N+5])/(Nfourier);

  tmp_err = 0.0;
  for(i=0; i<N*N; i++){
    tmp_err += cabs(oldGij[i] - Gij[i])*cabs(oldGij[i] - Gij[i]);
  }
  if(sqrt(tmp_err)/Nfourier < orderpara_eps
    && fabs(creal(parn - oldparn)) < orderpara_eps
    && fabs(creal(parm - oldparm)) < orderpara_eps
  ){
    break;
  }

  parn = parn*mixing + oldparn*(1.0-mixing);
  parm = parm*mixing + oldparm*(1.0-mixing);

  for(i=0; i<N*N; i++){
    oldGij[i] = Gij[i];
  }
  oldparn = parn;
  oldparm = parm;
}

  if(step==Niter){
    printf("\n\n");
    printf("# not converged within %d steps !!!\n",Niter);
    printf("# ");
/*
  }else if(
    creal(Nc[0]+DeltaNc[0])>1.0+1.0e10 ||
    creal(Nc[0]+DeltaNc[0])<0.0-1.0e10 ||
    creal(Nc[0]-DeltaNc[0])>1.0+1.0e10 ||
    creal(Nc[0]-DeltaNc[0])<0.0-1.0e10 ||
    creal(Nc[1]+DeltaNc[1])>1.0+1.0e10 ||
    creal(Nc[1]+DeltaNc[1])<0.0-1.0e10 ||
    creal(Nc[1]-DeltaNc[1])>1.0+1.0e10 ||
    creal(Nc[1]-DeltaNc[1])<0.0-1.0e10 ||
    creal(Nf[0]+DeltaNf[0])>1.0+1.0e10 ||
    creal(Nf[0]+DeltaNf[0])<0.0-1.0e10 ||
    creal(Nf[0]-DeltaNf[0])>1.0+1.0e10 ||
    creal(Nf[0]-DeltaNf[0])<0.0-1.0e10 ||
    creal(Nf[1]+DeltaNf[1])>1.0+1.0e10 ||
    creal(Nf[1]+DeltaNf[1])<0.0-1.0e10 ||
    creal(Nf[1]-DeltaNf[1])>1.0+1.0e10 ||
    creal(Nf[1]-DeltaNf[1])<0.0-1.0e10
  ){
    printf("\n\n");
    printf("# converged to <n> < 0 or <n> >1 (%d steps) !!! \n",step);
    printf("# ");
*/
  }else{
    printf("\n\n");
    printf("# converged in %d steps\n",step);
  }
  printf("%6d %13.10f ",step,Ene/Nfourier);
  printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
  printf("%13.10f ",creal(parn));
  printf("%13.10f ",creal(parm));
  printf("\n");

/*
  printf("\n\n");
  printf("# SzcA SzfA SzcB SzfB NcA NfA NcB NfB\n");
  printf("%13.10f ",creal(Nc[0]+DeltaNc[0]-Nc[1]-DeltaNc[1])*0.5);// SzcA
  printf("%13.10f ",creal(Nf[0]+DeltaNf[0]-Nf[1]-DeltaNf[1])*0.5);// SzfA
  printf("%13.10f ",creal(Nc[0]-DeltaNc[0]-Nc[1]+DeltaNc[1])*0.5);// SzcB
  printf("%13.10f ",creal(Nf[0]-DeltaNf[0]-Nf[1]+DeltaNf[1])*0.5);// SzfB
  printf("%13.10f ",creal(Nc[0]+DeltaNc[0]+Nc[1]+DeltaNc[1]));// NcA
  printf("%13.10f ",creal(Nf[0]+DeltaNf[0]+Nf[1]+DeltaNf[1]));// NfA
  printf("%13.10f ",creal(Nc[0]-DeltaNc[0]+Nc[1]-DeltaNc[1]));// NcB
  printf("%13.10f ",creal(Nf[0]-DeltaNf[0]+Nf[1]-DeltaNf[1]));// NfB
  printf("\n");
*/

  if(fabs(array[Ne-1].data-array[Ne].data) < 1.0e-16){
    printf("\n\n");
    printf("# !!! OPEN SHELL !!!\n");
  }

if(show_band==1){
/*
  printf("# tb1=%f tb2=%f tp=%f tq=%f\n",tb1,tb2,tp,tq);
  printf("\n\n");
  printf("# 1: band\n");
  printf("# number spin intkx intky kx ky band energy orig_posit\n");
  numk=0;
  for(spin=0; spin<2; spin++){
    for(intkx=0; intkx<intkmax; intkx++){
      for(intky=0; intky<intkmax; intky++){
        kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
        ky = doublekmax*((intky+0.0)*invintkmax-0.5);
//
        if(fabs(kx)+fabs(ky) < doublekmax*0.5+1.0e-10){// reduced BZ for Q=(pi,pi)
          for(i=0; i<N; i++){
            printf("%4d %2d %3d %3d %13.10f %13.10f %2d ",i+numk*N,spin,intkx,intky,kx,ky,i);
            printf("%13.10f %4d ",array[i+numk*N].data,array[i+numk*N].orig_pos);
            printf("[");
            for(j=0; j<N; j++){
              printf("%13.10f+i%13.10f ",creal(array[i+numk*N].vec[j]),cimag(array[i+numk*N].vec[j]));
            }
            printf("]\n");
          }
          numk++;
        }
      }
    }
  }
*/
/*
  printf("\n\n");
  printf("# 2: energy\n");
  printf("# i Eigenvalue orig_posit\n");
  for(i=0; i<Norbital; i++){
    printf("%4d %13.10f %4d ",i,array[i].data,array[i].orig_pos);
    printf("[");
    for(j=0; j<N; j++){
      printf("%13.10f+i%13.10f ",creal(array[i].vec[j]),cimag(array[i].vec[j]));
    }
    printf("]\n");
  }
*/
  Efermi = array[Ne-1].data;
  printf("\n\n");
  printf("# 3: dos\n");
  printf("# Efermi=%f\n",Efermi);
  printf("# Egap=%20.17f\n",array[Ne].data-Efermi);
  Emin=array[0].data-(array[Norbital-1].data-array[0].data)/Ndiv*2.0;
  Emax=array[Norbital-1].data+(array[Norbital-1].data-array[0].data)/Ndiv*2.0;
  Edelta=(Emax-Emin)/Ndiv;
  for(i=0; i<Ndiv; i++){
    intdos[i]=0;
  }
  for(i=0; i<Norbital; i++){
    j=(int)((array[i].data-Emin)/Edelta);
    intdos[j]++;
  }
  printf("\n");
  for(i=0; i<Ndiv; i++){
    dos[i]=intdos[i]/(double)(Norbital*Edelta);
    printf("%d %f %d %f\n",i,Emin+i*Edelta-Efermi,intdos[i],dos[i]);
  }
}

  free(array);
  free(Epsilon[0]);
  free(Epsilon);
  free(EpsilonQ[0]);
  free(EpsilonQ);
  free(Epsilon2Q[0]);
  free(Epsilon2Q);

  return 0;
}
