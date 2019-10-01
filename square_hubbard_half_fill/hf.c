#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#define N 2 // momenta k and k+Q
#define RBZ 2 // reduced BZ
#define Nspin 2 // spin up and down

double complex A[N*N];
double complex copyA[N*N];

double complex Guij[N*N];
double complex Gdij[N*N];
double complex oldGuij[N*N];
double complex oldGdij[N*N];

double complex parN[2];
double complex parDeltaN[2];
double complex oldparN[2];
double complex oldparDeltaN[2];

double w[N];
double complex work[3*N];
double rwork[3*N];
double doublekmax=2.0*M_PI;
int Niter=500;
double mixing=0.9;
double orderpara_eps=1.0e-12;

double Qx=M_PI;
double Qy=M_PI;
int intkmax=48;// L=intkmax, Ns=L^2
double filling=0.5;// 1/2 filling

double t=1.0;
double U=0.0;

double complex tmp;
double tmp_err;

int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);

struct myData {
  double data;
  double complex vec[N];
  double px;
  double py;
  int sigma;
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
  int init_cond = 0;
  int show_interval = 1;
  int show_band = 0;

  if (argv[optind] == NULL || argv[optind + 1] == NULL) {
    fprintf(stderr,"arguments missing\n");
    fprintf(stderr,"usage: %s -U [U] -k [kmax] -i [Niter] -m [mix] -f [filling] -c [init_cond] -s [show_interval] -b (show band)\n",argv[0]);
    exit(1);
  }
  while((c = getopt(argc, argv, "U:k:i:m:f:c:s:b")) != -1){
    switch (c){
      case 'U':
        U = strtod(optarg,NULL);
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
      case 'f':
        filling = strtod(optarg,NULL);
        break;
      case 'c':
        init_cond = strtol(optarg,NULL,10);
        break;
      case 's':
        show_interval = strtol(optarg,NULL,10);
        break;
      case 'b':
        show_band = 1;
        break;
      default:
        fprintf(stderr,"usage: %s -U [U] -k [kmax] -i [Niter] -m [mix] -f [filling] -c [init_cond] -s [show_interval] -b (show band)\n",argv[0]);
        exit(1);
    }
  }

  int i,j,k;
  char jobz = 'V';
  char uplo ='U';
  int n=N,lda=N,lwork=3*N,info;

  int intkx,intky;
  int numk;
  int spin;
  int Norbital=Nspin*N*intkmax*intkmax/RBZ;// 2*2*L^2/2 = 2*L^2 = (# eigenvalue)
  int Ne=(int)Norbital*filling;
  int Nfourier=intkmax*intkmax;
  double invintkmax=1.0/intkmax;
  double kx,ky;
  double Ene=0.0;
  double Ene_0=0.0;
  double Ene_U=0.0;

  int Ndiv=intkmax;
  int intdos[Ndiv];
  int intdosu[Ndiv];
  int intdosd[Ndiv];
  double Efermi=0.0;
  double Emin;
  double Emax;
  double Edelta;
  double dos[Ndiv];
  double dosu[Ndiv];
  double dosd[Ndiv];

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

  int step;

  for(i=0; i<N*N; i++){
    oldGuij[i]=0.0;
    oldGdij[i]=0.0;
  }

  if(init_cond == 1){
    // AF
    parN[0] = filling;
    parN[1] = filling;
    parDeltaN[0] = filling;
    parDeltaN[1] = -filling;
  }else if(init_cond == 2){
    // FM
    parN[0] = filling*2.0;
    parN[1] = 0.0;
    parDeltaN[0] = 0.0;
    parDeltaN[1] = 0.0;
  }else{
    // nonmag
    parN[0] = filling;
    parN[1] = filling;
    parDeltaN[0] = 0.0;
    parDeltaN[1] = 0.0;
  }

  oldparN[0] = parN[0];
  oldparN[1] = parN[1];
  oldparDeltaN[0] = parDeltaN[0];
  oldparDeltaN[1] = parDeltaN[1];

  printf("# filling: %f\n",filling);
  printf("# t: %f\n",t);
  printf("# U: %f\n",U);
  printf("# intkmax: %d\n",intkmax);
  printf("# init_cond: %d\n",init_cond);
  printf("# mixing: %f\n",mixing);
  printf("# orderpara_eps: %e\n",orderpara_eps);
  printf("\n\n");
//// parN should be filling=0.5
//// parDeltaN could be nonzero
  printf("# ene ene_0 ene_U ");
  printf("Ave(N) N[0] N[1] Ave(DeltaN) DeltaN[0] DeltaN[1] ");
  printf("\n");

  for(intkx=0; intkx<intkmax; intkx++){
    for(intky=0; intky<intkmax; intky++){
      kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
      ky = doublekmax*((intky+0.0)*invintkmax-0.5);
      Epsilon[intkx][intky] = -2.0*t*(cos(kx)+cos(ky));
      EpsilonQ[intkx][intky] = -2.0*t*(cos(kx+Qx)+cos(ky+Qy));
    }
  }

for(step=0; step<Niter; step++){

  numk=0;
  for(spin=0; spin<2; spin++){
    for(intkx=0; intkx<intkmax; intkx++){
      for(intky=0; intky<intkmax; intky++){
        kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
        ky = doublekmax*((intky+0.0)*invintkmax-0.5);
        if(fabs(kx)+fabs(ky) < doublekmax*0.5+1.0e-10){// reduced BZ for Q=(pi,pi)

          for(i=0; i<N*N; i++){
            A[i]=0.0;
            copyA[i]=0.0;
          }

          //  0x  1o
          //  2x  3x
          A[0] = Epsilon[intkx][intky] + U*parN[1-spin];
          A[3] = EpsilonQ[intkx][intky] + U*parN[1-spin];
          A[2] = U*parDeltaN[1-spin];
          // copy h.c.
          A[1]=conj(A[2]);

          for(i=0; i<N; i++){
            for(j=0; j<N; j++){
              copyA[i+j*N]=A[i*N+j];
            }
          }

          zheev_(&jobz, &uplo, &n, copyA, &lda, w, work, &lwork, rwork, &info);

          for(i=0; i<N; i++){
            array[i + numk*N].data = w[i];
//printf("%d %d %f %f %f\n",i,spin,kx,ky,w[i]);
            array[i + numk*N].orig_pos = i + numk*N;
            array[i + numk*N].sigma = spin;
            array[i + numk*N].px = kx;
            array[i + numk*N].py = ky;
            for(j=0; j<N; j++){
              array[i + numk*N].vec[j] = copyA[i*N+j];
            }
          }

          numk++;
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
  Ene_U += -U*Nfourier*(
    + conj(parN[0])*parN[1]
    + conj(parDeltaN[0])*parDeltaN[1]
    );
  Ene = Ene_0 + Ene_U;

if(step%show_interval == 0){
  printf("%6d %13.10f ",step,Ene/Nfourier);
  printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
  printf("%13.10f ",creal(parN[0]+parN[1])*0.5);
  printf("%13.10f ",creal(parN[0]));
  printf("%13.10f ",creal(parN[1]));
  printf("%13.10f ",creal(parDeltaN[0]+parDeltaN[1])*0.5);
  printf("%13.10f ",creal(parDeltaN[0]));
  printf("%13.10f ",creal(parDeltaN[1]));
  printf("\n");
}

  for(i=0; i<N*N; i++){
    Guij[i]=0.0;
    Gdij[i]=0.0;
  }
  for(k=0; k<Ne; k++){
    if(array[k].sigma==0){
      for(i=0; i<N; i++){
        tmp = conj(array[k].vec[i]);
        for(j=0; j<N; j++){
          Guij[i*N+j] += tmp * array[k].vec[j];
        }
      }
    }else{
      for(i=0; i<N; i++){
        tmp = conj(array[k].vec[i]);
        for(j=0; j<N; j++){
          Gdij[i*N+j] += tmp * array[k].vec[j];
        }
      }
    }
  }

  parN[0] = (Guij[0*N+0]+Guij[1*N+1])/(1.0*Nfourier);
  parN[1] = (Gdij[0*N+0]+Gdij[1*N+1])/(1.0*Nfourier);
  parDeltaN[0] = (Guij[0*N+1]+Guij[1*N+0])/(1.0*Nfourier);
  parDeltaN[1] = (Gdij[0*N+1]+Gdij[1*N+0])/(1.0*Nfourier);

//  printf("%6d ",step);
  tmp_err = 0.0;
  for(i=0; i<N*N; i++){
    tmp_err += cabs(oldGuij[i] - Guij[i])*cabs(oldGuij[i] - Guij[i]);
    tmp_err += cabs(oldGdij[i] - Gdij[i])*cabs(oldGdij[i] - Gdij[i]);
//    printf("%13.10f+I%13.10f ",creal(oldGuij[i]),cimag(oldGuij[i]));
//    printf("%13.10f+I%13.10f ",creal(oldGdij[i]),cimag(oldGdij[i]));
  }
//  printf("\n");
  if(sqrt(tmp_err)/Nfourier < orderpara_eps
    && fabs(creal(parN[0] - oldparN[0])) < orderpara_eps
    && fabs(creal(parN[1] - oldparN[1])) < orderpara_eps
    && fabs(creal(parDeltaN[0] - oldparDeltaN[0])) < orderpara_eps
    && fabs(creal(parDeltaN[1] - oldparDeltaN[1])) < orderpara_eps
  ){
    break;
  }

/*
  for(i=0; i<N*N; i++){
    oldGuij[i] = Guij[i]*mixing + oldGuij[i]*(1.0-mixing);
    oldGdij[i] = Gdij[i]*mixing + oldGuij[i]*(1.0-mixing);
  }
*/
  for(spin=0; spin<2; spin++){
    parN[spin] = parN[spin]*mixing + oldparN[spin]*(1.0-mixing);
    parDeltaN[spin] = parDeltaN[spin]*mixing + oldparDeltaN[spin]*(1.0-mixing);
  }

  for(i=0; i<N*N; i++){
    oldGuij[i] = Guij[i];
    oldGdij[i] = Gdij[i];
  }
  for(spin=0; spin<2; spin++){
    oldparN[spin] = parN[spin];
    oldparDeltaN[spin] = parDeltaN[spin];
  }

}

  if(step==Niter){
    printf("\n\n");
    printf("# not converged within %d steps !!!\n",Niter);
    printf("# ");
  }else if(
    creal(parN[0]+parDeltaN[0])>1.0+1.0e10 ||
    creal(parN[0]+parDeltaN[0])<0.0-1.0e10 ||
    creal(parN[0]-parDeltaN[0])>1.0+1.0e10 ||
    creal(parN[0]-parDeltaN[0])<0.0-1.0e10 ||
    creal(parN[1]+parDeltaN[1])>1.0+1.0e10 ||
    creal(parN[1]+parDeltaN[1])<0.0-1.0e10 ||
    creal(parN[1]-parDeltaN[1])>1.0+1.0e10 ||
    creal(parN[1]-parDeltaN[1])<0.0-1.0e10
  ){
    printf("\n\n");
    printf("# converged to <n> < 0 or <n> >1 (%d steps) !!! \n",step);
    printf("# ");
  }else{
    printf("\n\n");
    printf("# converged in %d steps\n",step);
  }
  printf("%6d %13.10f ",step,Ene/Nfourier);
  printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
  printf("%13.10f ",creal(parN[0]+parN[1])*0.5);
  printf("%13.10f ",creal(parN[0]));
  printf("%13.10f ",creal(parN[1]));
  printf("%13.10f ",creal(parDeltaN[0]+parDeltaN[1])*0.5);
  printf("%13.10f ",creal(parDeltaN[0]));
  printf("%13.10f ",creal(parDeltaN[1]));
  printf("\n");

  printf("\n\n");
  printf("# SzA SzB NA NB\n");
  printf("%13.10f ",creal(parN[0]+parDeltaN[0]-parN[1]-parDeltaN[1])*0.5);// SzA
  printf("%13.10f ",creal(parN[0]-parDeltaN[0]-parN[1]+parDeltaN[1])*0.5);// SzB
  printf("%13.10f ",creal(parN[0]+parDeltaN[0]+parN[1]+parDeltaN[1]));// NA
  printf("%13.10f ",creal(parN[0]-parDeltaN[0]+parN[1]-parDeltaN[1]));// NB
  printf("\n");

  if(fabs(array[Ne-1].data-array[Ne].data) < 1.0e-16){
    printf("\n\n");
    printf("# !!! OPEN SHELL !!!\n");
  }

if(show_band==1){
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
    intdosu[i]=0;
    intdosd[i]=0;
  }
  for(i=0; i<Norbital; i++){
    j=(int)((array[i].data-Emin)/Edelta);
    intdos[j]++;
    if(array[i].sigma==0){
      j=(int)((array[i].data-Emin)/Edelta);
      intdosu[j]++;
    }else{
      j=(int)((array[i].data-Emin)/Edelta);
      intdosd[j]++;
    }
  }
  printf("\n");
  for(i=0; i<Ndiv; i++){
    dos[i]=intdos[i]/(double)(Norbital*Edelta);
    dosu[i]=intdosu[i]/(double)(Norbital*Edelta);
    dosd[i]=intdosd[i]/(double)(Norbital*Edelta);
    printf("%d %f %d %f %d %f %d %f\n",i,Emin+i*Edelta-Efermi,intdos[i],dos[i],intdosu[i],dosu[i],intdosd[i],dosd[i]);
  }
}

  free(array);
  free(Epsilon[0]);
  free(Epsilon);
  free(EpsilonQ[0]);
  free(EpsilonQ);

  return 0;
}
