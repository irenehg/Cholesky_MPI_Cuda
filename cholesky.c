#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>


void Descomposicion(double *a,double *l, int n){
  int i, j, k;
  double sum;

  for(i=0; i<n ; i++){
    sum = 0;
    for(j=0; j<=i; j++){
      sum = 0;
      if(i==j){  //elemento de la diagonal principal
        for(k=0; k < i; k++)
          sum += (l[i*n+k]*l[i*n+k]);

        l[i*n+j]= sqrt(a[i*n+j]-sum);
      }
      else if(i>j){
        for(k=0; k < j; k++)
          sum += (l[j*n+k]*l[i*n+k]);

        l[i*n+j] = (a[i*n+j] - sum) / l[j*n+j];
      }
      else{
      	l[i*n+j] = 0.0;
      }
    }
  }
}


int main(int argc, char *argv[]){
  int n, i, j;

  n = 4000; //por defecto tam de la matriz

  if(argc == 2){    //si nos pasan parametro, es el tam de la matriz
    n = atoi(argv[1]);
  }

  //Partimos de una matrz A simetrica definida positiva (A*x=B)
  double *a = malloc(sizeof(double) * n * n);

  //Rellenamos valores de A con floats aleatorios en la triangular inferior
  for(i=0; i<n; i++){
    for(j=0; j<i; j++){
      a[i*n+j]=(rand() %10)/((rand() %9)+1.0);
    }
  }

  //Completamos la triangular superior de forma que sea simetrica
  for(i=0; i<n; i++){
    for(j=i+1; j<n; j++){
      a[i*n+j]=a[j*n+i];
    }
  }

  //Completamos con elementos de la diagonal de forma que nigun 'subdeterminante' sea <=0
  for (i=0; i< n; i++) {
    double s = 0.0;
    for (j=0; j< i; j++)
      s += a[i*n+j];

    for (j= i+1; j< n; j++)
      s += a[j*n+i];

    a[i*n+i] = s *5;
  }



  //Descomponemos A en matriz triangular inferior L y triangular superior U
  double *l = malloc(sizeof(double) * n * n);
  double *u = malloc(sizeof(double) * n * n);


  //inicializamos reloj
	struct timespec cgt1,cgt2;
	double ncgt; //para tiempo de ejecución

	clock_gettime(CLOCK_REALTIME,&cgt1);


  //Descomponemos A en Matriz L*U. Como U = Lt solo calculamos L
  Descomposicion(a,l,n);


  clock_gettime(CLOCK_REALTIME,&cgt2);
  ncgt=(double) (cgt2.tv_sec-cgt1.tv_sec)+ (double) ((cgt2.tv_nsec-cgt1.tv_nsec)/(1.e+9));


  // printf("\n Mostramos ahora la descomposición\n");
  // printf("\n A \n");
  // for(i = 0; i < n; i++){
  //   for(j=0; j<n;j++)
  //     printf("%f ",a[i*n+j]);
  //   printf ( "\n");
  // }
  //
  // printf("\n L \n");
  // for(i = 0; i < n; i++){
  //   for(j=0; j<n;j++)
  //     printf("%f ",l[i*n+j]);
  //   printf ( "\n");
  // }
  //
  // printf("\n U \n");
  // for(i = 0; i < n; i++){
  //   for(j=0; j<n;j++)
  //     printf("%f ",l[j*n+i]);
  //   printf ( "\n");
  // }

  printf("\n Para n=%d , tiempo = %11.9f segundos\n", n, ncgt);

  free(a);
  free(l);

}
