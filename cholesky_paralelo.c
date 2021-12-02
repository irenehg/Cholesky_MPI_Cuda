#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>


int main(int argc, char *argv[]){
  int n, i, j, k, p, c, r;
  double sum,aux;


  n = 4; //por defecto tam de la matriz

  if(argc != 1){    //si nos pasan parametro, es el tam de la matriz
    n = atoi(argv[1]);
  }

  //Partimos de una matrz A simetrica definida positiva (A*x=B)
  double *a = malloc(sizeof(double) * n * n);

  // //Rellenamos valores de A con floats aleatorios en la triangular inferior
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

  // a[0] = 4;
  // a[1] = -1;
  // a[2] = 0;
  // a[3] = 2;
  // a[4] = -1;
  // a[5] = 4;
  // a[6] = -1;
  // a[7] = 0;
  // a[8] = 0;
  // a[9] = -1;
  // a[10] = 4;
  // a[11] = 1;
  // a[12] = 2;
  // a[13] = 0;
  // a[14] = 1;
  // a[15] = 3;



  //Descomponemos A en matriz triangular inferior L y triangular superior U

  double *l = malloc(sizeof(double) * n * n);
  for(i=0; i < n; i++){
    for(j = 0; j < n; j++){
      if(j<=i)
        l[i*n+j]=a[i*n+j];
      else
        l[i*n+j] = 0.0;
    }
  }

  // Inicializa MPI

  int num_procesos, id_proceso;
	MPI_Init(&argc, &argv);

	// Inicialización del identificador de proceso y del número de procesos
	MPI_Comm_size(MPI_COMM_WORLD, &num_procesos);

	MPI_Comm_rank(MPI_COMM_WORLD, &id_proceso);

  //inicializamos reloj
	struct timespec cgt1,cgt2;
	double ncgt; //para tiempo de ejecución


  clock_gettime(CLOCK_REALTIME,&cgt1);

  for(k=0; k<n ; k++){ //columna
    if( id_proceso == 0){    //este proceso calcula todos elementos de la diagonal

      sum = 0;

      //elemento de la diagonal principal
      for(j=0; j < k; j++)
        l[k*n+k] -= (l[k*n+j]*l[k*n+j]); //depende de los anteriores de la fila

      l[k*n+k]= sqrt(l[k*n+k]);  //calculado :)

      //Ahora le manda a cada proceso que calcule los elementos que no estan en la diagonal
      for(p = 1; p < n-k; p++){   //este programa tiene tantos procesos como filas la matriz
        for(c = 0; c <= k; c++){
          for(r=c; r < n; r++){
            MPI_Send(&l[r*n+c], 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
          }
        }
        MPI_Recv(&l[(p+k)*n+k],1,MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

    }

    else if(id_proceso < n-k){ //si no soy el 0 y mi id entra dentro de las filas a calcular
      for(c = 0; c <= k; c++){
        for(r = c; r < n; r++){
          MPI_Recv(&l[r*n+c],1,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //recibo valores de L
        }
      }

      for(i=0; i < k; i++)
        l[(k+id_proceso)*n+k] -= l[k*n+i]*l[(k+id_proceso)*n+i];

      l[(k+id_proceso)*n+k] /= l[k*n+k];          // Calculo :)
      MPI_Send(&l[(k+id_proceso)*n+k], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // y mando!
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);



  if(id_proceso == 0){
    clock_gettime(CLOCK_REALTIME,&cgt2);
    ncgt=(double) (cgt2.tv_sec-cgt1.tv_sec)+ (double) ((cgt2.tv_nsec-cgt1.tv_nsec)/(1.e+9));

    printf("\n Para n=%d , tiempo = %11.9f segundos\n", n, ncgt);

    // printf("\n L\n");
    // for(i = 0; i < n; i++){
    //   for(j = 0; j < n; j++){
    //     printf("%f ", l[i*n+j]);
    //   }
    //   printf ( "\n");
    // }
    //
    //
    // printf("\n U\n");
    //
    // for(i = 0; i < n; i++){
    //   for(j = 0; j < n; j++){
    //     printf("%f ", l[j*n+i]);
    //   }
    //   printf ( "\n");
    // }
  }


  MPI_Finalize();


  free(a);
  free(l);

}
