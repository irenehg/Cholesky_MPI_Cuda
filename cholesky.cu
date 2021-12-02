#include <stdio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_functions.h>

using namespace std;
using std::vector;


//Kernel encargado de los elementos de la diagonal
__global__ void Diagonal(float* matriz_a, float* matriz_l, int n, int numCol) {

    int fila = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    // Solo sacamos elementos de la diagonal principal
    if (fila < n && col < n) {
        if (fila == col && col == numCol) {
            float diagSum = 0.0f;
            for (int k = 0; k < col; k++) {
                diagSum += (matriz_l[col * n + k]) * (matriz_l[col * n + k]);
            }
            matriz_l[fila * n + col] = sqrt((float)(matriz_a[fila * n + col] - diagSum));
        }
    }
}

//kernel encargado de los elementos de la columna
__global__ void Columna(float* matriz_a, float* matriz_l, int n, int num) {

    int fila = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    // Elementos de la columna por debajo de la diagonal principal
    if (col == num && fila > col && fila < n) {
        float sum = 0.0f;
        for (int k = 0; k < col; k++) {
            sum += (matriz_l[fila * n + k] * matriz_l[col * n + k]);
        }
        matriz_l[fila * n + col] = (matriz_a[fila * n + col] - sum) / matriz_l[col * n + col];
    }
}



int main(int argc, char *argv[]){
  int n, i, j, k, p, c, r;
  double sum,aux;


  n = 4; //por defecto tam de la matriz

  if(argc != 1){    //si nos pasan parametro, es el tam de la matriz
    n = atoi(argv[1]);
  }

  size_t bytes = n * n * sizeof(int);


  //Partimos de una matrz A simetrica definida positiva (A*x=B)
  vector<float> a(n * n);  //no he usado malloc como en la version de mpi porque me daba error

  //Rellenamos valores de A con floats aleatorios en la triangular inferior
  for(i=0; i<n; i++){
    for(j=0; j<i; j++){
      a[i*n+j]=(rand() %9)+1.0;
    }
  }

  //Completamos la triangular superior de forma que sea simetrica
  for(i=0; i<n; i++){
    for(j=i+1; j<n; j++){
      a[i*n+j]=a[j*n+i];
    }
  }

  //Completamos con elementos de la diagonal de forma que ningun 'subdeterminante' sea <=0
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



  // cout << "\nMostramos matriz origen:" << endl;
  // printf("\n A \n");
  // for(i = 0; i < n; i++){
  //   for(j=0; j<n;j++)
  //     printf("%f ",a[i*n+j]);
  //   printf ( "\n");
  // }


  //Descomponemos A en matriz triangular inferior L y triangular superior U
  vector<float> l(n * n, 0.0f);
  for(i=0; i < n; i++){
    for(j = 0; j < n; j++){
      if(j<=i)
        l[i*n+j]=a[i*n+j];
      else
        l[i*n+j] = 0.0;
    }
  }


  // Reservamos espacio en la GPU para ambas matrices
	float* d_A, * d_L;
  cudaMalloc(&d_A, bytes);
  cudaMalloc(&d_L, bytes);


  // Rellenamos las matrices creadas en GPU con valores
  cudaMemcpy(d_A, a.data(), bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_L, l.data(), bytes, cudaMemcpyHostToDevice);


  int threads = 64;
  int blocks = 240;
  dim3 threadsPerBlock(threads, threads);
  int numBlocks = 1;

  //inicializamos reloj
	struct timespec cgt1,cgt2;
	double ncgt; //para tiempo de ejecución

  clock_gettime(CLOCK_REALTIME,&cgt1);


  for (int i = 0; i < n; i++) {
      // Kernel que calcula elemento de la diagonal
      Diagonal << <numBlocks, threadsPerBlock >> > (d_A, d_L, n, i);
      // Kernel que calcula los elementos de esa columna
      Columna << <numBlocks, threadsPerBlock >> > (d_A, d_L, n, i);
  }

  clock_gettime(CLOCK_REALTIME,&cgt2);
  ncgt=(double) (cgt2.tv_sec-cgt1.tv_sec)+ (double) ((cgt2.tv_nsec-cgt1.tv_nsec)/(1.e+9));

  cudaMemcpy(l.data(), d_L, bytes, cudaMemcpyDeviceToHost);

  // cout << "Cholesky en GPU:" << endl;
  // cout << "\nMostramos ahora la descomposición:" << endl;
  // printf("\n L \n");
  // for(i = 0; i < n; i++){
  //   for(j=0; j<n;j++)
  //     printf("%f ",l[i*n+j]);
  //   printf ( "\n");
  // }
  cout << "\nTamaño: " << n << "\tTiempo GPU: " << ncgt << "s\n";
  cout << "\nBloques: "<< blocks << "\thebras por bloque: " << threads << "\n";
  cout << endl;

  // Liberamos memoria
  cudaFree(d_A);
  cudaFree(d_L);

}
