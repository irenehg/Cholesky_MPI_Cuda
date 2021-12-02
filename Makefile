all: cholcu cholp chol

cholcu: cholesky.cu
	nvcc -w cholesky.cu -Iinclude -o cholcu

cholp: cholesky_paralelo.c
	mpicc cholesky_paralelo.c -o cholp -lm -O2

chol:  cholesky.c
	gcc cholesky.c -o chol -lm -O2
