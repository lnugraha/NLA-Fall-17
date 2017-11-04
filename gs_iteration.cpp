#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void gs_method (int n, double **A, double *b, double epsilon, int maxit, int *numit, double *x );
int main ( )
{
    int n,i, size; double value;
    FILE *inputfile;
	char in_file[100]; 	
	printf("Enter the Matrix Name: ");
	scanf("%s", in_file);
	inputfile = fopen(in_file, "r");
	while(fscanf(inputfile, "%lf ", &value) != EOF){
		++n;
	}
	
	size = sqrt(n);
	printf("Matrix Size: %d \n", size);
	fclose(inputfile);
	// Define the solution b = [1 1 1 ... 1]'
	double *b;
	b = (double*) calloc(size,sizeof(double));
	for (i = 0; i < size; i++){
		b[i] = 1.0;
	}

	// Define the input A = (size X size), where A is an SPD matrix
	// Note: To ensure the convergency, A must be a diagonal dominant matrix (the matrix diagonal has the highest value)
	double **A;
   	A = (double**) calloc(size,sizeof(double*));
   	for(i = 0; i < size; i++)
		A[i] = (double*) calloc(size,sizeof(double));

	inputfile = fopen(in_file, "r");
	for (int row = 0; row < size; row++){
		for (int col = 0; col < size; col++){
			fscanf(inputfile, "%lf ", &A[row][col]);
		}
	}

	// Define the output x = [... ... ... ...]', which will be determined through iterations
	// Our initial guess is all zeroes, and will be updated for each iteration
	double *x;
	x = (double*) calloc(size,sizeof(double)); 
    for(i = 0; i < size; i++){
    	x[i] = 0.0;
	} 

    double eps   = 1.0e-6;	// Error Tolerance
    int maxiter  = 5000;	// Maximum Number of Iteration
	int cnt      = 0;		// Number of iteration used
	
    gs_method(size,A,b,eps,maxiter,&cnt,x);
    printf("Computed %d Iterations\n",cnt);
	
	/* Compute the Error */
	double sum = 0.0; 
	for(i = 0; i < size; i++){
		double d = x[i] - 1.0;
		sum += (d >= 0.0) ? d : -d;
    }
	printf("error : %.3e\n",sum);
	// See your solution
	/* for (i = 0; i < size; i++)
		printf("%.6lf \n", x[i]);
	*/
	fclose(inputfile);
	return 0;
}

void gs_method (int n, double **A, double *b, double epsilon, int maxiter, int *numit, double *x){
   double *dx = (double*) calloc(n,sizeof(double));
   int i,j,k;
   for(k=0; k < maxiter; k++){
      double sum = 0.0;
      for(i=0; i<n; i++){
         dx[i] = b[i];
         for(j=0; j<n; j++)
            dx[i] -= A[i][j]*x[j]; 
            
         dx[i] /= A[i][i]; x[i] += dx[i];
         sum += ( (dx[i] >= 0.0) ? dx[i] : -dx[i]);
      }
      if(sum <= epsilon) break;
   }
	if (k >= maxiter){
		printf("Fail to Converge! \n");		
		exit(1);
	}
   *numit = k+1; free(dx);
}
