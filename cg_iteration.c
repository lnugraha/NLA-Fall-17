/******************************************************************
* CG Iteration Method - Basic Mode				  *
* Input: A Symmetric, Positive, Definite SQUARE matrix		  *	
* Outputs: Number of Iterations to achieve convergence & Errors   *
* (Note) You can also print the solution			  *
* By: Leo Nugraha						  *
* Date: July 02, 2017						  *
******************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "matrix_vector_toolbox.c"

clock_t start,end;     		// CPU Timing - from time.h 
double cpu_time_used;		// Compute the time difference

void copy_vector(double *, double *, int);
void cg_method(int, double **, double *, double, int, int *, double *);

int main()
{
    int n,i,size; double value;
    FILE *inputfile;
    char in_file[100]; 	

    printf("Matrix Name: ");	// From .mtx to .txt
    scanf("%s", in_file);	// Matlab code from Matrix Market

	n = 1;
    inputfile = fopen(in_file, "r");			// Determine the size of your matrix
    while(fscanf(inputfile, "%lf ", &value) != EOF){	// Since the input matrix is a square matrix
		++n;					// Taking a square-root MUST yield an integer number
     }
	
    size = sqrt(n);
    printf("Matrix Size: %d \n", size);
    fclose(inputfile);

    // Define the solution b = [1 1 1 ... 1]'
    double *b;
    b = (double*) calloc(size,sizeof(double));
    for (i = 0; i < size; i++) b[i] = 1.0;
    
    // Define the input A = (size X size), where A is an SPD matrix
    // Note: To ensure the convergency, A must be a diagonal dominant matrix (the matrix diagonal has the highest value)
    double **A;
    A = (double**) calloc(size,sizeof(double*));
    for(i = 0; i < size; i++)
	A[i] = (double*) calloc(size,sizeof(double));

    inputfile = fopen(in_file, "r");
    for(int row = 0; row < size; row++){
		for(int col = 0; col < size; col++){
	    	fscanf(inputfile, "%lf ", &A[row][col]);
        }
    }

    // Define the output x = [... ... ... ...]', which will be determined through iterations
    // Our initial guess is all zeroes, and will be updated for each iteration
    double *x;
    x = (double*) calloc(size,sizeof(double)); 
    for(i = 0; i < size; i++) x[i] = 0.0; 

    double eps   = 1.0e-6;	// Error Tolerance
    int maxiter  = size*size;	// Maximum Number of Iteration
    int cnt      = 0;		// Number of iteration used
	
    start = clock();
    cg_method(size, A, b, eps, maxiter, &cnt, x);
    end = clock();

    printf("Computed %d Iterations\n",cnt);
    /* Compute the Error */
    double sum = 0.0; 
    for(i=0; i<size; i++){
	double d = x[i] - 1.0;
	sum += (d >= 0.0) ? d : -d;
    }
    printf("Error : %.3e\n",sum);

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf( "CPU Time: %16.8E seconds \n",cpu_time_used);

    // See your solution
    // for (i = 0; i < size; i++)
    //	printf("%.6lf \n", x[i]);
	  
   fclose(inputfile);
   return 0;
}

void cg_method(int size, double **A, double *b, double epsilon, int maxiter, int *numit, double *x)
{
	double alpha, beta, rdr, new_rdr, dMd;
	double	*r_ptr = (double*) calloc(size,sizeof(double));	//residual vector (r)
	double	*d_ptr = (double*) calloc(size,sizeof(double));	//direction vector (p)
	double	*e_ptr = (double*) calloc(size,sizeof(double));	//vector by matrix sln (Ax_k)
///////////////////////////////////////////////////
	mat_dot_vec(e_ptr, A, x, size);							// e_ptr = A*x, where A is the input matrix and x i sour guessed vector
	vec_add_scaled_vec(d_ptr, b, -1.0, e_ptr, size);		// d_ptr = b -1.0*e_pt ->  
//////////////////////////////////////////////////
	copy_vector(r_ptr, d_ptr, size);							// First Iter: p_0 = r_0
	rdr = dot_product(r_ptr, r_ptr, size);					// Find r_0T r_0
	
	while(rdr > epsilon)
	{
		mat_dot_vec(e_ptr, A, d_ptr, size);

		dMd = dot_product(d_ptr, e_ptr, size);			// dMd = p_kT * A * p_k
		alpha = rdr / dMd;							// alpha = (r_kT * r_k)/(p_kT * A * p_k)

		vec_add_scaled_vec(x, x, alpha, d_ptr, size);
		vec_add_scaled_vec(r_ptr, r_ptr, -alpha, e_ptr, size);

		new_rdr = dot_product(r_ptr, r_ptr, size);
		beta =  new_rdr / rdr;

		vec_add_scaled_vec(d_ptr, r_ptr, beta, d_ptr, size);
		if((*numit%10)==0){
	    	printf("Iteration,  pAp, rr, tt : %d  %16.8E  %16.8E  %16.8E \n", *numit, dMd, rdr, new_rdr); 
    	}
		rdr = new_rdr;
		*numit = *numit + 1;
		int k = *numit;
		if(k > maxiter)	
		{
			printf("Fail to Converge! \n");		
			exit(1);
		}
		
	}	// end while-loop
		free(r_ptr); free(d_ptr); free(e_ptr);
}

void copy_vector(double *copy, double *orig, int size)
{
	int i;
	for(i = 0; i < size; ++i)
		copy[i] = orig[i];
}
