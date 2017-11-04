/***************************
* Matrix vs Vector Toolbox *
* By: Leo Nugraha          *	
* July 02, 2017	           *	
***************************/
void mat_dot_vec(double *sln, double **m, double *x, int size)
{
    int i, j;	double sum;
    for(i=0; i<size; i++)
    {
        sum=0.0;
        for(j=0; j<size; j++){
            sum = sum + m[i][j]*x[j];        
		}
        sln[i] = sum;
    }
    return;
}

void vec_add_scaled_vec(double *sln, double *a, double c, double *b, int size)
{
    int i;
    for(i=0; i<size; i++)
        sln[i] = a[i]+c*b[i];
    return;    
}

// Handle dot products between 2 vectors
// out_vector = dot_product(vector_01, vector_02, length)
double dot_product(double *a, double *b, int size)
{
    int i; double sum = 0.0;
    for(i=0; i<size; i++){
        sum = sum+a[i]*b[i];
    }
    return sum;
}
