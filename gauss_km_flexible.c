#include <stdio.h>
#include <string.h>
#include <math.h>

/* Gauss Kernel */

void gauss_km_flexible(double *data_row, double *data_col, int *D, int *N_row, int *N_col, double *Kmat, double *sigma, int *isSymmetric)
{
    int i, j, d;

    for(i = 0; i < *N_row; i++)
    {
        int startJ = 0;
        
        if(*isSymmetric)
            startJ = i;
        
        for(j = startJ; j < *N_col; j++)
        {
            double euclidean_dist = 0;
            for(d = 0; d < *D; d++)
                euclidean_dist += pow(data_row[i + *N_row * d] - data_col[j + *N_col * d], 2.0);
            
            Kmat[i + *N_row * j] = exp(- euclidean_dist / (2.0 * pow(*sigma,2.0)));

            if((*isSymmetric) && i != j)
                Kmat[j + *N_row * i] = Kmat[i + *N_row * j];
        }
    }
}
