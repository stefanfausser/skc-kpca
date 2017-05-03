#include <stdio.h>
#include <string.h>
#include <math.h>

/* Polynomial Kernel */

void poly_km_flexible(double *data_row, double *data_col, int *D, int *N_row, int *N_col, double *Kmat, double *offset, double *degree, int *normalizeKernel, int *isDiag, int *isSymmetric)
{
    int i, j, d;

    if(*isDiag)
    {
        for(i = 0; i < *N_row; i++)
        {
            if(*normalizeKernel)
                Kmat[i] = 1.0;
            else
            {
                double summa = 0;
            
                for(d = 0; d < *D; d++)
                    summa += pow(data_row[i + *N_row * d], 2.0) / (double) *D;
            
                Kmat[i] = pow(summa + *offset, *degree);
            }
        }
        
        return;
    }
    
    for(i = 0; i < *N_row; i++)
    {
        int startJ = 0;
        
        if(*isSymmetric)
            startJ = i;
        
        for(j = startJ; j < *N_col; j++)
        {
            double summa = 0;
            double summa1 = 0, summa2 = 0;
            
            for(d = 0; d < *D; d++)
            {
                if(*normalizeKernel)
                {
                    summa += data_row[i + *N_row * d] * data_col[j + *N_col * d];

                    summa1 += pow(data_row[i + *N_row * d], 2.0);
                    summa2 += pow(data_col[j + *N_col * d], 2.0);
                }
                else
                    summa += data_row[i + *N_row * d] * data_col[j + *N_col * d] / (double) *D; 
            }
            
            if(*normalizeKernel)
                Kmat[i + *N_row * j] = pow(summa + *offset, *degree) / sqrt(pow(summa1 + *offset, *degree) * pow(summa2 + *offset, *degree));
            else
                Kmat[i + *N_row * j] = pow(summa + *offset, *degree);
            
            if((*isSymmetric) && i != j)
                Kmat[j + *N_row * i] = Kmat[i + *N_row * j];
        }
    }
}
