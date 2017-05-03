#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

/* Get distances between samples and clusters */

void getD_2(double *KmatOut, double *KmatDiag, int *K, int *Nin, int *Nout, double *g, double *amount_points, double *connection_strength, double *d, double *err, double *summed_diameter)
{
    /* declarations */
    int i, j, k;

    /* calculate distances in feature space */

    *err = 0;
    
    for(i = 0; i < *Nout; i++)
    {
        double minval = INFINITY;
        int kstar = 0;
        bool first = true;
        
        for(k = 0; k < *K; k++)
        {
            double summa = 0;
            for(j = 0; j < *Nin; j++)
                summa += connection_strength[k + *K * j] * KmatOut[i + *Nout * j];
            
            double h = -2.0 * summa / amount_points[k];

            d[k + *K * i] = KmatDiag[i] + h + g[k];

            double val = d[k + *K * i];
            if (val < minval || first)
            {
                first = false;
                kstar = k;
                minval = val;
            }            
        }
        
        *err += minval;

        summed_diameter[kstar] += minval;
    }
}
