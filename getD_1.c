#include <stdio.h>
#include <string.h>
#include <math.h>

/* Get distances between samples and clusters */

void getD_1(double *KmatIn, double *f, int *K, int *Nin, double *win, double *alpha, double *g, double *amount_points, double *connection_strength)
{
    /* declarations */
    int i, j, k;

    /* calculate distances in feature space */
    for(k = 0; k < *K; k++)
    {
        for(i = 0; i < *Nin; i++)
        {
            connection_strength[k + *K * i] = win[i] * f[k + *K * i] * alpha[k + *K * i];
            amount_points[k] += connection_strength[k + *K * i];
        }
    }

    for(k = 0; k < *K; k++)
    {
        double summa = 0;
        for(i = 0; i < *Nin; i++)
        {
            for(j = i; j < *Nin; j++)
            {
                if(i != j)
                    summa += 2.0 * connection_strength[k + *K * i] * connection_strength[k + *K * j] * KmatIn[i + *Nin * j];
                else
                    summa += connection_strength[k + *K * i] * connection_strength[k + *K * j] * KmatIn[i + *Nin * j];
            }
        }
        
        g[k] = summa / pow (amount_points[k], 2.0);
    }
}
