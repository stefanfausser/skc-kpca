#include <stdio.h>
#include <string.h>
#include <math.h>

/* Get distances between prototypes */

void kernel_dist_prototypes3(double *din, int *A, double *dout, int *K, int *N)
{
    /* declarations */
    unsigned k, k2, i;

    double amount_points[*K];

    for(k = 0; k < *K; k++)
    {
        amount_points[k] = 0;

        for(i = 0; i < *N; i++)
            amount_points[k] += A[k + *K * i];

        for(k2 = 0; k2 < *K; k2++)
            dout[k + *K * k2] = 0;
    }

    for(k = 0; k < *K; k++)
    {
        for(k2 = k + 1; k2 < *K; k2++)
        {
            double sum = 0;

            for(i = 0; i < *N; i++)
            {
                if(A[k + *K * i])
                    sum += din[k2 + *K * i];
                else if(A[k2 + *K * i])
                    sum += din[k + *K * i];                
            } /* for i */

            dout[k + *K * k2] = sum / (amount_points[k] + amount_points[k2]);
            
            // Symmetric matrix
            dout[k2 + *K * k] = dout[k + *K * k2];
        } /* for k2 */
    } /* for k */
}
