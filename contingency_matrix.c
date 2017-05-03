#include <stdio.h>
#include <string.h>
#include <math.h>

typedef enum {false, true} bool;

void contingency_matrix(double *f, int *sample_labels, int *K, int *N, int *M, int *c)
{
    int i, k, kstar;
    
    // Update the contingency matrix

    for(i = 0; i < *N; i++)
    {
        // Get cluster with index kstar with the highest sample-to-cluster assignment value
        kstar = 0;
        double maxval = f[0 + *K * i];
        for(k = 1; k < *K; k++ )
        {
            if(f[k + *K * i] > maxval)
            {
                maxval = f[k + *K * i];
                kstar = k;
            }
        }
        
        int m = sample_labels[i] - 1;
        if(m < 0)
        {
            printf("### sample_labels[i = %i] < 0\n", i);
            return;
        }

        c[m + *M * kstar]++;
    }
}

