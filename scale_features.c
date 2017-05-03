#include <stdio.h>
#include <string.h>

void scale_features(double *data, int *N, int *D, double *offset, double *divisor)
{
    int i, d;

    for(i = 0; i < *N; i++)
    {
        for(d = 0; d < *D; d++)
            data[i + *N * d] = (data[i + *N * d] + offset[d]) / divisor[d];
    }
}
