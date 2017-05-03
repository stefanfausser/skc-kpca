#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define SPEEDUP_1
#define SPEEDUP_2
#define SPEEDUP_3

// Real O(n + c) algorithm, where n is the number of samples for the approximation

void greedy_approximation(double *Kmat, double *f, int *A, double *w, double *alpha, int *K, int *nSamplesMax, int *N, int *best_samples, int *uniqueSampleIndices, double *best_samples_dist, int *number_samples, int *limitSamplesByHardAssignments, int *localMinimum, int *Nstart, int *ret)
{
    /* declarations */
    double g[*K];
    double amount_points[*K];
    double connection_strength[*K][*N];
    int i, j, k, ka, i_test, i_tmp, z;
    int sampleRepos[*K][*N];
    int sampleReposLen[*K];
    int number_samples_curr[*K];   
#ifdef SPEEDUP_3
    int sampleInCluster[*K][*N];
    int sampleInClusterLen[*K];
#endif

    *ret = 0;
    
    for(k = 0; k < *K; k++)
    {
        amount_points[k] = 0;
        sampleReposLen[k] = 0;
        number_samples[k] = 0;
        number_samples_curr[k] = 0;
        sampleInClusterLen[k] = 0;
        
        for(i = 0; i < *N; i++)
        {
            connection_strength[k][i] = w[i] * f[k + *K * i] * alpha[k + *K * i];
            amount_points[k] += connection_strength[k][i];

#ifdef SPEEDUP_3
            if(connection_strength[k][i])
            {
                sampleInCluster[k][sampleInClusterLen[k]] = i;
                sampleInClusterLen[k]++;
            }
#endif

            if(i < *Nstart)
                continue;

            // Add a unique sample only once to the repo:
            // If sample with index i is not unique, i.e. it is duplicated then ignore it and continue
            if(!uniqueSampleIndices[i])
                continue;

            if(*limitSamplesByHardAssignments)
            {
                // Ignore a sample not assigned to the cluster                
                if(!A[k + *K * i])
                    continue;
            }
            
            sampleRepos[k][sampleReposLen[k]] = i;
            sampleReposLen[k]++;
        }
    }

    for(k = 0; k < *K; k++)
    {
        double summa = 0;
        for(i = 0; i < *N; i++)
        {
            for(j = i; j < *N; j++)
            {
                if(i != j)
                    summa += 2.0 * connection_strength[k][i] * connection_strength[k][j] * Kmat[i + *N * j];
                else
                    summa += connection_strength[k][i] * connection_strength[k][j] * Kmat[i + *N * j];
            }
        }
        
        g[k] = summa / pow (amount_points[k], 2.0);
    }
    
    /* calculate distances in feature space */

    double err[*K];
    bool selectable[*K];
    double summa1g[*K];
    double summa2g[*K];
    
    for(k = 0; k < *K; k++)
    {
        err[k] = INFINITY;
        selectable[k] = true;
        summa1g[k] = 0;
        summa2g[k] = 0;
    }

    // Determine up to nSamplesMax winner-samples

    int t = 0;
    
    while(1)
    {
        // Determine the cluster with the highest error
        
        double maxErr = -INFINITY;
           
        int bestK = 0;
        
        for(k = 0; k < *K; k++)
        {
            if(err[k] > maxErr && selectable[k])
            {
                bestK = k;
                maxErr = err[k];
            }
        }
    
        // Choose a sample for cluster bestK
        
        k = bestK;
        
        ka = number_samples_curr[k];

        number_samples_curr[k]++;        
            
        double mindist = INFINITY;

        int bestz = 0;

        double bestsumma1 = 0, bestsumma2 = 0;
        
        for(z = 0; z < sampleReposLen[k]; z++)
        {
            i_test = sampleRepos[k][z];

            double amount_points2 = ka + 1;

#ifdef SPEEDUP_1
            double summa1 = summa1g[k];
            
            for(i_tmp = 0; i_tmp <= ka; i_tmp++)
            {
                if(i_tmp == ka)
                    i = i_test;
                else
                    i = best_samples[k + *K * i_tmp];
                
                if(i != i_test)
                    summa1 += 2.0 * Kmat[i + *N * i_test];
                else
                    summa1 += Kmat[i_test + *N * i_test];
            }

#else
            double summa1 = 0;
            for(i_tmp = 0; i_tmp <= ka; i_tmp++)
            {
                int j_tmp;
                
                for(j_tmp = i_tmp; j_tmp <= ka; j_tmp++)
                {
                    if(i_tmp == ka)
                        i = i_test;
                    else
                        i = best_samples[k + *K * i_tmp];

                    if(j_tmp == ka)
                        j = i_test;
                    else
                        j = best_samples[k + *K * j_tmp];

                    if(i != j)
                        summa1 += 2.0 * Kmat[i + *N * j];
                    else
                        summa1 += Kmat[i + *N * j];
                }
            }
#endif
                            
            double g2 = summa1 / pow (amount_points2, 2.0);

#ifdef SPEEDUP_2
            double summa2 = summa2g[k];

    #ifdef SPEEDUP_3
            for(i_tmp = 0; i_tmp < sampleInClusterLen[k]; i_tmp++)
            {
                i = sampleInCluster[k][i_tmp];
                
                summa2 += connection_strength[k][i] * Kmat[i + *N * i_test];
            }
    #else
            for(i = 0; i < *N; i++)
                summa2 += connection_strength[k][i] * Kmat[i + *N * i_test];
    #endif /* SPEEDUP_3 */
#else
            double summa2 = 0;
            for(i = 0; i < *N; i++)
            {
                int j_tmp;
                
                for(j_tmp = 0; j_tmp <= ka; j_tmp++)
                {
                    if(j_tmp == ka)
                        j = i_test;
                    else
                        j = best_samples[k + *K * j_tmp];

                    summa2 += connection_strength[k][i] * Kmat[i + *N * j];
                }
            }
#endif
            double h2 = -2.0 * summa2 / (amount_points2 * amount_points[k]);

            double d = g2 + h2 + g[k];
            if(d < mindist)
            {
                mindist = d;
                best_samples[k + *K * ka] = i_test;
                best_samples_dist[k + *K * ka] = d;
                bestz = z;
                bestsumma1 = summa1;
                bestsumma2 = summa2;
            }
        } /* for z */

        summa1g[k] = bestsumma1;
        summa2g[k] = bestsumma2;

        int k_tmp, nSamples = 0;
        for(k_tmp = 0; k_tmp < *K; k_tmp++)
            nSamples += number_samples[k_tmp];
        
        if(mindist < err[k])
        {
            if((number_samples_curr[k] - number_samples[k]) <= (*nSamplesMax - nSamples))
            {
                err[k] = mindist;
                number_samples[k] = number_samples_curr[k];
            }
            else
                selectable[k] = false;
        }
        else if(*localMinimum)
            selectable[k] = false;

        if(number_samples_curr[k] >= *nSamplesMax)
            selectable[k] = false;
        
        // Remove the sample from the repos for this cluster (prevents to re-select it)
        if(sampleReposLen[k] >= 2)
        {
            sampleReposLen[k]--;
            sampleRepos[k][bestz] = sampleRepos[k][sampleReposLen[k]];
        }
        else
        {
            // Should not happen, because of 'if(number_samples_curr[k] >= *nSamplesMax)'
            selectable[k] = false;
        }
        
        int nSelectable = 0;
        for(k_tmp = 0; k_tmp < *K; k_tmp++)
            nSelectable += selectable[k_tmp];
                
        if(!nSelectable)
            break; /* while */
            
        t++;
        
        if(t > 100000)
        {
            *ret = -1;
            return;
        }            
    } /* while(1) */
}
