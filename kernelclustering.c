#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdbool.h>

/* Kernel Clustering */

#define MIN_CLASS_FACTOR                1.0
#define MAX_ITERATIONS                  1000
#define EPSILON                         0.000001

static long randValLong(long min, long max);
static void getHardAssignments(double *f, double *f_old, int N, int K, double *d);
static int updateAlpha(double labelledUnlabelledFactor, double nUnlabelledSamples, double nLabelledSamples, double *nLabelledSamplesClass, int N, int Nstart, int K, int M, double *f, double *alpha, int *l, int *winner_label, double *winner_count);
static void getContingency(int N, int K, double *f, double *cMat, int *l, int Nstart);

void kernelclustering(int *seed, int *retval, int *method, double *Kmat, double *w, double *alpha, double *cMat, int *l, double *nLabelledSamplesClass, int *M, int *bSupervised, int *winner_label, double *winner_count, double *f, int *K, int *mu, int *T, double *lambda_start, double *lambda_end, int *N, double *d, double *errSoft, int *verbose, int *decreaseLambda, int *distToPseudoCentre, int *nBestSamples, double *labelledUnlabelledFactor, double *fuzzifier, double *epsilon, int *emptyCluster, int *getHardClusteringsKFCM, int *normalizeAssignmentsKFCM, int *hardAssignmentsApproxPseudoCentresRNG)
{
    /* declarations */
    double amount_points[*K];
    double connection_strength[*K][*N];
    double f_old[*N * *K];
    double g2[*mu];
    
    int i, j, t, k, k2;
    double lambda = *lambda_start, summa;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srandom(*seed);

    double nLabelledSamples = 0;
    double nUnlabelledSamples = 0;
    
    int Nstart = 0;

    if(*distToPseudoCentre)
    {
        int nSummedBestSamples = 0;
        
        // Calculate g2[k] once! (independent on changes in cluster assignments)
        
        int startI = 0;
        
        for(k = 0; k < *mu; k++)
        {
            int endI = startI + nBestSamples[k];

            summa = 0;
            for(i = startI; i < endI; i++)
            {
                for(j = i; j < endI; j++)
                {
                    if(i != j)
                        summa += 2.0 * Kmat[i + *N * j];
                    else
                        summa += Kmat[i + *N * j];
                }
            }

            g2[k] = summa / pow ((double) nBestSamples[k], 2.0);
                                
            startI = endI;

            nSummedBestSamples += nBestSamples[k];
        }
        
        Nstart = nSummedBestSamples;
    }

    *emptyCluster = 0;
    *retval = 0;
    
    int converged = 0;

    int checkEmptyClusters = 0;
    
    if(*verbose)
        printf("Nstart: (%i)\n", Nstart);
    
    if(!(*bSupervised))
    {
        for(k = 0; k < *K; k++)
        {
            for(i = 0; i < *N; i++)
                alpha[k + *K * i] = 1.0;
        }
    }
    else
    {
        // j = 0 => number of unlabelled samples
        nUnlabelledSamples = nLabelledSamplesClass[0];
        nLabelledSamples = 0;
        for(j = 1; j < *M; j++)
            nLabelledSamples += nLabelledSamplesClass[j];
    }

    if(*verbose)
        printf("(kernelclustering) method (%i), lambda (%lf), fuzzifier (%lf), epsilon (%lf)\n", *method, lambda, *fuzzifier, *epsilon);

    t = 0;
        
    while(1)
    {
        t++;
        
        if(*method == 0)
        {
            if(*decreaseLambda)
                lambda = *lambda_start * pow(*lambda_end / *lambda_start, (double) t / *T);

            if(*verbose)
                printf("(kernelclustering) lambda (%lf)\n", lambda);
        }
        
        if(*verbose)
            printf("(kernelclustering) t = %i\n", t);

        if(*verbose)
            printf("(kernelclustering) calculate alpha\n");
        
        // updating / setting alpha is part of the Expection step (E-step)
        // It is NOT part of the Maximization step (M-step), because alpha is not set as to minimize the distances.
        // Instead setting alpha may increase the distances and, thus, changes the error function.
            
        if(*bSupervised)
        {
            int ret = updateAlpha(*labelledUnlabelledFactor, nUnlabelledSamples, nLabelledSamples, nLabelledSamplesClass, *N, Nstart, *K, *M, f, alpha, l, winner_label, winner_count);

            if(ret)
            {
                *retval = ret;
                
                return;
            }
        }
        
        for(k = 0; k < *K; k++)
        {
            amount_points[k] = 0;

            for(i = 0; i < *N; i++)
            {
                connection_strength[k][i] = w[i] * f[k + *K * i] * alpha[k + *K * i];
                amount_points[k] += connection_strength[k][i];
            }
            
            if(*method == 2 || checkEmptyClusters)
            {
                if (amount_points[k] == 0)
                {
                    printf( "### fatal error: cluster (%i) has no assigned points !\n",k );
                    *emptyCluster = 1;
                }
            }
        }

        if(*emptyCluster)
            return;

        if(*verbose)
            printf("t (%i), Calculate distances\n", t);
        
        *errSoft = 0;
            
        /* calculate distances in feature space */
        for(k = 0; k < *K; k++)
        {
            summa = 0;
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

            double g = summa / pow (amount_points[k], 2.0);

            if(*distToPseudoCentre)
            {
                // Calculate the distance from (approximate) pseudo-centres (now treated as weighted sample) to pseudo-centre k
                
                int startI = 0;
                for(k2 = 0; k2 < *mu; k2++)
                {
                    int endI = startI + nBestSamples[k2];

                    summa = 0;
                    for(i = 0; i < *N; i++)
                    {
                        for(j = startI; j < endI; j++)
                            summa += connection_strength[k][i] * Kmat[i + *N * j];
                    }

                    double h2 = -2.0 * summa / ((double) nBestSamples[k2] * amount_points[k]);

                    // Same distance for all samples belonging to the same (approximate) pseudo-centre
                    for(j = startI; j < endI; j++)
                        d[k + *K * j] = g2[k2] + h2 + g;
                    
                    *errSoft += f[k + *K * startI] * d[k + *K * startI];
                    
                    startI = endI;
                } /* for k2 */
            }
            
            for(i = Nstart; i < *N; i++)
            {
                summa = 0;
                for(j = 0; j < *N; j++)
                    summa += connection_strength[k][j] * Kmat[i + *N * j];

                double h = -2.0 * summa / amount_points[k];

                d[k + *K * i] = Kmat[i + *N * i] + h + g;
                
                *errSoft += f[k + *K * i] * d[k + *K * i];                
            }
        } /* for k */
       
        /* output current error */
        
        if(*verbose)
            printf("kernelclustering: (soft) err = %lf\n",*errSoft);
        
        if(converged)
            break; /* Leave main function */
            
        if(*verbose)
            printf("t (%i), calculate assignments\n", t);
            
        if(*method == 0)
        {
            int Nstart2 = 0;            
            
            if((*hardAssignmentsApproxPseudoCentresRNG) && Nstart > 0)
            {
                getHardAssignments(f, NULL, Nstart, *K, d);

                Nstart2 = Nstart;
            }
            
            /* recalculate ranks */

            for(i = Nstart2; i < *N; i++)
            {
                int rank[*K];
                
                /* ranking goes from 0 to K-1, where 0 is ranked highest
                and K-1 is ranked lowest */

                for(k = 0; k < *K; k++)
                    rank[k] = 0;

                for(k = 0; k < *K; k++)
                {
                    for(k2 = k + 1; k2 < *K; k2++)
                    {
                        if(d[k2 + *K * i] < d[k + *K * i])
                            rank[k]++;
                        else
                            rank[k2]++;
                    }
                }

                double sum = 0;
                for(k = 0; k < *K; k++)
                {
                    f[k + *K * i] = exp(-1.0 * (double) rank[k] / lambda);
                    sum += f[k + *K * i];
                }
                
                for(k = 0; k < *K; k++)
                    f[k + *K * i] /= sum;
            }
        }
        else if(*method == 1)
        {
            memcpy(f_old, f, sizeof(double) * *K * *N);

            /* update entries in f => redefine clusters ! */
            
            for(i = 0; i < *N; i++)
            {
                double sum = 0;
                for(k = 0; k < *K; k++)
                {
                    summa = 0;
                    for(k2 = 0; k2 < *K; k2++)
                        summa += pow(d[k + *K * i] / d[k2 + *K * i], 2.0 / (*fuzzifier - 1.0));
                    
                    f[k + *K * i] = 1.0 / summa;
                    
                    if(*normalizeAssignmentsKFCM)
                        sum += f[k + *K * i];
                }
                
                if(*normalizeAssignmentsKFCM)
                {
                    for(k = 0; k < *K; k++)
                        f[k + *K * i] /= sum;                    
                }
            }
        }
        else if(*method == 2)
        {
            /* update entries in f => redefine clusters ! */

            getHardAssignments(f, f_old, *N, *K, d);
        }
                
        /* test break condition(s) */
        
        if(t >= *T)
            converged = 1;
        else if(*method == 1)
        {
            // Maximum absolute difference between memberships
            double dist = -INFINITY;
            for(i = 0; i < *N; i++)
            {
                for(k = 0; k < *K; k++)
                {
                    double dtmp = fabs(f_old[k + *K * i] - f[k + *K * i]);
                    if(dtmp > dist)
                        dist = dtmp;
                }
            }

            if(*verbose)
                printf("kernelclustering: (assignment) err = %lf\n", dist);
    
            if(dist < *epsilon)
            {
                if(*verbose)
                    printf("kernelclustering: break condition: SE between assignments f (%lf) < epsilon (%lf) at iteration t (%i)\n", dist, *epsilon, t);
                
                converged = 1;
            }
        }
        else if(*method == 2)
        {
            int diff = 0;
            
            for(k = 0; k < *K; k++)
            {
                for(i = 0; i < *N; i++)
                {
                    if(f_old[k + *K * i] != f[k + *K * i])
                        diff++;
                }
            }

            if(!diff)
            {
                if(*verbose)
                    printf("break condition: assignment of prototypes have not changed. t = %i\n",t );
                
                converged = 1;
            }
        }
        
        if(converged && (*method == 0 || (*method == 1 && *getHardClusteringsKFCM)))
        {
            getHardAssignments(f, NULL, *N, *K, d);

            /* get hard assignments */
            
            checkEmptyClusters = 1;
        }
    } /* while(1) */
    
    if(*bSupervised)
    {
        // Get the contingency matrix

        getContingency(*N, *K, f, cMat, l, 0);
    }
}

static long randValLong (
    long min,
    long max)
{
    long val = random () % (max - min + 1) + min;
    
    return val;
}

static void getContingency(int N, int K, double *f, double *cMat, int *l, int Nstart)
{
    int i, k;
    
    for(i = Nstart; i < N; i++)
    {
        double maxval = f[0 + K * i];
        int kstar = 0;
        
        for(k = 1; k < K; k++)
        {
            if(f[k + K * i] > maxval)
            {
                maxval = f[k + K * i];
                kstar = k;
            }
        }

        // Sample i is (best) hard assigned to cluster kstar and has class-label l[i] >= 0
        cMat[kstar + K * l[i]]++;
    }            
}

static void getHardAssignments(double *f, double *f_old, int N, int K, double *d)
{
    int i, k;
    
    for(i = 0; i < N; i++)
    {
        for(k = 0; k < K; k++)
        {
            if(f_old != NULL)
                f_old[k + K * i] = f[k + K * i];

            f[k + K * i] = 0;            
        }
    }

    for(i = 0; i < N; i++)
    {
        double minval = d[0 + K * i];
        int kstar = 0;
        for(k = 1; k < K; k++)
        {
            double val = d[k + K * i];

            if (val < minval)
            {
                kstar = k;
                minval = val;
            }
        }

        f[kstar + K * i] = 1;        
    } /* for i */
}

static int updateAlpha(double labelledUnlabelledFactor, double nUnlabelledSamples, double nLabelledSamples, double *nLabelledSamplesClass, int N, int Nstart, int K, int M, double *f, double *alpha, int *l, int *winner_label, double *winner_count)
{
    int i, j, k, k2;
    double cMat[K * M];
    
    for(k = 0; k < K; k++)
    {
        for(j = 0; j < M; j++)
            cMat[k + K * j] = 0;

        // Initialize alpha
        for(i = 0; i < N; i++)
            alpha[k + K * i] = 1.0;
    }

    // Get the contingency matrix

    getContingency(N, K, f, cMat, l, Nstart);

    for(k = 0; k < K; k++)
    {
        winner_label[k] = 0;
        winner_count[k] = 0;
    }

    int t2;

    int hasConverged = 0;
    
    int F[K];
    int Ftilde[K];
    
    for(k = 0; k < K; k++)
    {
        j = randValLong(1,(M) - 1);
        F[k] = j;
    }
    
    for(t2 = 0; t2 < MAX_ITERATIONS; t2++)
    {
        for(k = 0; k < K; k++)
            Ftilde[k] = F[k];

        /* calculate energy */
        double enLast = 0;
        for(j = 1; j < M; j++)
        {
            double summ = 0;
            for(k = 0; k < K; k++)
            {
                if(F[k] == j)
                    summ += cMat[k + K * j];
            }

            enLast += sqrt(summ);
        }
        
        for(k = 0; k < K; k++)
        {
            int bestJ = F[k];
            double bestEn = 0;
            for(j = 1; j < M; j++)
            {
                Ftilde[k] = j;
                
                double summ = 0;
                double summ2 = 0;
                double summ3 = 0;
                double summ4 = 0;
                
                for(k2 = 0; k2 < K; k2++)
                {
                    if(F[k2] == j)
                        summ += cMat[k2 + K * j];

                    if(F[k2] == F[k])
                        summ2 += cMat[k2 + K * F[k]];

                    if(Ftilde[k2] == j)
                        summ3 += cMat[k2 + K * j];
                
                    if(Ftilde[k2] == F[k])
                        summ4 += cMat[k2 + K * F[k]];
                }
        
                double enDelta = - sqrt(summ) - sqrt(summ2) + sqrt(summ3) + sqrt(summ4);
                
                if(F[k] == j && enDelta != 0)
                {
                    printf("### enDelta = (%lf) != 0, summ (%lf), summ2 (%lf), summ3 (%lf), summ4 (%lf)\n", enDelta, summ, summ2, summ3, summ4);
                    
                    return -3;
                }

                if(enDelta > bestEn)
                {
                    bestEn = enDelta;
                    bestJ = j;
                }
            }
            
            F[k] = bestJ;
            Ftilde[k] = bestJ;
        }

        /* calculate energy */
        double en = 0;
        for(j = 1; j < M; j++)
        {
            double summ = 0;
            for(k = 0; k < K; k++)
            {
                if(F[k] == j)
                    summ += cMat[k + K * j];
            }

            en += sqrt(summ);
        }
        
        if(fabs(en - enLast) < EPSILON)
        {
            hasConverged = 1;
            break;                        
        }                    
    } /* for t2 */
    
    if(!hasConverged)
    {
        printf("### kernelclustering: The class-to-cluster assignments have not converged (lower epsilon or rise the maximum number of iterations)\n");
        
        for(k = 0; k < K; k++)
            printf("%i ", F[k]);
        printf("\n");
                            
        return -4;
    }
    
    for(k = 0; k < K; k++)
    {
        if(winner_label[k])
            continue;
        
        winner_label[k] = F[k];
        winner_count[k] = cMat[k + K * winner_label[k]];
    }
    
    for(k = 0; k < K; k++)
    {
        double winnerImpact = 0;
        
        winnerImpact = 1.0 + labelledUnlabelledFactor * (double) nUnlabelledSamples / (double) nLabelledSamplesClass[winner_label[k]];
        
        for(i = Nstart; i < N; i++)
        {
            if(l[i] == 0 || winner_label[k] == 0)
                alpha[k + K * i] = 1.0;
            else if(winner_label[k] == l[i])
                alpha[k + K * i] = winnerImpact;
            else if(winner_label[k] != l[i])
                alpha[k + K * i] = 0;
        }
    } /* for k */            
    
    return 0;
}        
