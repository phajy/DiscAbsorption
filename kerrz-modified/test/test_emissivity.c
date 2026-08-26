#include <kerrz.h>
#include <stdio.h>
#include <math.h>

#define NUM_TRACES 500 * 1000

int main() {
    krz_RETCODE retcode = 0;
    krz_KerrMetric metric = krz_kerrMetric(1.0, 0.9982);
    krz_CoronaModel corona = krz_ringCorona(5.0, 3.0);

    // Initilise a thread pool
    krz_ThreadPool pool;
    retcode = krz_ThreadPool_init(&pool, 0);
    if (retcode != RETCODE_SUCCESS)
        return retcode;

    krz_EmissivityCache cache;
    retcode = krz_EmissivityCache_init(&cache, NUM_TRACES);
    if (retcode != RETCODE_SUCCESS)
        goto CLEANUP1;

    // Parallel calculation for the emissivity profile
    retcode = krz_emissivity(&pool, &cache, metric, corona);
    if (retcode != RETCODE_SUCCESS)
        goto CLEANUP2;

    // Interpolate the emissivity at a particular radius
    double r = 5.0;
    double phi = M_PI;
    krz_EmissivityTrace em = krz_interpolate_emissivity(&cache, r, phi);
    printf("em: %lf\n", em.em);
    printf("t : %lf\n", em.t);
    printf("g : %lf\n", em.g);

CLEANUP2:
    krz_EmissivityCache_deinit(&cache);
CLEANUP1:
    krz_ThreadPool_deinit(&pool);
    return retcode;
}
