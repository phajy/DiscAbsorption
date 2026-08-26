#include <kerrz.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_IMPACTS 1000

int main() {
    krz_KerrMetric metric = krz_kerrMetric(1.0, 0.9982);

    double *mem = malloc(sizeof(double) * 2 * NUM_IMPACTS);
    if (mem == NULL) {
        return 1;
    }

    double *alpha = mem;
    double *beta = (mem + NUM_IMPACTS);

    krz_shadow(metric, 9.99 * DEG_TO_RAD, alpha, beta, NUM_IMPACTS);

    FILE *fd = fopen("shadow.test", "w");
    for (int i = 0; i < NUM_IMPACTS; ++i) {
        fprintf(fd, "%.5f, %.5f\n", alpha[i], beta[i]);
    }
    fclose(fd);
    return 0;
}
