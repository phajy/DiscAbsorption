#include <kerrz.h>
#include <stdio.h>
#include <math.h>

int main() {
    krz_RETCODE retcode = 0;
    krz_KerrMetric metric = krz_kerrMetric(1.0, 0.9982);
    krz_RingCorona ring = {.height = 5.0, .radius = 3.0};
    krz_FourVector x_obs = {.t = 0, .r = 1e6, .th = 0.4, .ph = 0};

    krz_ContinuumRing continuum;
    retcode = krz_traceContinuumRing(&continuum, metric, x_obs, ring);
    if (retcode != RETCODE_SUCCESS) {
        return 1;
    }

    // Interpolate a point of the continuum transfer function.
    krz_ContinuumRingPoint point = krz_interpolate_continuum(&continuum, 0.3);
    printf("g    = %lf\n", point.energyshift);
    printf("dt   = %lf\n", point.delta_t);
    printf("lens = %lf\n", point.dcosd_dcosth);

    krz_ContinuumRing_deinit(&continuum);
    printf("Success\n");
}
