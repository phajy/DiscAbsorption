#include <kerrz.h>
#include <math.h>
#include <stdio.h>

#define FIXTURE_PATH "test/fixtures/emissivity_lamppost.fits"

static int finite_positive_sum(const krz_Lineprofile *lp) {
    double sum = 0;
    for (size_t i = 0; i < lp->num_g; ++i) {
        if (!isfinite(lp->flux[i]))
            return 0;
        sum += lp->flux[i];
    }
    return sum > 0;
}

int main(void) {
    krz_RETCODE retcode;
    krz_ThreadPool pool;
    krz_tool_Lineprofile tool = krz_tool_Lineprofile_defaults();
    krz_Lineprofile profile = {0};

    tool.emissivity_path = FIXTURE_PATH;
    tool.nradii = 20;
    tool.nangles = 50;
    tool.ng = 100;
    tool.ngstar = 200;
    tool.nrsteps = 200;
    tool.r_out = 100;

    retcode = krz_ThreadPool_init(&pool, 1);
    if (retcode != RETCODE_SUCCESS)
        return 1;

    retcode = krz_tool_Lineprofile_run(&pool, tool, &profile);
    if (retcode != RETCODE_SUCCESS) {
        fprintf(stderr, "krz_tool_Lineprofile_run failed: %d\n", retcode);
        krz_ThreadPool_deinit(&pool);
        return 1;
    }

    if (profile.num_g == 0 || !finite_positive_sum(&profile)) {
        krz_Lineprofile_deinit(&profile);
        krz_ThreadPool_deinit(&pool);
        return 1;
    }

    printf("line profile bins: %zu\n", profile.num_g);

    krz_Lineprofile_deinit(&profile);
    krz_ThreadPool_deinit(&pool);
    return 0;
}
