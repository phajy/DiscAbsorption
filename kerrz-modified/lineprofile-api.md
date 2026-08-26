# File-based line profile C API

Thin library wrapper around the CLI `lineprof --emissivity-profile` path: load a precomputed emissivity FITS file, build transfer functions, integrate, return `(g, flux)` in memory. No corona tracing inside the library.

## What happens inside

`krz_tool_Lineprofile_run`:

1. Open `tool.emissivity_path` (kerrz `emissivity` FITS) → axisymmetric ε(r)
2. Build Cunningham transfer-function table for `(metric, x_obs, rin, rout)`
3. Integrate → `g_grid` and `flux`

You precompute emissivity yourself (CLI or the helper script below).

## Files touched (library)

| File | Change |
|---|---|
| `wrappers/kerrz.h` | `RETCODE_INVALID_ARGUMENT`, `krz_tool_Lineprofile`, `krz_Lineprofile` |
| `wrappers/interface.zig` | `defaults` / `run` / `deinit` |
| `test/test_lineprofile.c` | Smoke test |
| `test/fixtures/emissivity_lamppost.fits` | Small fixture for CI |
| `test/build_helpers.zig` | Registers the test |

No changes to core corona / emissivity algorithms.

## Minimal C example

```c
#include <kerrz.h>

int main(void) {
    krz_ThreadPool pool;
    krz_ThreadPool_init(&pool, 0);

    krz_tool_Lineprofile tool = krz_tool_Lineprofile_defaults();
    tool.emissivity_path = "emissivity/emis_a0.998_h5.fits";
    tool.x_obs.th = 30.0 * DEG_TO_RAD;

    krz_Lineprofile profile = {0};
    if (krz_tool_Lineprofile_run(&pool, tool, &profile) != RETCODE_SUCCESS) {
        /* handle error */
    }

    /* profile.g_grid[0 .. num_g-1], profile.flux[0 .. num_g-1] */

    krz_Lineprofile_deinit(&profile);
    krz_ThreadPool_deinit(&pool);
    return 0;
}
```

Build:

```bash
zig build --release=fast lib
```

## Precomputing emissivity grids (local, untracked)

```bash
chmod +x scripts/precompute_emissivity.sh
./scripts/precompute_emissivity.sh
```

Writes `emissivity/emis_a{spin}_h{height}.fits` for the spins/heights listed at the top of the script. Edit those arrays (and `NPHOTONS`) as needed. Keep `emissivity/` and the script untracked if you do not want large FITS in git.

Use a small positive spin (e.g. `0.001`) instead of exact `0.0` — kerrz currently panics in the elliptic integrals for Schwarzschild (`a=0`).

## Memory ownership

`krz_tool_Lineprofile_run` allocates `g_grid` and `flux`. Call `krz_Lineprofile_deinit` when finished. Thread pool is separate (`krz_ThreadPool_init` / `deinit`).

`emissivity_path == NULL` or a bad FITS path returns `RETCODE_INVALID_ARGUMENT`.
