#ifndef _KERRZ_H
#define _KERRZ_H

#include <unistd.h>

#define DEG_TO_RAD 3.14159265358979323846 / 180.0

// Return codes that a kerrz function may return.
typedef enum _krz_RETCODE {
    RETCODE_SUCCESS = 0,
    RETCODE_ALLOCATION_FAILED,
    RETCODE_THREAD_ERROR,
    RETCODE_CONVERGENCE_FAILED,
    RETCODE_INVALID_ARGUMENT
} krz_RETCODE;

// Status codes that are associated with a geodesic.
typedef enum _krz_STATUS {
    STATUS_NONE,
    STATUS_EVENT_HORIZON,
    STATUS_INTERSECTED_DISC,
    STATUS_AT_INFINITY,
} krz_STATUS;

// The radial case of the geodesic, as defined by the radial turning points.
typedef enum _krz_RADIAL_CASE {
    RADIAL_CASE_I = 1,
    RADIAL_CASE_II = 2,
    RADIAL_CASE_III = 3,
    RADIAL_CASE_IV = 4,
} krz_RADIAL_CASE;

// The angular case of the geodesic, as defined by the angular (theta) turning
// points.
typedef enum _krz_ANGULAR_CASE {
    ANGULAR_CASE_NORMAL = 1,
    ANGULAR_CASE_VORTICAL = 2,
} krz_ANGULAR_CASE;

// A four-vector in the Boyer--Lindquist coordinates.
typedef struct _krz_FourVector {
    double t;
    double r;
    double th;
    double ph;
} krz_FourVector;

// A thread (worker) pool container.
typedef struct _krz_ThreadPool {
    void *pool;
    size_t num_threads;
} krz_ThreadPool;

// Initialise a thread pool. `krz_ThreadPool_deinit` must be called to cleanup
// resources allocated by this function.
//
// If `n_threads == 0`, will use as many threads as there are CPU cores
// available.
//
// The return value indicates success if 0, else returns an error code from
// `krz_RETCODE`.
krz_RETCODE krz_ThreadPool_init(krz_ThreadPool *pool, size_t n_threads);

// Cleanup a thread pool.
void krz_ThreadPool_deinit(krz_ThreadPool *pool);

// The Kerr geometry metric. Use `krz_kerrMetric` (lowercase k) to construct a
// metric from the mass and spin.
typedef struct _krz_KerrMetric {
    double M;
    double a;
    double horizon_radius;
    double horizon_radius_negative;
    double isco;
} krz_KerrMetric;

// Construct a Kerr metric. Given a mass and spin (in dimensionless units),
// calculates the fields of the `krz_KerrMetric`.
krz_KerrMetric krz_kerrMetric(double M, double a);

// A ring-like coronal model, parameterised by a height and radius in
// gravitational units.
typedef struct _krz_RingCorona {
    double height;
    double radius;
} krz_RingCorona;

// An enumeration of all coronal models.
typedef enum _krz_CORONA_MODEL { RING_CORONA } krz_CORONA_MODEL;

// A union of all coronal models.
typedef struct _krz_CoronaModel {
    union {
        krz_RingCorona ring;
    } as;
    krz_CORONA_MODEL tag;
} krz_CoronaModel;

// Initilise a ring coronal model.
static krz_CoronaModel krz_ringCorona(double height, double radius) {
    krz_CoronaModel out;
    out.tag = RING_CORONA;
    out.as.ring.height = height;
    out.as.ring.radius = radius;
    return out;
}

// Null-geodesic initial conditions. This must be used in the following way:
// ```c
// krz_KerrMetric metric = krz_kerrMetric(1.0, 0.998);
// krz_FourVector x_init = {0};
// x_init.r = 1e4;
// x_init.th = 0.3;
//
// krz_InitialConditions conds;
// // For example:
// krz_fromImpactParameters(&conds, metric, x_init, 2.0, -2.0);
// ```
typedef struct _krz_InitialConditions {
    // The geodesic energy.
    double E;
    // The geodesic angular momentum.
    double L;
    // The geodesic Carter constant.
    double Q;
    // The reduced geodesic angular momentum L / E.
    double lambda;
    // The reduced geodesic Carter constant Q / E^2.
    double eta;
    // The initial four position.
    krz_FourVector x_init;
    // The signs of the momenta and the (maximum) number of windings the
    // geodesic should follow.
    double theta_sign, radial_sign, windings;
} krz_InitialConditions;

typedef struct _krz_OrthonormalFrame {
    krz_FourVector x;
    // The g_(mu, nu) metric components: [tt, rr, θθ, φφ, tφ]
    double metric_components[5];
    // This is an implementation-specific matrix that contains the local basis
    // (or frame) vectors.
    double matrix[16];
} krz_OrthonormalFrame;

// Calculate a stationary frame at a point in the spacetime.
krz_OrthonormalFrame
krz_stationaryFrame(krz_KerrMetric metric, krz_FourVector x_init);

// Calculate a locally non-rotating frame at a point in the spacetime.
krz_OrthonormalFrame
krz_lnrFrame(krz_KerrMetric metric, krz_FourVector x_init);

// Calculate a frame with a particular four-velocity at a point in the
// spacetime.
krz_OrthonormalFrame
krz_frame(krz_KerrMetric metric, krz_FourVector x_init, krz_FourVector v_frame);

// Calculate the four-velocity of a circular orbit in the equatorial at a
// particular radius. Uses a free-falling prescription for the sub-ISCO radii.
krz_FourVector krz_circularOrbitVelocity(krz_KerrMetric metric, double r);

// Calculate initial conditions from impact parameters.
krz_InitialConditions krz_fromImpactParameters(
    krz_KerrMetric metric,
    krz_FourVector x_init,
    double alpha,
    double beta
);

// Calculate initial conditions from local sky angles projected in a tangent
// space.
krz_InitialConditions krz_fromSkyAngles(
    krz_KerrMetric metric,
    krz_OrthonormalFrame frame,
    double theta,
    double phi
);

// The result of a kerrz geodesic trace.
typedef struct _krz_TraceResult {
    // The status of the geodesic.
    krz_STATUS status;
    // The integrated Mino time at this point.
    double mino_time;
    // Coordinates of the geodesic at this point.
    krz_FourVector x_final;
    // The turning point winding number.
    double winding;
} krz_TraceResult;

// Trace a single geodesic from a set of initial conditions to a particular
// angle.
krz_TraceResult krz_traceToAngle(
    krz_KerrMetric metric,
    krz_InitialConditions init_conds,
    double angle
);

// Trace a single geodesic from a set of initial conditions to a particular
// radius.
krz_TraceResult krz_traceToRadius(
    krz_KerrMetric metric,
    krz_InitialConditions init_conds,
    double radius
);

// The result of a kerrz continuum trace for the lamppost model. Since the
// lamppost is entirely axis-symmetric, this is a single unique geodesic.
typedef struct _krz_ContinuumLamppost {
    // The resulting photon.
    krz_TraceResult res;
    // The emission angle (`delta`) that solves this geodesic.
    double angle_delta;
    // The lensing factor, `cos(delta) / cos(theta)`, where delta is the polar
    // angle on the local sky of the lamppost, and `theta` is the observer
    // inclination.
    double dcosd_dcosth;
    // The impact parameters of this photon.
    double alpha, beta;
} krz_ContinuumLamppost;

// Trace a single geodesic from a set of initial conditions to a particular
// radius.
krz_ContinuumLamppost krz_traceContinuumLamppost(
    krz_KerrMetric metric,
    krz_FourVector x_obs,
    double height
);

// The result of a kerrz continuum trace for the ring-like coronal model
typedef struct _krz_ContinuumRing {
    void *cache;
} krz_ContinuumRing;


void krz_ContinuumRing_deinit(krz_ContinuumRing *output);

// Trace a single geodesic from a set of initial conditions to a particular
// radius.
krz_RETCODE krz_traceContinuumRing(
    krz_ContinuumRing *output,
    krz_KerrMetric metric,
    krz_FourVector x_obs,
    krz_RingCorona ring
);

typedef struct _krz_ContinuumRingPoint {
    // The corona-to-observer light travel time.
    double delta_t;
    // The corona-to-observer energyshift.
    double energyshift;
    // The lensing factor term.
    double dcosd_dcosth;
} krz_ContinuumRingPoint;

krz_ContinuumRingPoint krz_interpolate_continuum(krz_ContinuumRing *continuum, double phi);

typedef struct _krz_PathBuilder {
    void *state;
} krz_PathBuilder;

// Initialise a path builder for the given metric parameters and initial
// conditions.
krz_RETCODE krz_PathBuilder_init(
    krz_PathBuilder *pb,
    krz_KerrMetric metric,
    krz_InitialConditions init_conds
);

// Free the resources allocated for the `krz_PathBuilder`.
void krz_PathBuilder_deinit(krz_PathBuilder *pb);

// Calculate the geodesic from a path builder up to a given Mino time.
krz_TraceResult krz_at_mino_time(krz_PathBuilder *pb, double mino_time);

typedef struct _krz_TurningPoints {
    double r_0, r_1;
    double theta_0, theta_1;
} krz_TurningPoints;

// Return the Mino time to the turning points.
krz_TurningPoints krz_mino_time_to_turning_points(krz_PathBuilder *pb);

// Used to store a single evaluation of a (Cunningham) transfer function table.
typedef struct _krz_TransferFunctionCache {
    void *cache;
} krz_TransferFunctionCache;

// Used to free resources allocated by the transfer function cache.
void krz_TransferFunctionCache_free(krz_TransferFunctionCache *cache);

// Calculate (Cunningham) transfer functions for a particular metric and
// observer position between some inner and outer radius on the accretion disc.
// Calculates `num_radii` worth of transfer functions. By default, these are
// linearly spaced between `r_min` and `r_max`.
krz_STATUS krz_transfer_functions(
    krz_TransferFunctionCache *cache,
    krz_KerrMetric metric,
    krz_FourVector observer,
    double r_min,
    double r_max,
    size_t num_radii
);

// A point on the observer plane that is interpolated from a (Cunningham)
// transfer function.
typedef struct _krz_InterpolatedPoint {
    // Impact parameters.
    double alpha, beta;
    // Disc to observer time with the observer distance subtracted.
    double delta_t;
} krz_InterpolatedPoint;

// Interpolate a point on the image plane given a set of coordinates on the
// accretion disc.
krz_InterpolatedPoint krz_interpolate_disc_coordinates(
    krz_TransferFunctionCache *cache,
    double r,
    double phi
);

// Memory and other pre-allocations and pre-calculations used to calculate
// emissivity profiles.
typedef struct _krz_EmissivityCache {
    void *context;
} krz_EmissivityCache;

// Initialise an EmissivityCache for up to `num_traces` traces.
krz_RETCODE krz_EmissivityCache_init(
    krz_EmissivityCache *ctx,
    size_t num_traces
);

// Free the resources allocated with `krz_EmissivityCache_init`.
void krz_EmissivityCache_deinit(krz_EmissivityCache *ctx);

// The values that can be interpolated from an emissivity calculation.
typedef struct _krz_EmissivityTrace {
    // The emissivity itself.
    double em;
    // The corona-to-disc light travel time.
    double t;
    // The corona-to-disc energyshift.
    double g;
    // The incident angle between the photon and the disc surface. A value of
    // `0` is directly parallel with the disc normal at that radius.
    double local_theta;
} krz_EmissivityTrace;

// Calculate an emissivity profile for a particular coronal model, making use
// of the thread pool if needed.
krz_RETCODE krz_emissivity(
    krz_ThreadPool *pool,
    krz_EmissivityCache *ctx,
    krz_KerrMetric metric,
    krz_CoronaModel model
);

// Calculate the emissivity as in `krz_emissivity` but explicitly for a
// ring-like corona.
krz_RETCODE krz_emissivity_ring(
    krz_ThreadPool *pool,
    krz_EmissivityCache *ctx,
    krz_KerrMetric metric,
    krz_RingCorona ring
);

// Interpolate the emissivity profile at a particular radius and azimuthal
// coordinate on the disc.
krz_EmissivityTrace krz_interpolate_emissivity(krz_EmissivityCache *cache, double r, double phi);

// Calculate the impact parameters that trace the shadow of the black hole.
double krz_shadow(
    krz_KerrMetric metric,
    double observer_incl,
    double *alpha,
    double *beta,
    size_t num_pts
);

// Below are bindings to all of the different kerrz.tools. These are WIP, and
// hence only a subset are supported.

// Return codes that a kerrz function may return.
typedef enum _krz_TOOL {
    TOOL_TRANSFER_FUNCTION = 0,
} krz_TOOL;

// Return codes that a kerrz function may return.
typedef enum _krz_TRANSFER_FUNCTION_HEURISTIC {
    TF_HEURISTIC_NONE = 0,
    TF_HEURISTIC_ARCLEN,
    TF_HEURISTIC_IMPACT,
} krz_TRANSFER_FUNCTION_HEURISTIC;

typedef struct _krz_tool_TransferFunction {
    krz_TOOL tag;
    krz_KerrMetric metric;
    krz_FourVector x_obs;
    double r_target;
    struct {
        size_t max_points;
        size_t refine_N;
        size_t refine_M;
        size_t optimise;
        double minimum_guess;
        krz_TRANSFER_FUNCTION_HEURISTIC heuristic;
    } options;
} krz_tool_TransferFunction;

// Obtain default values for the transfer function tool.
krz_tool_TransferFunction krz_tool_TransferFunction_defaults();

typedef struct _krz_CunninghamTrace {
    // The offset angle on the image plane.
    double image_angle;
    // The offset radius on the image plane. Together with `image_angle` they
    // describe the impact parameters.
    double image_radius;
    // The alpha impact parameter.
    double alpha;
    // The beta impact parameter.
    double beta;
    // The redshift value of this trace.
    double g;
    // The normalised redshift value of this trace.
    double g_star;
    // The azimuthal coordinate of this trace.
    double phi;
    // Light-travel time with the observer distance subtracted.
    double delta_t;
    // The |d(r, g)/d(alpha, beta)| Jacobian term, specifically the
    // absolute value of the determinant.
    double jacobian;
    // Cunningham's actual transfer function, or the value thereof.
    double f;
    // The root solver error.
    double r_err;
} krz_CunninghamTrace;

typedef struct _krz_CunninghamTransferFunction {
    double g_min;
    double g_max;
    krz_CunninghamTrace *traces;
    size_t num_traces;
} krz_CunninghamTransferFunction;

// Free resources associated with the krz_CunninghamTransferFunction.
void krz_CunninghamTransferFunction_deinit(krz_CunninghamTransferFunction *ctf);

// Run the transfer function tool.
krz_RETCODE krz_tool_TransferFunction_run(
    krz_tool_TransferFunction tool,
    krz_CunninghamTransferFunction *result
);

// Parameters for integrating a line profile from a precomputed emissivity
// FITS file (the same format written by `kerrz emissivity`).
typedef struct _krz_tool_Lineprofile {
    krz_KerrMetric metric;
    krz_FourVector x_obs;
    const char *emissivity_path;
    double r_in;
    double r_out;
    size_t nradii;
    size_t nangles;
    size_t ng;
    size_t ngstar;
    size_t nrsteps;
    int normalise;
} krz_tool_Lineprofile;

krz_tool_Lineprofile krz_tool_Lineprofile_defaults(void);

typedef struct _krz_Lineprofile {
    double *g_grid;
    double *flux;
    size_t num_g;
} krz_Lineprofile;

void krz_Lineprofile_deinit(krz_Lineprofile *lp);

krz_RETCODE krz_tool_Lineprofile_run(
    krz_ThreadPool *pool,
    krz_tool_Lineprofile tool,
    krz_Lineprofile *result
);

#endif
