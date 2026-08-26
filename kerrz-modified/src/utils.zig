const std = @import("std");
const zfits = @import("zfits");
const options = @import("options");

/// Add information about the kerrz version to the HDU header, if it does not
/// already exist.
pub fn addKerrzFITSInfo(hdu: *zfits.Hdu) !void {
    if (hdu.record_keys.get("KZVERSN") == null) {
        const version_string = std.fmt.comptimePrint("{d}.{d}.{d}-{s}", .{
            options.version.major,
            options.version.minor,
            options.version.patch,
            options.version.build.?[0..8],
        });
        try hdu.appendHeaderRecord("KZVERSN", .{
            .value = .{ .string = version_string },
            .comment = "kerrz version",
        });
    }
}

/// Add observer information to HDU.
pub fn addObserverInformation(hdu: *zfits.Hdu, x_obs: anytype) !void {
    try hdu.appendHeaderRecord("R_OBS", .{
        .comment = "The observer radial distance in rg",
        .value = .{ .float = @floatCast(x_obs.r.x) },
    });
    try hdu.appendHeaderRecord("INCL_OBS", .{
        .comment = "The observer inclination in degrees",
        .value = .{ .float = @floatCast(std.math.radiansToDegrees(x_obs.th.x)) },
    });
    try hdu.appendHeaderRecord("AZM_OBS", .{
        .comment = "The observer azimuth in degrees",
        .value = .{ .float = @floatCast(std.math.radiansToDegrees(x_obs.ph.x)) },
    });
}
