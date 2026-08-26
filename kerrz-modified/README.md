
<img width="30%" src=".forgejo/kerrz-logo.svg"> GPL 3.0

kerrz is a general relativistic 'ray-tracing' library and command line tool.
Note that 'ray-tracing' is in quotation marks, as kerrz actually computes the
(semi-)analytic geodesic formulation of Gralla and Lupsasca, 2020
(arXiv1910.12881v3), and so spends very little time actually ray-tracing in the
conventional sense.

Additionally, kerrz is written with a Zig automatic-differentation library that
implements simple forward-mode autodiff. Everything that kerrz traces can
therefore have derivatives and Jacobian terms calculated for them too.

This program is being developed primarily for research purposes.

## Quickstart

- The easiest way to install kerrz is via the Python package manager:
  ```bash
  pip install kerrz-cli
  ```
  This will make `kerrzcli` available on your path. It has a different name from
  the binaries distributed in other channels to make provenance easy to
  differentiate, and to avoid clobbering the executable name if installing other
  Python related kerrz tools. Feel free to alias.

- Alternatively, you can grab binaries, header / module files, or libraries
  from the [releases](https://codeberg.org/astro-group/kerrz/releases).

To get started, try
```
kerrz image
```
to render an image of a black hole. Open the `output.pgm` file to view the output.

For more, use the builtin help, either `kerrz --help image`, or simply:
```
$ kerrz --help
kerrz GPL 3.0 version: 0.3.0+f2ec0ab48c51027cfd0465d93e66bc85382ed7ec
https://git.sr.ht/~fjebaker/kerrz

Usage: kerrz {command} [--flags and positional arguments]

General flags:
    [--help]                  Print this help message or help for a specific command.
    [--version]               Print version information.

Commands:
- calc        A simple calculator.
- continuum   Calculate continuum profiles.
- emissivity  Calculate emissivity profiles for a corona model.
- grid        Generate grids.
- image       Render simple images of black holes.
- impulse     Calculate impulse responses.
- lineprof    Calculate line profiles.
- sky         Render the sky of a point in the spacetime.
- table       Calculate non-trivial tables.
- tf          Calculate transfer functions.
- trace       Trace a single geodesic and print information.
- mapper      Print information about the mapper.
- completion  Generate shell completion helpers.
```

## Compiling

kerrz is written in Zig using compiler version 0.15.2:
```
zig build --release=fast
```
All dependencies are Zig libraries, so static cross-compilation works fine.
You can even compile [kerrz to WASM](https://cosroe.com/kerrz).

