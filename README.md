# oomph-thin-film-control

An academic research code for simulating a gravity-driven thin liquid film on
an inclined plane and applying an LQR controller to stabilise its flat free
surface. It is built on [oomph-lib](https://oomph-lib.github.io/).

## Requirements

- A working build of oomph-lib, including the linear-algebra dependencies used
  by its build (the supplied Makefile is configured for MUMPS/SuperLU).
- A C++17 compiler and `make`.
- Python 3 with NumPy and Matplotlib to produce the supplied plots.

Set `OOMPHLIB` to the root of the oomph-lib checkout/build before compiling:

```sh
export OOMPHLIB=/path/to/oomph-lib
```

The Makefile contains machine-specific linker paths. The macOS section assumes
Homebrew GCC and LAPACK; the non-macOS section was written for a particular
MPI cluster. Adapt the compiler, include, library, and MPI paths there to
match the oomph-lib installation on your machine.

## Build and run

From the repository root, create the output directory and build the program:

```sh
mkdir -p output
make
```

Run a small test case first:

```sh
./main --nx 8 --ny 8 --nz 2 --nxcontrol 8 --nycontrol 8 \
  --mcontrol 4 --pcontrol 4 --tburn 1 --tcontrol 1 \
  --dtburn 0.1 --dtcontrol 0.1
```

The default problem is considerably larger (`40 x 40 x 3` fluid elements,
with a `40 x 40` control grid) and runs for 3000 time units. Run it with:

```sh
./main
```

The driver currently requires `--dtburn` and `--dtcontrol` to have identical
values. Use `./main --help` to see the available parameters. In particular,
domain size (`--lx`, `--ly`), physical parameters (`--re`, `--ca`,
`--theta`), mesh resolution, control-grid resolution, number of actuators, and
run lengths can all be changed from the command line. The initial perturbation
and control method are currently selected in `src/main.cpp`.

## Results and plots

Each run writes `output/surface_*.dat`: one free-surface/control-grid snapshot
per timestep. Rows contain spanwise position, streamwise position, film height,
two flux components, and applied forcing. When control is initialised, the
program also writes controller and model matrices such as `K.dat`, `A_wr.dat`,
`B_wr.dat`, and `F.dat`.

Generate the supplied diagnostic plots from the repository root:

```sh
python3 gen_plots.py
```

This writes `output/hmax.png`, which shows the maximum deviation from the flat
film, and plots of the control/model matrices. `gen_plots.py` assumes the
default `40 x 40` control grid. If you used another control-grid size, change
its `nx` and `ny` values at the top of the script before plotting the matrix
fields. The `post-processing/` directory also contains optional Gnuplot scripts
for surface and mesh animations.

To discard generated build products and output data:

```sh
make clean
```
