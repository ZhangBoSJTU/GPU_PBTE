# `GPU_PBTE`：An Efficient Solver for Three and Four Phonon Scattering Rates on Graphics Processing Units

## Authors:

* Zhang Bo (Shanghai Jiao Tong University)
  * zhangbo.sjtu@qq.com
* Xiaokun Gu (Shanghai Jiao Tong University)
  * xiaokun.gu@sjtu.edu.cn

## What is `GPU_PBTE` ?

* `GPU_PBTE` stands for Graphics Processing Units accelerated Phonon Boltzmann Translation Equation solver.
* `GPU_PBTE` can take three-phonon scattering process, four-phonon scattering process and isotope scattering process into consideration and has been significantly accelerated by using GPUs. 
   
## Prerequisites

* You need to have a GPU card, `CUDA` toolkit and NVIDIA HPC SDK.
* It works for Linux (with GCC) operating systems. 
* Intel MKL(Math Kernel Library) is need for this program.

## Compile `GPU_PBTE`

* In order to link against Atsushi Togo's [spglib](http://spglib.sourceforge.net/), `external/lib` is included in the root directory.
`nvhpc/20.11` is known to work.
* Environment
  ```bash
  module purge
  module load intel-mkl
  module load nvhpc/20.11-gcc-4.8.5 
  export LD_LIBRARY_PATH="./external/spglib-1.7.4-intel16/lib:$LD_LIBRARY_PATH"
  ```
* Go to the `PBTE` root directory and type `make help` to show help documentation and get started!

## Run `GPU_PBTE`

* Go to the `Test-Si` directory.
* Modify the `input/input.dat` file as needed.
* Modify the scripts as needed. 
* `Job1` is to mesh the Brillouin zone and obtain the eigenvector, frequency, group velocity and other physical quantities based on the diagonalized dynamical matrix.
* `Job2` is to calculate the phonon scattering processes.
* `job3` is to solve the thermal conductivity.

## Input files

Exactly five files are required for a `GPU_PBTE` run: `input.txt`, `HFC.dat` , `CFC.dat`,  `QFC.dat ` and `jobtype.txt`. Their contents are detailed below.

### The `input.txt` file

The contents of this file describe the system to be studied and specify a set of parameters and flags controlling execution. Here is an example of the silicon input.

```bash
.true. 20             # .true. —— RTA method    .false. —— Iterative method
.true.                # three-phonon sacttering processes
.true.                # four-phonon sacttering processes
.true.                # isotope sacttering processes
10.2621               # lattice constant
0 0.5 0.5             # vectors for unit cell
0.5 0 0.5
0.5 0.5 0
2                     # number of basis
28.085 28.085         # atomic masses
1 1                   # atom type
1.529e-4 1.529e-4     # isotope g2
0.0   0.0   0.0       # positions for bases
0.25  0.25  0.25
16 16 16              # number of q points to sample
1
300                   # temperature list
.true. .true. .true.  # which directions are computed
1e20                  # sample size
1 1.0 1.0
```

- `&crystal` namelist:
    - `a0` (real): lattice constant
    - `a1,a2,a3` (real): vectors for unit cell
    - `nBasis` (integer): number of basis
    - `mass` (real): atomic masses
    - `types` (integer): atom type
    - `g2` (integer): isotope parameter g2
    - `pos` (real): positions for bases
- `mesh` namelist:
    - `nQx, nQy, nQz` (integer): number of q points to sample
- `&flags` namelist:
    - `RTA` (logical): True for RTA method and false for iterative method.
    - `isCubic` (logical): whether calculate three-phonon sacttering processes or not.
    - `isQuart` (logical): whether calculate four-phonon sacttering processes or not.
    - `isIso` (logical): whether calculate isotope sacttering processes or not.
    - `isKx, isKy isKz` (logical): which directions are computed.

### The `HFC.dat` file

This file contains the second derivatives of the system's energy with respect to the Cartesian coordinates of the nuclei. The first line of the file declares the total number of HFC.

The following is an example of part of HFC.dat:

```txt 
2130
-1 -1 0 1 2 1 1    0.000674324651356
-1 -1 0 1 2 1 2   -0.000680759987194
-1 -1 0 1 2 1 3    0.000524309131218
-1 -1 0 1 2 2 1   -0.000680759987194
-1 -1 0 1 2 2 2    0.000674324651356
-1 -1 0 1 2 2 3    0.000524309131218
-1 -1 0 1 2 3 1    0.000524309131218
      ...                 ...
```

### The `CFC.dat` file

Similarly, this file contains the third-order interatomic force constant matrix, but uses a sparse description to save space. All constants are implicitily refered to a central unit cell i taken as the origin of coordinates. The first line of the file declares the total number of CFC.

- A line with the Cartesian coordinates labling the unit cell
- 27 force constant data

The following is an example of part of CFC.dat:

```txt
842
-1 0 0 -1 0 0 1 2 2
   0.033127780731527   -0.066057656888135   -0.066057656888135   -0.066057656888135    0.064979717826671    0.090454846582396   -0.066057656888135    0.090454846582396    0.064979717826671   -0.064979717826671    0.066057656888135    0.090454846582396    0.066057656888135   -0.033127780731527   -0.066057656888135    0.090454846582396   -0.066057656888135   -0.064979717826671   -0.064979717826671    0.090454846582396    0.066057656888135    0.090454846582396   -0.064979717826671   -0.066057656888135    0.066057656888135   -0.066057656888135   -0.033127780731527 
-1 0 0 0 0 0 1 2 1
  -0.033127780731527    0.066057656888135    0.066057656888135    0.064979717826671   -0.066057656888135   -0.090454846582396    0.064979717826671   -0.090454846582396   -0.066057656888135    0.066057656888135   -0.064979717826671   -0.090454846582396   -0.066057656888135    0.033127780731527    0.066057656888135   -0.090454846582396    0.064979717826671    0.066057656888135    0.066057656888135   -0.090454846582396   -0.064979717826671   -0.090454846582396    0.066057656888135    0.064979717826671   -0.066057656888135    0.066057656888135    0.033127780731527 
0 -1 0 0 -1 0 1 2 2
  -0.033127780731527    0.066057656888135   -0.066057656888135    0.066057656888135   -0.064979717826671    0.090454846582396   -0.066057656888135    0.090454846582396   -0.064979717826671    0.064979717826671   -0.066057656888135    0.090454846582396   -0.066057656888135    0.033127780731527   -0.066057656888135    0.090454846582396   -0.066057656888135    0.064979717826671   -0.064979717826671    0.090454846582396   -0.066057656888135    0.090454846582396   -0.064979717826671    0.066057656888135   -0.066057656888135    0.066057656888135   -0.033127780731527 
... ... 
```

### The `QFC.dat` file

Similarly, this file contains the forth-order interatomic force constant matrix, but uses a sparse description to save space. All constants are implicitily refered to a central unit cell i taken as the origin of coordinates. The first line of the file declares the total number of QFC.

- A line with the Cartesian coordinates labling the unit cell
- 81 lines, each line contains a force constant

The following is an example of part OF QFC.dat:

```txt
58
0 0 0 0 0 0 0 0 0 1 1 1 1
  -0.220186929483593
   0.000000000000000
   0.000000000000000
   0.000000000000000
   0.273929132908144
   0.000000000000000
   0.000000000000000
          ...
   0.000000000000000
   0.000000000000000
   0.273929132908144
   0.000000000000000
   0.000000000000000
   0.000000000000000
  -0.220186929483593
```

## Output files

Many files and directories are created during a successful run of `GPU_PBTE`. They contain not only the thermal conductivity and related quantities, but also a set of intermediate results that may be useful to diagnose problems.

* `Job 1`
  - `grpvelx.dat, grpvely.dat, grpvelz.dat `: phonon group velocity
  - `mesh.dat `: __q__-point
  - `x.z`: where x is a 1-based sequential index as a flag to record whether this __q__-point is calculated or not

* `Job 2`
  - `Q_x_y.dat `: phonon scattering rates for x __q__-point at y-th temperature for RTA method
  - `A_x_y.dat `: phonon scattering rates for x __q__-point at y-th temperature for iterative method

* `Job 3`
  - `rt_y.dat `: phonon relaxation time y-th temperature
  - `TC_k_inf.dat`: thermal conductivity at k temperature

## Website:

* https://gitlab.com/xiaokun.gu/GPU_PBTE