# Parallel Speedy.f90 towards Multi-core Processors
This is an implementation of parallel [Speedy.f90](https://github.com/samhatfield/speedy.f90) software for multi-core processors.

Original project is from: https://github.com/samhatfield/speedy.f90

> speedy.f90 is an intermediate complexity atmospheric general circulation model written in modern Fortran. It is based on SPEEDY, developed by Fred Kucharski, Franco Molteni and Martin P. King.
> Zenodo DOI: 10.5281/zenodo.5816982

In this project, we implement a parallel speedy.f90 using OpenMP, add supports for building on multiple platforms including Intel® CPU, AMD® CPU, with different compilers:
- GNU Fortran compiler `gfortran`
- Intel Fortran compiler `ifort`
- Intel OneAPI Fortran compiler `ifx`
- Cray Fortran compiler `ftn`

## Installation
Steps for installation are basically the same as the original speedy.f90, but with extensions of builing options.
1. Install the NetCDF library and get the location the install directory.
2. Set the NETCDF environment variable to point to the directory containing the NetCDF include and lib directories.
3. Run build.sh to build speedy.f90: `bash build.sh --target <TARGET>`, some optional targets include:
   
   | Parameter | Compiler       | Platforms  | OpenMP |
   | --------- | -------------- | ---------- | ------ |
   | general   | gfortran       | Intel, AMD | Off    |
   | genomp    | gfortran       | Intel, AMD | On     |
   | intel     | Intel compiler | Intel      | Off    |
   | intelomp  | Intel compiler | Intel      | On     |
   
  For example, to install on Intel CPUs with OpenMP enabled:
  ```bash
  bash build.sh --target intelomp
  ```

## Run
Edit the parameters in [namelist.mnl](https://github.com/cheny16/FDD3260_speedy.f90/blob/main/namelist.nml) if necessary, then `bash run.sh`.

## Possible runtime errors
- "Segmentation fault" when omp enabled: try `ulimit -s unlimited`.
- Multiple threads not monitored: `setenv OMP_DYNAMIC false` or `export OMP_DYNAMIC=FALSE`.
