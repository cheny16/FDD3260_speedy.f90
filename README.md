# Parallel Speedy.f90 towards Multi-core Processors
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cheny16/FDD3260_speedy.f90/general.yml) ![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cheny16/FDD3260_speedy.f90/genomp.yml?label=omp%20build)



This is an implementation of parallel [speedy.f90](https://github.com/samhatfield/speedy.f90) software for multi-core processors.

> Original project is from: https://github.com/samhatfield/speedy.f90
> 
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

Most commonly used parameters in our experiments:
```yaml
! nsteps_out = model variables are output every nsteps_out timesteps
! nstdia     = model diagnostics are printed every nstdia timesteps
&params
nsteps_out = 100
nstdia     = 180
/
! start_datetime = the start datetime of the model integration
! end_datetime   = the end datetime of the model integration
&date
start_datetime%year   = 1982
start_datetime%month  = 1
start_datetime%day    = 1
start_datetime%hour   = 0
start_datetime%minute = 0
end_datetime%year     = 1982
end_datetime%month    = 2
end_datetime%day      = 1
end_datetime%hour     = 0
end_datetime%minute   = 0
/
```

## Possible runtime errors
- "Segmentation fault" when omp enabled: try `ulimit -s unlimited`.
- Multiple threads not monitored: `setenv OMP_DYNAMIC false` or `export OMP_DYNAMIC=FALSE`.
