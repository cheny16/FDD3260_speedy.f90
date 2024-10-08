#!/bin/bash

# Check that NetCDF variable is set
if [ -z "$NETCDF" ]; then
    echo "\$NETCDF variable not set"
    echo "Aborting build"
    exit 1
fi

# Name of makefile
if [ "$1" = "--target" ]; then
    if [ "$2" = "intel" ]; then
	MAKE=intel.makefile
    fi
    if [ "$2" = "intelomp" ]; then
	MAKE=intel_omp.makefile
    fi
    if [ "$2" = "intelacc" ]; then
	MAKE=intel_acc.makefile
    fi
    if [ "$2" = "general" ]; then
	MAKE=gfortran.makefile
    fi
    if [ "$2" = "genomp" ]; then
	MAKE=gfortran_omp.makefile
    fi
fi

# Define directory names
ROOT=`pwd`
SRC=$ROOT/source
BIN=$ROOT/bin

# Make binary directory
rm -rf $BIN
mkdir $BIN
cd $BIN

# Copy source files
cp $SRC/*.f90    .
cp $SRC/$MAKE .

# Compile SPEEDY and delete source files
make -f $MAKE -s clean
echo 'Compiling SPEEDY'

make -f $MAKE -s && { echo "Compilation successful"; } || { echo "Compilation failed"; exit 1; }
rm *.f90 *.o $MAKE *.mod
