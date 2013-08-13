#!/bin/bash

sourceDir=$1
installDir=$2
if [ -z "$sourceDir" ] || [ ! -d "$sourceDir" ] || [ -z "$installDir" ]
then
  echo "USAGE: $0 /path/to/KmernatorSource /path/to/install/version [moduletag]"
  exit 1
fi
mod=$3
if [ -z "$mod" ]
then
  mod=prod
fi

set -e

mkdir -p build
cd build

echo "Building in `pwd`"

unset CFLAGS
unset CXXFLAGS

if [ "$NERSC_HOST" == "edison" ]
then
  # edison is Cray XC30

  module rm PrgEnv-intel
  module rm boost
  module rm cmake
  module rm bzip2
  module rm cray-mpich
  module load PrgEnv-gnu
  module load cray-mpich
  module load cmake
  module load bzip2
  module load samtools

  module list

  mpilib=$(basename $MPICH_DIR/lib/libmpich_gnu*.a); mpilib=${mpilib%.a} ; mpilib2=${mpilib#lib}
  export KMERNATOR_C_FLAGS="-Wall -Wno-unused-local-typedefs"

  echo "Executing cmake (see `pwd`/cmake.log)"
  cmake $sourceDir -DCMAKE_C_COMPILER=`which cc` -DCMAKE_CXX_COMPILER=`which CC` -DMPI_COMPILER=`which CC` \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    -DKMERNATOR_C_FLAGS="$KMERNATOR_C_FLAGS" -DKMERNATOR_CXX_FLAGS="$KMERNATOR_C_FLAGS" \
    -DCMAKE_SHARED_LIBRARY_LINK_C_FLAGS="" -DCMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS="" \
    -DMPI_LIBRARY="$mpilib2" -DMPI_INCLUDE_PATH="${MPICH_DIR}/include" \
    -DBOOST_MPI_USER_CONFIG="using mpi : : <find-static-library>$mpilib ;" -DMPI_EXTRA_LIBRARY="" \
    > cmake.log 2>&1

fi

if [ "$NERSC_HOST" == "genepool" ]
then
  # genepool is linux with openmpi

  module purge
  export GP_INFINIBAND=0
  module load PrgEnv-gnu
  module load openmpi
  module load cmake
  module load git
  module load boost/1.53.0

  module list

  echo "Executing cmake $sourceDir (see `pwd`/cmake.log)"
  cmake $sourceDir  \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    > cmake.log 2>&1

fi



echo "Building executables (see `pwd`/make.log)"
make VERBOSE=1 -j18 >> make.log 2>&1 
echo "Testing (see `pwd`/make-test.log)"
make test >> make-test.log 2>&1 
echo "Installing (see `pwd`/make-install.log)"
make install >> make-install.log 2>&1

ver=$(cat git-version)

instdir=${installDir}/${ver}
echo "Installing to $instdir"

mkdir -p $instdir
rsync -a include bin $instdir/

echo "Set $ver as default $mod module"
echo $ver > $installDir/$mod

