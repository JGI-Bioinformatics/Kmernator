#!/bin/bash

sourceDir=$1
installDir=$2
if [ -z "$sourceDir" ] || [ ! -d "$sourceDir" ] || [ -z "$installDir" ]
then
  echo "USAGE: $0 /path/to/KmernatorSource /path/to/install/version [moduletag] [BuildType = Release]
where BuildType can be Release or Debug"
  exit 1
fi
moduletag=$3
if [ -z "$moduletag" ]
then
  moduletag=dev
fi
buildType=$4

set -e
buildDir=build
buildDirective=
installSuffix=
if [ -n "$buildType" ]
then
  buildDirective="-DCMAKE_BUILD_TYPE=$buildType"
  buildDir=${buildDir}_$buildType
  installSuffix=_${buildType}
fi

mkdir -p $buildDir
cd $buildDir

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

  echo "Executing cmake (see `pwd`/cmake.log) `date`" | tee -a cmake.log
  cmake $sourceDir -DCMAKE_C_COMPILER=`which cc` -DCMAKE_CXX_COMPILER=`which CC` -DMPI_COMPILER=`which CC` \
    $buildDirective \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    -DKMERNATOR_C_FLAGS="$KMERNATOR_C_FLAGS" -DKMERNATOR_CXX_FLAGS="$KMERNATOR_C_FLAGS" \
    -DCMAKE_SHARED_LIBRARY_LINK_C_FLAGS="" -DCMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS="" \
    -DMPI_LIBRARY="$mpilib2" -DMPI_INCLUDE_PATH="${MPICH_DIR}/include" \
    -DBOOST_MPI_USER_CONFIG="using mpi : : <find-static-library>$mpilib ;" -DMPI_EXTRA_LIBRARY="" \
    2>&1 | tee -a cmake.log

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
  module load samtools

  module list

  echo "Executing cmake $sourceDir (see `pwd`/cmake.log) `date`" | tee -a cmake.log
  cmake $sourceDir  \
    $buildDirective \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    2>&1 | tee -a cmake.log 

fi



echo "Building executables (see `pwd`/make.log) `date`" | tee -a make.log
make VERBOSE=1 -j18 2>&1 | tee -a make.log 
echo "Testing (see `pwd`/make-test.log) `date`" | tee -a make-test.log
make test 2>&1 | tee -a make-test.log
echo "Installing (see `pwd`/make-install.log) `date`" | tee -a make-install.log
make install 2>&1 | tee -a  make-install.log 

ver=$(cat git-version)

instdir=${installDir}/${ver}${installSuffix}
echo "Installing to $instdir"

mkdir -p $instdir
rsync -a include bin $instdir/

echo "Set $ver as default $moduletag module"
echo $ver > $installDir/$moduletag

