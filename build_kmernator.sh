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

SKIP_TEST=${SKIP_TEST:=0}

procs=$(($(lscpu -p 2>/dev/null | tail -1 | awk -F, '{print $2}')+1))

logcmd()
{
  logfile=$1
  shift
  echo "Starting $@ at $(date)" >> $logfile
  $@ 2>&1 | tee -a $logfile
  [ "${PIPESTATUS[*]}" == "0 0" ]
}

set -e
buildDir=${BUILD_DIR:=build}
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

execute_modules=""
build_modules=""
cmakeOptions=""

if [ "$NERSC_HOST" == "edison" ] || [ "${NERSC_HOST}" == "cori" ]
then
  # edison is Cray XC30

  execute_modules="PrgEnv-gnu"
  build_modules="cmake boost zlib bzip2"

  module rm PrgEnv-intel
  module rm PrgEnv-gnu
  module rm boost
  module rm cmake
  module rm bzip2
  module rm zlib

  module load $execute_modules $build_modules
  module list

  mpilib=$(basename $MPICH_DIR/lib/libmpich_gnu*.a); mpilib=${mpilib%.a} ; mpilib2=${mpilib#lib}
  export KMERNATOR_C_FLAGS="-Wno-unused-local-typedefs"

  set -x
  cmake -DCMAKE_C_COMPILER=`which cc` -DCMAKE_CXX_COMPILER=`which CC` -DMPI_COMPILER=`which CC` \
    $buildDirective \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    -DKMERNATOR_C_FLAGS="${KMERNATOR_C_FLAGS}" -DKMERNATOR_CXX_FLAGS="${KMERNATOR_C_FLAGS}" \
    -DCMAKE_SHARED_LIBRARY_LINK_C_FLAGS='' -DCMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS='' \
    -DMPI_LIBRARY="$mpilib2" -DMPI_INCLUDE_PATH="${MPICH_DIR}/include" \
    -DBOOST_MPI_USER_CONFIG="using mpi : : <find-static-library>$mpilib ;" -DMPI_EXTRA_LIBRARY='' \
    $sourceDir
  set +x

elif [ "$NERSC_HOST" == "genepool" ]
then
  # genepool is linux with openmpi

  module purge
  export GP_INFINIBAND=0
  execute_modules="PrgEnv-gnu openmpi boost"
  build_modules="cmake git samtools"

  module load $execute_modules $build_modules

  set -x
  cmake \
    $buildDirective \
    -DCMAKE_INSTALL_PREFIX:PATH=`pwd` \
    $sourceDir
  set +x

else

  set -x
  cmake $buildDirective -DCMAKE_INSTALL_PREFIX:PATH=`pwd` $sourceDir
  set +x

fi

> .deps
for module in ${execute_modules}
do
  
  echo $module >> .deps
done

module list || /bin/true

echo "Building executables (see `pwd`/make.log) `date`" | tee -a make.log
logcmd make.log make VERBOSE=1 -j$procs

if [ $SKIP_TEST -eq 0 ]
then
  OLD_TMP=$TMPDIR
  export TMPDIR=$(mktemp -d)
  echo "Testing (see `pwd`/make-test.log) `date`" | tee -a make-test.log
  logcmd make-test.log make test
  export TMPDIR=$OLD_TMP
fi

echo "Installing (see `pwd`/make-install.log) `date`" | tee -a make-install.log
logcmd make-install.log make install

ver=$(cat git-version)

instdir=${installDir}/${ver}${installSuffix}
echo "Installing to $instdir"

mkdir -p $instdir
rsync -a include bin $instdir/
cp .deps $instdir

echo "Set $ver as default $moduletag module"
echo $ver > $installDir/$moduletag

