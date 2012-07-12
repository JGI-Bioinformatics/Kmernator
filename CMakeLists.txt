#
# Kmernator/CMakeLists.txt
#
# Copyright 2010 The Regents of the University of California.
# All rights reserved.
#
# The United States Government has rights in this work pursuant
# to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
# W-7405-ENG-48 between the United States Department of Energy
# and the University of California.
#
# Redistribution and use in source and binary forms are permitted
# provided that: (1) source distributions retain this entire
# copyright notice and comment, and (2) distributions including
# binaries display the following acknowledgement:  "This product
# includes software developed by the University of California,
# JGI-PSF and its contributors" in the documentation or other
# materials provided with the distribution and in all advertising
# materials mentioning features or use of this software.  Neither the
# name of the University nor the names of its contributors may be
# used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE.
#

project( Kmernator )
cmake_minimum_required(VERSION 2.6)
enable_testing()

if( EXISTS "${CMAKE_SOURCE_DIR}/.git" )

  add_custom_target(GET_GIT_VERSION ALL
                    COMMAND git describe --tags > ${CMAKE_CURRENT_BINARY_DIR}/git-version
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    COMMENT "Retrieving git version")
                    
  add_custom_target(REPLACE_GIT_VERSION ALL 
                    COMMAND ${CMAKE_SOURCE_DIR}/make_version.pl ${CMAKE_SOURCE_DIR}/version.h ${CMAKE_CURRENT_BINARY_DIR}/git-version 
                    COMMENT "re-building version.h")

endif()

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_SHARED_LIBS OFF)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_SINGLETHREADED OFF)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(BIT_CXX_FLAGS "-DIS_64_BIT")
endif()

#
# OpenMP
#
find_package( OpenMP )
if(DISABLE_OPENMP)
  set(OpenMP_CXX_FLAGS "-DDISABLE_OPENMP")
else()
  if(OPENMP_FOUND)
     set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DENABLE_OPENMP")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}") 
  else()
     set(OpenMP_CXX_FLAGS "-DDISABLE_OPENMP")
  endif()
endif()

#
# MPI
#
find_package( MPI )
if (DISABLE_MPI)
  set(KMERNATOR_MPI_CXX_FLAGS = "-DDISABLE_MPI")
else()
  if ( MPI_FOUND )
    set(KMERNATOR_MPI_CXX_FLAGS "${KMERNATOR_MPI_CXX_FLAGS} ${MPI_COMPILE_FLAGS} -DENABLE_MPI")
    set(KMERNATOR_MPI_LINK_FLAGS "${KMERNATOR_MPI_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    link_directories(${MPI_LINK_PATH})
    set(KMERNATOR_P_LIBS boost_mpi boost_serialization ${MPI_LIBRARIES})
  else()
      set(DISABLE_MPI)
      set(KMERNATOR_MPI_CXX_FLAGS "-DDISABLE_MPI")
  endif()
endif()

set(Boost_DEBUG 1)
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package ( Boost 1.41.0 COMPONENTS mpi program_options iostreams thread serialization system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  include_directories(${BOOST_ROOT}/libs)
  link_directories(${Boost_LIBRARY_DIRS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else()
  set(BOOST_MAJOR_VER 1)
  set(BOOST_MINOR_VER 50)
  set(BOOST_MINOR_SUBVER 0)
  set(BOOST_INSTALL_VER ${BOOST_MAJOR_VER}_${BOOST_MINOR_VER}_${BOOST_MINOR_SUBVER})
  set(BOOST_FULL_VERSION ${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}.${BOOST_MINOR_SUBVER})
  set(BOOST_DOWNLOAD_URL http://sourceforge.net/projects/boost/files/boost/${BOOST_FULL_VERSION}/boost_${BOOST_INSTALL_VER}.tar.bz2)
  set(BOOST_ROOT ${CMAKE_CURRENT_BINARY_DIR}/boost_${BOOST_INSTALL_VER})
  if(NOT EXISTS ${BOOST_ROOT})
    set(BOOST_LIBRARY_LIST iostreams,program_options,thread,test,system)
    if(NOT DISABLE_MPI)
      set(BOOST_LIBRARY_LIST mpi,serialization,${BOOST_LIBRARY_LIST})
    endif()

    set(BOOST_TBZ2 ${CMAKE_SOURCE_DIR}/opt/boost_${BOOST_INSTALL_VER}.tar.bz2)
    if(NOT EXISTS ${BOOST_TBZ2})
      message("Retrieving Boost from repository: ${BOOST_INSTALL_VER}")
      file(DOWNLOAD ${BOOST_DOWNLOAD_URL} ${BOOST_TBZ2})
    endif()
    message("Unpaking Boost: ${BOOST_INSTALL_VER}")
    execute_process(COMMAND tar -xvjf ${BOOST_TBZ2} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    if(NOT DISABLE_MPI)
      if(NOT BOOST_MPI_USER_CONFIG)
        file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "using mpi ;\n")
      else()
        file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "${BOOST_MPI_USER_CONFIG}\n")
      endif()
    endif()
    file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "using gcc : : ${CMAKE_CXX_COMPILER} ;\n")
  endif()

  find_package ( Boost ${BOOST_FULL_VERSION} COMPONENTS mpi program_options iostreams thread serialization system)
  if(NOT Boost_FOUND)
    execute_process(COMMAND ./bootstrap.sh --with-libraries=${BOOST_LIBRARY_LIST} WORKING_DIRECTORY ${BOOST_ROOT})
    execute_process(COMMAND ./b2 --debug-configuration variant=release link=static threading=multi runtime-link=static WORKING_DIRECTORY ${BOOST_ROOT}) 
  endif()

  message("Using local copy of boost: ${BOOST_ROOT}")
  include_directories( ${BOOST_ROOT} )
  link_directories( ${BOOST_ROOT}/stage/lib )

endif()

set(KMERNATOR_BOOST_LIBS         boost_program_options
                                 boost_iostreams
                                 boost_thread
				 boost_system
                                 )
set(KMERNATOR_BOOST_TEST_LIBS    boost_test_exec_monitor
                                 boost_unit_test_framework
                                 )
  
set(SAMTOOLS_VERSION 0.1.18)
set(SAMTOOLS_DIR ${CMAKE_CURRENT_BINARY_DIR}/samtools-${SAMTOOLS_VERSION})
if(NOT EXISTS ${SAMTOOLS_DIR})
  set(SAMTOOLS_TBZ2 ${CMAKE_SOURCE_DIR}/opt/samtools-${SAMTOOLS_VERSION}.tar.bz2)
  if(NOT EXISTS ${SAMTOOLS_TBZ2})
    set(SAMTOOLS_DOWNLOAD_URL http://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2)
    message("Retrieving samtools from repository: ${SAMTOOLS_DOWNLOAD_URL}")
    file(DOWNLOAD ${SAMTOOLS_DOWNLOAD_URL} ${SAMTOOLS_TBZ2})
  endif()
  message("Unpacking samtools: ${SAMTOOLS_VERSION}")
  execute_process(COMMAND tar -xvjf ${SAMTOOLS_TBZ2} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/samtools_CMakeLists.txt DESTINATION ${SAMTOOLS_DIR})
  file(RENAME ${SAMTOOLS_DIR}/samtools_CMakeLists.txt ${SAMTOOLS_DIR}/CMakeLists.txt)
endif()

include_directories( ${SAMTOOLS_DIR} )
add_subdirectory( ${SAMTOOLS_DIR} binary_dir )

include_directories( src )
add_subdirectory( src )

add_subdirectory( test )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/10.fastq )

add_subdirectory( apps )

