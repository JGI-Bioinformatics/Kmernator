#
# Kmernator/CMakeLists.txt
#
#
# Kmernator Copyright (c) 2012, The Regents of the University of California, 
# through Lawrence Berkeley National Laboratory (subject to receipt of any 
# required approvals from the U.S. Dept. of Energy).  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# (1) Redistributions of source code must retain the above copyright notice, this 
# list of conditions and the following disclaimer.
# 
# (2) Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#  
# (3) Neither the name of the University of California, Lawrence Berkeley 
# National Laboratory, U.S. Dept. of Energy nor the names of its contributors may 
# be used to endorse or promote products derived from this software without 
# specific prior written permission.
#  
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
# You are under no obligation whatsoever to provide any bug fixes, patches, or 
# upgrades to the features, functionality or performance of the source code 
# ("Enhancements") to anyone; however, if you choose to make your Enhancements 
# available either publicly, or directly to Lawrence Berkeley National 
# Laboratory, without imposing a separate written license agreement for such 
# Enhancements, then you hereby grant the following license: a  non-exclusive, 
# royalty-free perpetual license to install, use, modify, prepare derivative 
# works, incorporate into other computer software, distribute, and sublicense 
# such enhancements or derivative works thereof, in binary and source code form.
# 


set(CMAKE_USER_MAKE_RULES_OVERRIDE override.cmake)
# by default, build the RELEASE version
if(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release.")
else()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Building the default, release version")
endif()

project( Kmernator CXX C)


cmake_minimum_required(VERSION 2.8)
enable_testing()

set(Kmernator_VERSION_MAJOR "1")
set(Kmernator_VERSION_MINOR "0")
set(Kmernator_VERSION_TAG "0")

set(Kmernator_VERSION_FILE_TEMPLATE ${CMAKE_SOURCE_DIR}/src/version.h.in)
set(Kmernator_VERSION_FILE ${CMAKE_SOURCE_DIR}/src/version.h)
set(Kmernator_GIT_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/git-version)

function(git_describe _var)
  find_package(Git QUIET)
  execute_process(COMMAND "${GIT_EXECUTABLE}" describe --tags
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  RESULT_VARIABLE res
                  OUTPUT_VARIABLE out
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message("Found Git version: ${out}")
  set(${_var} "${out}" PARENT_SCOPE)
endfunction()

if( EXISTS "${CMAKE_SOURCE_DIR}/.git" )
  #git_describe(Kmernator_VERSION_TAG)
  add_custom_target(GET_GIT_VERSION ALL
                    COMMAND git describe --tags > ${Kmernator_GIT_VERSION_FILE}
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    COMMENT "Retrieving git version")
#  file(READ ${Kmernator_GIT_VERSION_FILE} Kmernator_VERSION_TAG)

else()
  add_custom_target(GET_GIT_VERSION COMMAND /bin/true)
  file(WRITE ${Kmernator_GIT_VERSION_FILE} "${Kmernator_VERSION_TAG}")
endif()

  add_custom_target(REPLACE_VERSION_H ALL
                    COMMAND ${CMAKE_COMMAND}
			-DKmernator_GIT_VERSION_FILE=${Kmernator_GIT_VERSION_FILE} 
			-DKmernator_VERSION_FILE_TEMPLATE=${Kmernator_VERSION_FILE_TEMPLATE} 
			-DKmernator_VERSION_FILE=${Kmernator_VERSION_FILE} 
			-DKmernator_VERSION_MAJOR=${Kmernator_VERSION_MAJOR} 
			-DKmernator_VERSION_MINOR=${Kmernator_VERSION_MINOR} 
			-DKmernator_VERSION_TAG=${Kmernator_VERSION_TAG} 
			-P ${CMAKE_SOURCE_DIR}/makeVersionFile.cmake
                    COMMENT "re-building version.h")
  add_dependencies(REPLACE_VERSION_H GET_GIT_VERSION)

#configure_file(${Kmernator_VERSION_FILE_TEMPLATE} ${Kmernator_VERSION_FILE})

#set(Kmernator_VERSION "${Kmernator_VERSION_MAJOR}.${Kmernator_VERSION_MINOR}.${Kmernator_VERSION_TAG}")
#message("Building ${PROJECT_NAME} version ${Kmernator_VERSION}")
message("BuildType: ${CMAKE_BUILD_TYPE}")

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_SHARED_LIBS OFF)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_SINGLETHREADED OFF)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(BIT_CXX_FLAGS "-DIS_64_BIT")
endif()

#
# zlib
#
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
  include_directories(${ZLIB_INCLUDE_DIR})
  link_directories(${ZLIB_LIBRARIES})
else()
  error("Can not compile without zlib")
endif()

#
# OpenMP
#

if(DISABLE_OPENMP)
  set(OpenMP_CXX_FLAGS "-DDISABLE_OPENMP")
else()
  find_package( OpenMP )
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
if (DISABLE_MPI)
  set(KMERNATOR_MPI_CXX_FLAGS = "-DDISABLE_MPI")
else()
  find_package( MPI )
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

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package ( Boost 1.50.0 COMPONENTS mpi program_options iostreams thread serialization system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  include_directories(${BOOST_ROOT}/libs)
  link_directories(${Boost_LIBRARY_DIRS})
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
      file(DOWNLOAD ${BOOST_DOWNLOAD_URL} ${BOOST_TBZ2} STATUS DL_STATUS LOG DL_LOG)
      list(GET STATUS_DL 0 STATUS_DL0)
      if((NOT EXISTS ${BOOST_TBZ2}) OR (NOT STATUS_DL0 EQUAL 0))
        execute_process(COMMAND wget -O ${BOOST_TBZ2} ${BOOST_DOWNLOAD_URL} RESULT_VARIABLE DL_STATUS0)
        if((NOT EXISTS ${BOOST_TBZ2}) OR (NOT DL_STATUS0 EQUAL 0))
          message(FATAL_ERROR "Could not download boost from ${BOOST_DOWNLOAD_URL} (${DL_STATUS0}) : ${DL_LOG}")
        endif()
      endif()
    endif()
    message("Unpacking Boost locally: ${BOOST_INSTALL_VER}")
    execute_process(COMMAND tar -xjf ${BOOST_TBZ2} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    RESULT_VARIABLE SUCCESS)
    if(NOT SUCCESS EQUAL 0)
      message(FATAL_ERROR "Could not unpack boost")
    endif()

    if(NOT DISABLE_MPI)
      if(NOT BOOST_MPI_USER_CONFIG)
        file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "using mpi ;\n")
      else()
        file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "${BOOST_MPI_USER_CONFIG}\n")
      endif()
    endif()
    file(APPEND ${BOOST_ROOT}/tools/build/v2/user-config.jam "using gcc : : ${CMAKE_CXX_COMPILER} ;\n")

    execute_process(COMMAND ./bootstrap.sh --with-libraries=${BOOST_LIBRARY_LIST} WORKING_DIRECTORY ${BOOST_ROOT})
    execute_process(COMMAND ./b2 variant=release link=static threading=multi runtime-link=static WORKING_DIRECTORY ${BOOST_ROOT}) 
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
    file(DOWNLOAD ${SAMTOOLS_DOWNLOAD_URL} ${SAMTOOLS_TBZ2} STATUS DL_STATUS LOG DL_LOG)
    list(GET STATUS_DL 0 STATUS_DL0)
    if((NOT EXISTS ${SAMTOOLS_TBZ2}) OR (NOT STATUS_DL0 EQUAL 0))
      execute_process(COMMAND wget -O ${SAMTOOLS_TBZ2} ${SAMTOOLS_DOWNLOAD_URL} RESULT_VARIABLE DL_STATUS0)
      if((NOT EXISTS ${SAMTOOLS_TBZ2}) OR (NOT DL_STATUS0 EQUAL 0))
        message(FATAL_ERROR "Could not download samtools from ${SAMTOOLS_DOWNLOAD_URL}: ${DL_LOG}")
      endif()
    endif()
  endif()
  message("Unpacking samtools: ${SAMTOOLS_VERSION}")
  execute_process(COMMAND tar -xjf ${SAMTOOLS_TBZ2} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} RESULT_VARIABLE SUCCESS)
  if(NOT SUCCESS EQUAL 0)
    message(FATAL_ERROR "Could not unpack samtools")
  endif()
endif()

file(COPY ${CMAKE_SOURCE_DIR}/samtools_CMakeLists.txt DESTINATION ${SAMTOOLS_DIR})
file(RENAME ${SAMTOOLS_DIR}/samtools_CMakeLists.txt ${SAMTOOLS_DIR}/CMakeLists.txt)

include_directories( ${SAMTOOLS_DIR} )
add_subdirectory( ${SAMTOOLS_DIR} binary_dir )


set(MCT_VERSION1 1.6)
set(MCT_VERSION ${MCT_VERSION1}.0)
set(MCT_UNPACK_DIR ${CMAKE_CURRENT_BINARY_DIR}/mct-${MCT_VERSION})
set(MCT_DIR ${CMAKE_CURRENT_BINARY_DIR}/mct)
if(NOT EXISTS ${MCT_DIR})
  set(MCT_TGZ ${CMAKE_SOURCE_DIR}/opt/mct-${MCT_VERSION}.tar.gz)
  if(NOT EXISTS ${MCT_TGZ})
    set(MCT_DOWNLOAD_URL http://launchpad.net/libmct/${MCT_VERSION1}/${MCT_VERSION}/+download/mct-${MCT_VERSION}.tar.gz)
    message("Retrieving mct from repository: ${MCT_DOWNLOAD_URL}")
    file(DOWNLOAD ${MCT_DOWNLOAD_URL} ${MCT_TGZ} STATUS DL_STATUS LOG DL_LOG)
    list(GET STATUS_DL 0 STATUS_DL0)
    if((NOT EXISTS ${MCT_TGZ}) OR (NOT STATUS_DL0 EQUAL 0))
      execute_process(COMMAND wget -O ${MCT_TGZ} ${MCT_DOWNLOAD_URL} RESULT_VARIABLE DL_STATUS0)
      if((NOT EXISTS ${MCT_TGZ}) OR (NOT DL_STATUS0 EQUAL 0))
        message(FATAL_ERROR "Could not download mct from ${MCT_DOWNLOAD_URL}: ${DL_LOG}")
      endif()
    endif()
  endif()
  if(NOT EXISTS ${MCT_UNPACK_DIR})
    message("Unpacking mct: ${MCT_VERSION}")
    execute_process(COMMAND tar -xzf ${MCT_TGZ} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} RESULT_VARIABLE SUCCESS)
    if(NOT SUCCESS EQUAL 0)
      message(FATAL_ERROR "Could not unpack mct")
    endif()
  endif()
  execute_process(COMMAND make install prefix=${CMAKE_CURRENT_BINARY_DIR}/mct WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mct-${MCT_VERSION} RESULT_VARIABLE SUCCESS)
  if(NOT SUCCESS EQUAL 0)
    message(FATAL_ERROR "Could not install mct to ${CMAKE_CURRENT_BINARY_DIR}/mct")
  endif()
  
endif()

include_directories( ${MCT_DIR}/include )

include_directories( src )

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${BIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

add_subdirectory( src )

add_subdirectory( test )

add_subdirectory( apps )

