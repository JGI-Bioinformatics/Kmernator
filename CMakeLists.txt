# $Header: /repository/PI_annex/robsandbox/KoMer/CMakeLists.txt,v 1.20 2010-03-08 18:11:55 regan Exp $
#
# KoMer/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)
enable_testing()

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_SHARED_LIBS OFF)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_SINGLETHREADED OFF)

find_package( OpenMP )
if(OPENMP_FOUND)
  set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DENABLE_OPENMP")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}") 
else()
  set(OpenMP_CXX_FLAGS "-DDISABLE_OPENMP")
endif()

find_package ( Boost 1.41.0 COMPONENTS regex mpi)
if(Boost_FOUND)
  include_directories(${BOOST_INCLUDE_DIRS})
  link_directories(${BOOST_LIB_DIRS})
else()
  set(BOOST_INSTALL_VER boost-1.40.0.cmake2)
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_INSTALL_VER})
    message("Unpacking Distro Version of Boost")
    execute_process(COMMAND tar -xzf ${CMAKE_CURRENT_SOURCE_DIR}/opt/${BOOST_INSTALL_VER}.tar.gz )
  endif()
  add_subdirectory( ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_INSTALL_VER} ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_INSTALL_VER} EXCLUDE_FROM_ALL)
  add_custom_target( boost_subset ALL)
  add_dependencies( boost_subset boost_mpi-mt-static 
                                 boost_regex-mt-static 
                                 boost_test_exec_monitor-mt-static
                                 boost_unit_test_framework-mt-static
                                 boost_program_options-mt-static)
  include_directories( ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_INSTALL_VER} )
  link_directories( ${CMAKE_CURRENT_BINARY_DIR}/lib )
endif()

include_directories( src )
add_subdirectory( src )
link_directories( src )

add_subdirectory( test )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/10.fastq )

add_subdirectory( apps )
#
# $Log: CMakeLists.txt,v $
# Revision 1.20  2010-03-08 18:11:55  regan
# bugfix
#
# Revision 1.19  2010-03-04 06:37:04  regan
# fixed to support building without openmp
#
# Revision 1.18  2009-12-24 00:38:45  cfurman
# CompareSpectrum application added
#
# Revision 1.17  2009-12-18 19:05:07  regan
# added compile-time thread-safety option to kmer classes
#
# Revision 1.16  2009-11-26 09:03:32  regan
# refactored and stuff
#
# Revision 1.15  2009-11-10 18:28:09  regan
# *** empty log message ***
#
# Revision 1.13  2009-11-06 04:06:16  regan
# added boost options modules
#
# Revision 1.12  2009-10-31 23:38:51  regan
# fixed build to read env CXX_FLAGS
#
# Revision 1.11  2009-10-30 19:27:42  regan
# added iterator goodness, but KmerMap::Iterator still does not work
#
# Revision 1.10  2009-10-28 18:42:57  regan
# added debug flags, fixed tests, bugs
#
# Revision 1.9  2009-10-26 17:42:52  regan
# added debug build
#
# Revision 1.8  2009-10-23 21:14:25  regan
# fixed test data
#
# Revision 1.7  2009-10-23 07:06:55  regan
# more unit testing
#   ReadSetTest
#   KmerTest
#
# Revision 1.6  2009-10-22 07:04:01  regan
# added a few unit tests
# minor refactor
#
# Revision 1.5  2009-10-21 06:51:33  regan
# bug fixes
# build lookup tables for twobitsequence
#
# Revision 1.4  2009-10-21 00:02:00  cfurman
# working on kmers....
#
# Revision 1.3  2009-10-20 20:56:25  cfurman
# Got it to compile!
#
# Revision 1.2  2009-10-20 17:25:49  regan
# added CVS tags
#
#
