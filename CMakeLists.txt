# $Header: /repository/PI_annex/robsandbox/KoMer/CMakeLists.txt,v 1.10 2009-10-28 18:42:57 regan Exp $
#
# KoMer/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)
enable_testing()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_BUILD_TYPE Debug)
enable_language(CXX)
set(CMAKE_CXX_FLAGS -g)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_SHARED_LIBS OFF)

set(Boost_USE_MULTITHREADED ON)
set(Booost_USE_SINGLETHREADED OFF)

find_package ( Boost 1.40.0 COMPONENTS regex mpi)
if(Boost_FOUND)
  include_directories(${BOOST_INCLUDE_DIRS})
  link_directories(${BOOST_LIB_DIRS})
else()
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2)
    message("Unpacking Distro Version of Boost")
    execute_process(COMMAND tar -xzf ${CMAKE_CURRENT_SOURCE_DIR}/opt/boost-1.40.0.cmake2.tar.gz )
  endif()
  add_subdirectory( ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 EXCLUDE_FROM_ALL)
  add_custom_target( boost_subset ALL)
  add_dependencies( boost_subset boost_mpi-mt-static 
                                 boost_regex-mt-static 
                                 boost_test_exec_monitor-mt-static
                                 boost_unit_test_framework-mt-static)
  include_directories( ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 )
  link_directories( ${CMAKE_CURRENT_BINARY_DIR}/lib )
endif()

include_directories( src )
add_subdirectory( src )
link_directories( src )

add_subdirectory( test )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/10.fastq )


#
# $Log: CMakeLists.txt,v $
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
