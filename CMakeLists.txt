# $Header: /repository/PI_annex/robsandbox/KoMer/CMakeLists.txt,v 1.4 2009-10-21 00:02:00 cfurman Exp $
#
# KoMer/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_SHARED_LIBS OFF)

set(Boost_USE_MULTITHREADED ON)
set(Booost_USE_SINGLETHREADED OFF)

find_package ( Boost 1.40.0 COMPONENTS regex mpi)
if(Boost_FOUND)
  include_directories(${BOOST_INCLUDE_DIRS})
else()
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2)
    message("Unpacking Distro Version of Boost")
    execute_process(COMMAND tar -xzf ${CMAKE_CURRENT_SOURCE_DIR}/opt/boost-1.40.0.cmake2.tar.gz )
  endif()
  add_subdirectory( ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 EXCLUDE_FROM_ALL)
  add_custom_target( boost_subset ALL)
  add_dependencies( boost_subset boost_mpi-mt-static boost_regex-mt-static )
endif()


add_subdirectory( src )

add_subdirectory( test )
enable_testing()
add_test(my_test ${CMAKE_CURRENT_BINARY_DIR}/test/ktest2 ${CMAKE_CURRENT_SOURCE_DIR}/test/10.fastq)

#
# $Log: CMakeLists.txt,v $
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
