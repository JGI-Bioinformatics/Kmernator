# $Header: /repository/PI_annex/robsandbox/KoMer/CMakeLists.txt,v 1.2 2009-10-20 17:25:49 regan Exp $
#
# KoMer/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)

find_package ( Boost 1.40.0 )
if(Boost_FOUND)
  include_directories(${BOOST_INCLUDE_DIRS})
else()
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2)
    execute_process(COMMAND tar -xzf ${CMAKE_CURRENT_SOURCE_DIR}/opt/boost-1.40.0.cmake2.tar.gz )
  endif()
  ADD_SUBDIRECTORY( ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 ${CMAKE_CURRENT_BINARY_DIR}/boost-1.40.0.cmake2 )
endif()
set(Boost_USE_STATIC_LIBS ON)

add_subdirectory( src )

add_subdirectory( test )

#
# $Log: CMakeLists.txt,v $
# Revision 1.2  2009-10-20 17:25:49  regan
# added CVS tags
#
#
