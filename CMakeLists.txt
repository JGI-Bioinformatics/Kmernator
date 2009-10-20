# KoMer CmakeLists.txt
PROJECT( KoMer )
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

ADD_SUBDIRECTORY( src )

ADD_SUBDIRECTORY( test )

