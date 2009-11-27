# $Header: /repository/PI_annex/robsandbox/KoMer/apps/CMakeLists.txt,v 1.1 2009-11-27 23:28:07 cfurman Exp $
#
 
cmake_minimum_required(VERSION 2.6)
 

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

 
add_executable( CompareSpectrums CompareSpectrums )
target_link_libraries( CompareSpectrums Sequence Kmer
                              ReadSet
                              boost_program_options-mt-static )


#
# $Log: CMakeLists.txt,v $
# Revision 1.1  2009-11-27 23:28:07  cfurman
# CompareSpectrum application added
#
 
