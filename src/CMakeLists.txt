#
# Kmernator/src/CMakeLists.txt
#
# Copyright 2010 The Regents of the University of California.
# All rights reserved.
#
# Author: Rob Egan
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

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)

IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${BIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set_source_files_properties( TwoBitSequence
                            PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                          )
set_source_files_properties( Sequence
                            PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                          )
set_source_files_properties( ReadSet
                            PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                          )

set_source_files_properties( Kmer
                            PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                          )
set_source_files_properties( Log
                            PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                          )
                          
                           
add_library( TwoBitSequence TwoBitSequence )
add_library( Sequence Sequence TwoBitSequence )
add_library( Kmer Kmer TwoBitSequence)
add_library( ReadSet ReadSet )
add_library( Log Log )
add_dependencies( Log REPLACE_VERSION_H)

