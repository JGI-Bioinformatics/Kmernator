#
# Kmernator/src/CMakeLists.txt
#
# Copyright 2010 The Regents of the University of California.
# All rights reserved.
#
# Author: Rob Egan, Craig Furman
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

cmake_minimum_required(VERSION 2.6)

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

#
# $Log: CMakeLists.txt,v $
# Revision 1.11  2010-04-16 22:44:18  regan
# merged HEAD with changes for mmap and intrusive pointer
#
# Revision 1.10.2.1.2.1  2010-04-16 05:30:00  regan
# checkpoint.. broke it
#
# Revision 1.10.2.1  2010-04-14 21:29:38  regan
# fixed openmp bug
#
# Revision 1.10  2009-12-21 07:54:18  regan
# minor parallelization of reading files step
#
# Revision 1.9  2009-12-18 19:05:09  regan
# added compile-time thread-safety option to kmer classes
#
# Revision 1.8  2009-11-26 09:03:30  regan
# refactored and stuff
#
# Revision 1.7  2009-10-31 23:41:46  regan
# fixed build to read env CXX_FLAGS
#
# Revision 1.6  2009-10-28 18:42:59  regan
# added debug flags, fixed tests, bugs
#
# Revision 1.5  2009-10-23 20:32:50  cfurman
# more kmer changes
#
# Revision 1.4  2009-10-22 07:04:06  regan
# added a few unit tests
# minor refactor
#
# Revision 1.3  2009-10-21 06:51:34  regan
# bug fixes
# build lookup tables for twobitsequence
#
# Revision 1.2  2009-10-20 17:25:50  regan
# added CVS tags
#
#
