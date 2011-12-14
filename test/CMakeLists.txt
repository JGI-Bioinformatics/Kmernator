#
# Kmernator/test/CMakeLists.txt
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
enable_testing()

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${BIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set_source_files_properties( TwoBitSequenceTest
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( TwoBitSequenceTest TwoBitSequenceTest )
target_link_libraries( TwoBitSequenceTest TwoBitSequence
                                          Log
                                          ${KMERNATOR_BOOST_LIBS}
                                          ${KMERNATOR_BOOST_TEST_LIBS}
                                          )
add_test( TwoBitSequenceTest TwoBitSequenceTest )

set_source_files_properties( ReadSetTest
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
                            
add_executable( ReadSetTest ReadSetTest )
target_link_libraries( ReadSetTest Sequence
                                   ReadSet 
                                   Kmer
                                   Log
                                   ${KMERNATOR_BOOST_LIBS}
                                          ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_test( ReadSetTest ReadSetTest )
                            
set_source_files_properties( KmerTest
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )
                                                                                                           
add_executable( KmerTest KmerTest )
target_link_libraries( KmerTest TwoBitSequence Kmer Sequence ReadSet Log
                                ${KMERNATOR_BOOST_LIBS}
                                ${KMERNATOR_BOOST_TEST_LIBS}
                                )
add_test( KmerTest KmerTest )

set_source_files_properties( ktest2
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )

add_executable( ktest2 ktest2 )
target_link_libraries( ktest2 Sequence Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
add_test(ktest2 ktest2 --mask-simple-repeats 0 21 10.fastq)

set_source_files_properties( HashTester
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( HashTester HashTester)
target_link_libraries( HashTester 
                              Sequence 
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )

set_source_files_properties( MmapTempFileTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( MmapTempFileTest MmapTempFileTest )
target_link_libraries( MmapTempFileTest
                                   Log
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_test( MmapTempFileTest MmapTempFileTest )

set_source_files_properties( UtilTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( UtilTest UtilTest )
target_link_libraries( UtilTest
                                   Log
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_test( UtilTest UtilTest )

if( MPI_FOUND )
add_executable( TestMPI TestMPI )
set_source_files_properties( TestMPI
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( TestMPI
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}"
                       )
target_link_libraries (TestMPI
                            Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
add_test(testmpi testMPI.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/testMPI.sh test/ )
    
add_executable( TestKmerMatchMPI TestKmerMatchMPI)
set_source_files_properties( TestKmerMatchMPI
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( TestKmerMatchMPI
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}"
                       )
target_link_libraries (TestKmerMatchMPI
                            Log
                            Kmer
                            Sequence
                            ReadSet
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
add_test(testKmerMatchMPI testKmerMatchMPI.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/testKmerMatchMPI.sh test/ )

endif()

execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/10.fastq )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fasta test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.qual test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus1.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2-diff.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus3.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runFilterTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mercount.m21 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mergraph.m21.D2 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runMeraculousTests.sh test/ )

add_test(filter runFilterTests.sh)
add_test(meraculous runMeraculousTests.sh)

#
# $Log: CMakeLists.txt,v $
# Revision 1.22  2010-05-18 20:50:21  regan
# merged changes from PerformanceTuning-20100506
#
# Revision 1.21.12.2  2010-05-12 18:25:27  regan
# added tests
#
# Revision 1.21.12.1  2010-05-07 22:59:41  regan
# refactored base type declarations
#
# Revision 1.21  2010-05-01 21:57:51  regan
# merged head with serial threaded build partitioning
#
# Revision 1.20.2.3  2010-04-30 23:33:44  regan
# replaced read cache with LRU 3rd party cache
#
# Revision 1.20.2.2  2010-04-26 22:52:44  regan
# more testing
#
# Revision 1.20.2.1  2010-04-23 23:39:23  regan
# added a quick program to test hashes
#
# Revision 1.20  2010-04-21 23:39:39  regan
# got kmermap mmap store and restore working
#
# Revision 1.19  2010-04-16 22:44:25  regan
# merged HEAD with changes for mmap and intrusive pointer
#
# Revision 1.18.2.2.2.1  2010-04-16 05:30:02  regan
# checkpoint.. broke it
#
# Revision 1.18.2.2  2010-04-05 05:42:57  regan
# checkpoint mmaping input files
#
# Revision 1.18.2.1  2010-04-05 03:32:14  regan
# moved read file reader
#
# Revision 1.18  2010-03-10 13:18:03  regan
# fixed linking
#
# Revision 1.17  2010-02-22 14:41:18  regan
# minor bug fixes
#
# Revision 1.16  2009-12-18 19:05:12  regan
# added compile-time thread-safety option to kmer classes
#
# Revision 1.15  2009-11-07 00:28:38  cfurman
# ReadSet now takes fasta, fastq or  fasta+qual files.
#
# Revision 1.14  2009-11-06 04:06:18  regan
# added boost options modules
#
# Revision 1.13  2009-11-02 18:48:44  regan
# fixed linking issues
#
# Revision 1.12  2009-10-31 23:41:48  regan
# fixed build to read env CXX_FLAGS
#
# Revision 1.11  2009-10-30 00:51:37  regan
# bug fix and working on executable
#
# Revision 1.10  2009-10-28 18:43:02  regan
# added debug flags, fixed tests, bugs
#
# Revision 1.9  2009-10-28 00:04:28  regan
# added more bugs
#
# Revision 1.8  2009-10-23 20:32:52  cfurman
# more kmer changes
#
# Revision 1.7  2009-10-23 07:06:57  regan
# more unit testing
#   ReadSetTest
#   KmerTest
#
# Revision 1.6  2009-10-23 01:24:56  cfurman
# ReadSet test created
#
# Revision 1.5  2009-10-22 07:04:03  regan
# added a few unit tests
# minor refactor
#
# Revision 1.4  2009-10-21 06:51:37  regan
# bug fixes
# build lookup tables for twobitsequence
#
# Revision 1.3  2009-10-20 20:56:29  cfurman
# Got it to compile!
#
# Revision 1.2  2009-10-20 17:25:53  regan
# added CVS tags
#
#
