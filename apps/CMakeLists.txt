#
# Kmernator/apps/CMakeLists.txt
#
# Copyright 2010 The Regents of the University of California.
# All rights reserved.
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

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${CMAKE_CXX_FLAGS}")

set_source_files_properties( CompareSpectrums
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( CompareSpectrums CompareSpectrums )
target_link_libraries( CompareSpectrums 
                              Sequence 
                              Kmer
                              ReadSet
                              Log
                              boost_program_options-mt-static 
                              boost_iostreams-mt-static
                              boost_thread-mt-static
                              )


set_source_files_properties( FilterReads
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( FilterReads FilterReads )
target_link_libraries( FilterReads 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )

if( MPI_FOUND )

set_source_files_properties( FilterReads-P
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( FilterReads-P FilterReads-P )
set_target_properties( FilterReads-P
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}"
                       )
target_link_libraries( FilterReads-P 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
endif()

set_source_files_properties( TnfDistance
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( TnfDistance TnfDistance )
target_link_libraries( TnfDistance 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
   
set_source_files_properties( Fastq2Fasta
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( Fastq2Fasta Fastq2Fasta )
target_link_libraries(        Fastq2Fasta
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )

set_source_files_properties( Fastq2FastaQual
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( Fastq2FastaQual Fastq2FastaQual )
target_link_libraries(        Fastq2FastaQual
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )


set_source_files_properties( FixPair
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( FixPair FixPair )
target_link_libraries(        FixPair
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )


add_executable( SplitSequenceOnTheFly SplitSequenceOnTheFly )
if( MPI_FOUND )
set_source_files_properties( SplitSequenceOnTheFly
                             PROPERTIES COMPILE_FLAGS "${KMERNATOR_MPI_CXX_FLAGS}"
                             )

target_link_libraries(        SplitSequenceOnTheFly
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
else()
target_link_libraries(        SplitSequenceOnTheFly
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )

endif()

add_executable( RandomlySample RandomlySample )
target_link_libraries(        RandomlySample
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )

#
# $Log: CMakeLists.txt,v $
# Revision 1.9  2010-08-18 17:47:49  regan
# added TNF distance application
#
# Revision 1.8.16.1  2010-07-13 17:39:12  regan
# created TnfDistance application
#
# Revision 1.8  2010-05-05 06:28:40  regan
# merged changes from FixPairOutput-20100504
#
# Revision 1.7.4.1  2010-05-04 21:33:44  regan
# checkpoint
#
# Revision 1.7  2010-05-01 21:58:00  regan
# merged head with serial threaded build partitioning
#
# Revision 1.6.4.1  2010-04-30 23:33:42  regan
# replaced read cache with LRU 3rd party cache
#
# Revision 1.6  2010-04-16 22:44:23  regan
# merged HEAD with changes for mmap and intrusive pointer
#
# Revision 1.4.2.2  2010-04-05 05:42:56  regan
# checkpoint mmaping input files
#
# Revision 1.4.2.1  2010-04-05 03:32:13  regan
# moved read file reader
#
# Revision 1.4  2010-03-16 18:58:10  regan
# created a simple fastq2fasta program
#
# Revision 1.3  2010-01-13 23:32:35  regan
# new app
#
# Revision 1.2  2010-01-08 18:34:51  regan
# refactored a bit
# enabled openmp parallelizations
#
# Revision 1.1  2009-11-27 23:28:07  cfurman
# CompareSpectrum application added
#
 
