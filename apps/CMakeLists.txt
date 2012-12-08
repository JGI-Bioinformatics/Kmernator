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

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${BIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set_source_files_properties( CompareSpectrums
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( CompareSpectrums CompareSpectrums )
target_link_libraries( CompareSpectrums 
                              Sequence 
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
install(TARGETS CompareSpectrums DESTINATION bin)

set_source_files_properties( ContigExtender
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( ContigExtender ContigExtender )
target_link_libraries( ContigExtender 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
install(TARGETS ContigExtender DESTINATION bin)


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
install(TARGETS FilterReads DESTINATION bin)

if( MPI_FOUND )

set_source_files_properties( FilterReads-P
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( FilterReads-P FilterReads-P )
set_target_properties( FilterReads-P
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( FilterReads-P 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS FilterReads-P DESTINATION bin)

set_source_files_properties( EstimateSize-P
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( EstimateSize-P EstimateSize-P )
set_target_properties( EstimateSize-P
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( EstimateSize-P 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS EstimateSize-P DESTINATION bin)
install(FILES EstimateSize.R DESTINATION bin)

set_source_files_properties( DistributedNucleatingAssembler
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( DistributedNucleatingAssembler DistributedNucleatingAssembler )
set_target_properties( DistributedNucleatingAssembler
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( DistributedNucleatingAssembler 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS DistributedNucleatingAssembler DESTINATION bin)

set_source_files_properties( MeraculousCounter
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( MeraculousCounter MeraculousCounter )
set_target_properties( MeraculousCounter
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( MeraculousCounter 
                              Sequence
                              Kmer
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS MeraculousCounter DESTINATION bin)

set_source_files_properties( BamSort-P
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( BamSort-P BamSort-P )
set_target_properties( BamSort-P
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( BamSort-P 
                              libbam
                              ${ZLIB_LIBRARIES} 
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS BamSort-P DESTINATION bin)

set_source_files_properties( MPIBroadcastStream
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
add_executable( MPIBroadcastStream MPIBroadcastStream )
set_target_properties( MPIBroadcastStream
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries( MPIBroadcastStream 
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS MPIBroadcastStream DESTINATION bin)

                              
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
install(TARGETS TnfDistance DESTINATION bin)
   
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
install(TARGETS Fastq2Fasta DESTINATION bin)

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
install(TARGETS Fastq2FastaQual DESTINATION bin)


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
install(TARGETS FixPair DESTINATION bin)



if( MPI_FOUND )
add_executable( SplitSequenceOnTheFly SplitSequenceOnTheFly )

set_source_files_properties( SplitSequenceOnTheFly
                             PROPERTIES COMPILE_FLAGS "${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( SplitSequenceOnTheFly
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )

target_link_libraries(        SplitSequenceOnTheFly
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
install(TARGETS SplitSequenceOnTheFly DESTINATION bin)


endif()

add_executable( RandomlySample RandomlySample )
target_link_libraries(        RandomlySample
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
install(TARGETS RandomlySample DESTINATION bin)

add_executable( pipebuf pipebuf )
set_source_files_properties( pipebuf PROPERTIES COMPILE_FLAGS "-s -Wall -O6 -fexpensive-optimizations -fomit-frame-pointer -D_GNU_SOURCE=1" )
install(TARGETS pipebuf DESTINATION bin)

add_executable( FilterAndSortSam FilterAndSortSam )
set_source_files_properties(FilterAndSortSam
                        PROPERTIES COMPILE_FLAGS " ${CMAKE_C_FLAGS} ")
target_link_libraries( FilterAndSortSam libbam z)
install(TARGETS FilterAndSortSam DESTINATION bin)

