#
# Kmernator/apps/CMakeLists.txt
#
#
# Kmernator Copyright (c) 2012, The Regents of the University of California, 
# through Lawrence Berkeley National Laboratory (subject to receipt of any 
# required approvals from the U.S. Dept. of Energy).  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# (1) Redistributions of source code must retain the above copyright notice, this 
# list of conditions and the following disclaimer.
# 
# (2) Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#  
# (3) Neither the name of the University of California, Lawrence Berkeley 
# National Laboratory, U.S. Dept. of Energy nor the names of its contributors may 
# be used to endorse or promote products derived from this software without 
# specific prior written permission.
#  
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
# You are under no obligation whatsoever to provide any bug fixes, patches, or 
# upgrades to the features, functionality or performance of the source code 
# ("Enhancements") to anyone; however, if you choose to make your Enhancements 
# available either publicly, or directly to Lawrence Berkeley National 
# Laboratory, without imposing a separate written license agreement for such 
# Enhancements, then you hereby grant the following license: a  non-exclusive, 
# royalty-free perpetual license to install, use, modify, prepare derivative 
# works, incorporate into other computer software, distribute, and sublicense 
# such enhancements or derivative works thereof, in binary and source code form.
# 
 
cmake_minimum_required(VERSION 2.8)
 
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
add_dependencies(CompareSpectrums REPLACE_VERSION_H)
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
add_dependencies(ContigExtender REPLACE_VERSION_H)
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
add_dependencies(FilterReads REPLACE_VERSION_H)
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
add_dependencies(FilterReads-P REPLACE_VERSION_H)
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
add_dependencies(EstimateSize-P REPLACE_VERSION_H)
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
add_dependencies(DistributedNucleatingAssembler REPLACE_VERSION_H)
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
add_dependencies(MeraculousCounter REPLACE_VERSION_H)
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
add_dependencies(BamSort-P REPLACE_VERSION_H)
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
add_dependencies(MPIBroadcastStream REPLACE_VERSION_H)
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
add_dependencies(TnfDistance REPLACE_VERSION_H)
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
add_dependencies(Fastq2Fasta REPLACE_VERSION_H)
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
add_dependencies(Fastq2FastaQual REPLACE_VERSION_H)
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
add_dependencies(FixPair REPLACE_VERSION_H)
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
add_dependencies(SplitSequenceOnTheFly REPLACE_VERSION_H)
install(TARGETS SplitSequenceOnTheFly DESTINATION bin)


endif()

add_executable( RandomlySample RandomlySample )
target_link_libraries(        RandomlySample
                              Sequence
                              ReadSet
                              Log
                              ${KMERNATOR_BOOST_LIBS}
                              )
add_dependencies(RandomlySample REPLACE_VERSION_H)
install(TARGETS RandomlySample DESTINATION bin)

add_executable( pipebuf pipebuf )
set_source_files_properties( pipebuf PROPERTIES COMPILE_FLAGS "-D_GNU_SOURCE=1" )
add_dependencies(pipebuf REPLACE_VERSION_H)
install(TARGETS pipebuf DESTINATION bin)

add_executable( FilterAndSortSam FilterAndSortSam )
set_source_files_properties(FilterAndSortSam
                        PROPERTIES COMPILE_FLAGS " ${CMAKE_C_FLAGS} ")
target_link_libraries( FilterAndSortSam libbam z)
add_dependencies(FilterAndSortSam REPLACE_VERSION_H)
install(TARGETS FilterAndSortSam DESTINATION bin)

