#
# Kmernator/test/CMakeLists.txt
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

cmake_minimum_required(VERSION 2.6)
enable_testing()

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

set_source_files_properties( BgzfStreamTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( BgzfStreamTest BgzfStreamTest )
target_link_libraries( BgzfStreamTest
                                   Log
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   z
                                   )
add_test( BgzfStreamTest BgzfStreamTest )

if( MPI_FOUND )
add_test(testSplitSequenceOnTheFly runSplitSequenceOnTheFlyTests.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runSplitSequenceOnTheFlyTests.sh test/ )

add_executable( TestMPI TestMPI )
set_source_files_properties( TestMPI
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( TestMPI
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
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
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
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

add_executable( SamUtilsTest SamUtilsTest)
set_source_files_properties( SamUtilsTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( SamUtilsTest
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries (SamUtilsTest
                            Log
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              libbam
                              z
                              )
add_test(SamUtilsTest SamUtilsTest.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/SamUtilsTest.sh test/ )

endif()

execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.1.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.2.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fasta test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.qual test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus1.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2-diff.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus3.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength-both.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runFilterTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mercount.m21 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mergraph.m21.D2 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runMeraculousTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runRandomlySampleTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10k.bam test/ )

add_test(filter runFilterTests.sh)
add_test(meraculous runMeraculousTests.sh)
add_test(randomlysample runRandomlySampleTests.sh)
