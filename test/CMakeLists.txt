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
                                          ${KMERNATOR_BOOST_LIBS}
                                          ${KMERNATOR_BOOST_TEST_LIBS}
                                          )
add_dependencies(TwoBitSequenceTest REPLACE_VERSION_H)
add_test( TwoBitSequenceTest TwoBitSequenceTest )

set_source_files_properties( ReadSetTest
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
                            
add_executable( ReadSetTest ReadSetTest )
target_link_libraries( ReadSetTest Sequence
                                   ReadSet 
                                   Kmer
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_dependencies(ReadSetTest REPLACE_VERSION_H)
add_test( ReadSetTest ReadSetTest )
  

set_source_files_properties( ReadSetStdTest
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
                            
add_executable( ReadSetStdTest ReadSetStdTest )
target_link_libraries( ReadSetStdTest Sequence
                                   ReadSet 
                                   Kmer
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_dependencies(ReadSetStdTest REPLACE_VERSION_H)
add_test( ReadSetStdTest ReadSetStdTest )
                          
set_source_files_properties( KmerTest
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )
                                                                                                           
add_executable( KmerTest KmerTest )
target_link_libraries( KmerTest TwoBitSequence Kmer Sequence ReadSet
                                ${KMERNATOR_BOOST_LIBS}
                                ${KMERNATOR_BOOST_TEST_LIBS}
                                )
add_dependencies(KmerTest REPLACE_VERSION_H)
add_test( KmerTest KmerTest )

set_source_files_properties( ktest2
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )

add_executable( ktest2 ktest2 )
target_link_libraries( ktest2 Sequence Kmer
                              ReadSet
                              ${KMERNATOR_BOOST_LIBS}
                              )
add_dependencies(ktest2 REPLACE_VERSION_H)
add_test(ktest2 ktest2 --mask-simple-repeats 0 21 10.fastq)

set_source_files_properties( HashTester
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( HashTester HashTester)
target_link_libraries( HashTester 
                              Sequence 
                              Kmer
                              ReadSet
                              ${KMERNATOR_BOOST_LIBS}
                              )
add_dependencies(HashTester REPLACE_VERSION_H)

set_source_files_properties( MmapTempFileTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( MmapTempFileTest MmapTempFileTest )
target_link_libraries( MmapTempFileTest
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_dependencies(MmapTempFileTest REPLACE_VERSION_H)
add_test( MmapTempFileTest MmapTempFileTest )

set_source_files_properties( UtilTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( UtilTest UtilTest )
target_link_libraries( UtilTest
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   )
add_dependencies(UtilTest REPLACE_VERSION_H)
add_test( UtilTest UtilTest )

set_source_files_properties( BgzfStreamTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
                             )
add_executable( BgzfStreamTest BgzfStreamTest )
target_link_libraries( BgzfStreamTest
                                   ${KMERNATOR_BOOST_LIBS}
                                   ${KMERNATOR_BOOST_TEST_LIBS}
                                   z
                                   )
add_dependencies(BgzfStreamTest REPLACE_VERSION_H)
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
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
add_dependencies(TestMPI REPLACE_VERSION_H)
add_test(testmpi testMPI.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/testMPI.sh test/ )
install(TARGETS TestMPI DESTINATION bin)
  
add_executable( TestForkDaemonMPI TestForkDaemonMPI )
set_source_files_properties( TestForkDaemonMPI
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( TestForkDaemonMPI
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries (TestForkDaemonMPI
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
add_dependencies(TestForkDaemonMPI REPLACE_VERSION_H)
add_test(TestForkDaemonMPI TestForkDaemonMPI.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/TestForkDaemonMPI.sh test/ )
install(TARGETS TestForkDaemonMPI DESTINATION bin)
    
add_executable( TestKmerMatchMPI TestKmerMatchMPI)
set_source_files_properties( TestKmerMatchMPI
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( TestKmerMatchMPI
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries (TestKmerMatchMPI
                            Kmer
                            Sequence
                            ReadSet
                              ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              )
add_dependencies(TestKmerMatchMPI REPLACE_VERSION_H)
add_test(testKmerMatchMPI testKmerMatchMPI.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/testKmerMatchMPI.sh test/ )
install(TARGETS TestKmerMatchMPI DESTINATION bin)

add_executable( SamUtilsTest SamUtilsTest)
set_source_files_properties( SamUtilsTest
                             PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${KMERNATOR_MPI_CXX_FLAGS}"
                             )
set_target_properties( SamUtilsTest
                       PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KMERNATOR_MPI_LINK_FLAGS}"
                       )
target_link_libraries (SamUtilsTest
                             ${KMERNATOR_BOOST_LIBS}
                              ${KMERNATOR_P_LIBS}
                              libbam
                              z
                              )
add_dependencies(SamUtilsTest REPLACE_VERSION_H)
add_test(SamUtilsTest SamUtilsTest.sh)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/SamUtilsTest.sh test/ )
install(TARGETS SamUtilsTest DESTINATION bin)

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
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-0.85.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength-both.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runFilterTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mercount.m21 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/phix.mergraph.m21.D2 test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runMeraculousTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/runRandomlySampleTests.sh test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10k.bam test/ )

execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.1.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10-cs18.2.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fasta test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.qual test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus1.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus2-diff.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/consensus3.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-0.85.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength.std.fastq test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/1000-Filtered-readlength-both.std.fastq test/ )

add_test(filter runFilterTests.sh)

add_test(meraculous runMeraculousTests.sh)

add_test(randomlysample runRandomlySampleTests.sh)
