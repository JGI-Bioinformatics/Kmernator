# $Header: /repository/PI_annex/robsandbox/KoMer/test/CMakeLists.txt,v 1.17 2010-02-22 14:41:18 regan Exp $
#
# KoMer/test/CmakeLists.txt
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


add_executable( TwoBitSequenceTest TwoBitSequenceTest )
target_link_libraries( TwoBitSequenceTest TwoBitSequence
                                          boost_test_exec_monitor-mt 
                                          boost_unit_test_framework-mt )
add_test( TwoBitSequenceTest TwoBitSequenceTest )


add_executable( ReadSetTest ReadSetTest )
target_link_libraries( ReadSetTest Sequence
                                   ReadSet 
                                   boost_test_exec_monitor-mt 
                                   boost_unit_test_framework-mt )
add_test( ReadSetTest ReadSetTest )
                            
set_source_files_properties( KmerTest
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )
                                                                                                           
add_executable( KmerTest KmerTest )
target_link_libraries( KmerTest TwoBitSequence Kmer Sequence ReadSet
                                boost_test_exec_monitor-mt 
                                boost_unit_test_framework-mt
                                 )
add_test( KmerTest KmerTest )

set_source_files_properties( ktest2
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )

add_executable( ktest2 ktest2 )
target_link_libraries( ktest2 Sequence Kmer
                              ReadSet
                              boost_program_options-mt-static )
add_test(ktest2 ktest2 21 10.fastq)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fastq test/10.fastq )

# v2.8: file(COPY test/500.fasta test/500.qual DESTINATION test)
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.fasta test/ )
execute_process(COMMAND cp -p ${CMAKE_SOURCE_DIR}/test/10.qual test/ )

#
# $Log: CMakeLists.txt,v $
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
