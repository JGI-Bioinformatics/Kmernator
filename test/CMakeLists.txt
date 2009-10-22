# $Header: /repository/PI_annex/robsandbox/KoMer/test/CMakeLists.txt,v 1.5 2009-10-22 07:04:03 regan Exp $
#
# KoMer/test/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)
enable_testing()

add_executable( TwoBitSequenceTest TwoBitSequenceTest )
target_link_libraries( TwoBitSequenceTest TwoBitSequence 
                                          boost_test_exec_monitor-mt 
                                          boost_unit_test_framework-mt )
add_test( TwoBitSequenceTest TwoBitSequenceTest )

add_executable( ktest2 ktest2 )
target_link_libraries( ktest2 Sequence 
                              ReadSet )
add_test(ktest2 ktest2 ${CMAKE_SOURCE_DIR}/test/10.fastq)


#
# $Log: CMakeLists.txt,v $
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
