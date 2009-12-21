# $Header: /repository/PI_annex/robsandbox/KoMer/src/CMakeLists.txt,v 1.10 2009-12-21 07:54:18 regan Exp $
#
# KoMer/src/CmakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)

IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

set_source_files_properties( ReadSet
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )

set_source_files_properties( Kmer
                            PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                          )
                          
add_library( TwoBitSequence TwoBitSequence )
add_library( Sequence Sequence TwoBitSequence )
add_library( Kmer Kmer TwoBitSequence)
add_library( ReadSet ReadSet )

#
# $Log: CMakeLists.txt,v $
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
