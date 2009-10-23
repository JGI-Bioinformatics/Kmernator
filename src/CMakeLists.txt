# $Header: /repository/PI_annex/robsandbox/KoMer/src/CMakeLists.txt,v 1.5 2009-10-23 20:32:50 cfurman Exp $
#
# KoMer/src/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)

add_library( TwoBitSequence TwoBitSequence )
add_library( Sequence Sequence TwoBitSequence )
add_library( Kmer Kmer TwoBitSequence)
add_library( ReadSet ReadSet )

#
# $Log: CMakeLists.txt,v $
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
