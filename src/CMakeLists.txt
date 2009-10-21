# $Header: /repository/PI_annex/robsandbox/KoMer/src/CMakeLists.txt,v 1.3 2009-10-21 06:51:34 regan Exp $
#
# KoMer/src/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)

add_library( Sequence Sequence TwoBitSequence )
add_library( ReadSet ReadSet )

#
# $Log: CMakeLists.txt,v $
# Revision 1.3  2009-10-21 06:51:34  regan
# bug fixes
# build lookup tables for twobitsequence
#
# Revision 1.2  2009-10-20 17:25:50  regan
# added CVS tags
#
#
