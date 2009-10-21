# $Header: /repository/PI_annex/robsandbox/KoMer/test/CMakeLists.txt,v 1.4 2009-10-21 06:51:37 regan Exp $
#
# KoMer/test/CmakeLists.txt
#

project( KoMer )
cmake_minimum_required(VERSION 2.6)

add_executable( ktest2 ktest2 )
target_link_libraries(ktest2 Sequence ReadSet)

#
# $Log: CMakeLists.txt,v $
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
