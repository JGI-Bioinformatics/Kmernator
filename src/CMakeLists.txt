# KoMer CmakeLists.txt
project( KoMer )
cmake_minimum_required(VERSION 2.6)

add_library( Sequence Sequence )
add_library( ReadSet ReadSet )

add_executable( ktest2 ktest2 )
