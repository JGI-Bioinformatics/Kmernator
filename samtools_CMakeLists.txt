#
# Kmernator/samtools/CMakeLists.txt
#
# Copyright 2012 The Regents of the University of California.
# All rights reserved.
#
# Author: Rob Egan
#
# The United States Government has rights in this work pursuant
# to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
# W-7405-ENG-48 between the United States Department of Energy
# and the University of California.
#
# Redistribution and use in source and binary forms are permitted
# provided that: (1) source distributions retain this entire
# copyright notice and comment, and (2) distributions including
# binaries display the following acknowledgement:  "This product
# includes software developed by the University of California,
# JGI-PSF and its contributors" in the documentation or other
# materials provided with the distribution and in all advertising
# materials mentioning features or use of this software.  Neither the
# name of the University nor the names of its contributors may be
# used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE.
#

cmake_minimum_required(VERSION 2.6)

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)

IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

SET(CMAKE_CXX_FLAGS "-fno-strict-aliasing ${BIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

add_library(bam STATIC bam)
add_library(bgzf STATIC bgzf)
add_library(kstring STATIC kstring)
add_library(bam_aux STATIC bam_aux)
add_library(bam_import STATIC bam_import)
add_library(sam STATIC sam)
add_library(bam_index STATIC bam_index)
add_library(razf STATIC razf)
add_library(faidx STATIC faidx)
add_library(bedidx STATIC bedidx)
add_library(knetfile STATIC knetfile)
add_library(bam_sort STATIC bam_sort)
add_library(sam_header STATIC sam_header)
add_library(bam_reheader STATIC bam_reheader)
add_library(kprobaln STATIC kprobaln)
add_library(bam_cat STATIC bam_cat)
add_library(bam_md STATIC bam_md)
add_library(bam_pileup STATIC bam_pileup)

add_library(libbam 
    STATIC 
        bam sam bgzf bam_sort sam_header bam_aux kstring 
        faidx bam_import bam_pileup razf bam_index) 
# kstring   bam_aux bam_import bam_index razf faidx bedidx 
#knetfile bam_sort sam_header bam_reheader kprobaln bam_cat bam_md 
#)
