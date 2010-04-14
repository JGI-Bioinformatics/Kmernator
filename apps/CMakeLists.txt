# $Header: /repository/PI_annex/robsandbox/KoMer/apps/CMakeLists.txt,v 1.5 2010-04-14 15:58:34 regan Exp $
#
 
cmake_minimum_required(VERSION 2.6)
 

# compiler flags
IF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
  SET(CXX_FLAGS $ENV{CXX_FLAGS})
ENDIF(ENV{CXX_FLAGS} AND NOT CXX_FLAGS)
IF(CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "Default compiler flags" FORCE)
ENDIF(CXX_FLAGS)

set_source_files_properties( CompareSpectrums
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( CompareSpectrums CompareSpectrums )
target_link_libraries( CompareSpectrums 
                              Sequence 
                              Kmer
                              ReadSet
                              boost_program_options-mt-static )


set_source_files_properties( FilterReads
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( FilterReads FilterReads )
target_link_libraries( FilterReads 
                              Sequence
                              Kmer
                              ReadSet
                              boost_program_options-mt-static )
                              
add_executable( Fastq2Fasta Fastq2Fasta )
target_link_libraries( Fastq2Fasta
                              Sequence
                              ReadSet
                              boost_program_options-mt-static
                              )
set_source_files_properties( Fastq2FastaQual
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( Fastq2FastaQual Fastq2FastaQual )
target_link_libraries( Fastq2FastaQual
                              Sequence
                              ReadSet
                              boost_program_options-mt-static
                              boost_iostreams-mt-static
                              )


#
# $Log: CMakeLists.txt,v $
# Revision 1.5  2010-04-14 15:58:34  regan
# fixed and added fastq2fastaqual
#
# Revision 1.4  2010-03-16 18:58:10  regan
# created a simple fastq2fasta program
#
# Revision 1.3  2010-01-13 23:32:35  regan
# new app
#
# Revision 1.2  2010-01-08 18:34:51  regan
# refactored a bit
# enabled openmp parallelizations
#
# Revision 1.1  2009-11-27 23:28:07  cfurman
# CompareSpectrum application added
#
 
