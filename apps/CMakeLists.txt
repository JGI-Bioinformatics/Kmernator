# $Header: /repository/PI_annex/robsandbox/KoMer/apps/CMakeLists.txt,v 1.9 2010-08-18 17:47:49 regan Exp $
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
                              boost_program_options-mt-static 
                              boost_iostreams-mt-static
                              boost_thread-mt-static
                              )


set_source_files_properties( FilterReads
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( FilterReads FilterReads )
target_link_libraries( FilterReads 
                              Sequence
                              Kmer
                              ReadSet
                              boost_program_options-mt-static
                              boost_iostreams-mt-static
                              boost_thread-mt-static
                              )
                           
set_source_files_properties( TnfDistance
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( TnfDistance TnfDistance )
target_link_libraries( TnfDistance 
                              Sequence
                              Kmer
                              ReadSet
                              boost_program_options-mt-static
                              boost_iostreams-mt-static
                              boost_thread-mt-static
                              )
   
set_source_files_properties( Fastq2Fasta
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( Fastq2Fasta Fastq2Fasta )
target_link_libraries( Fastq2Fasta
                              Sequence
                              ReadSet
                              boost_program_options-mt-static
                              boost_iostreams-mt-static
                              boost_thread-mt-static
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
                              boost_thread-mt-static
                              )


set_source_files_properties( FixPair
                             PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS}
                             )
add_executable( FixPair FixPair )
target_link_libraries( FixPair
                              Sequence
                              ReadSet
                              boost_program_options-mt-static
                              boost_iostreams-mt-static
                              boost_thread-mt-static
                              )


#
# $Log: CMakeLists.txt,v $
# Revision 1.9  2010-08-18 17:47:49  regan
# added TNF distance application
#
# Revision 1.8.16.1  2010-07-13 17:39:12  regan
# created TnfDistance application
#
# Revision 1.8  2010-05-05 06:28:40  regan
# merged changes from FixPairOutput-20100504
#
# Revision 1.7.4.1  2010-05-04 21:33:44  regan
# checkpoint
#
# Revision 1.7  2010-05-01 21:58:00  regan
# merged head with serial threaded build partitioning
#
# Revision 1.6.4.1  2010-04-30 23:33:42  regan
# replaced read cache with LRU 3rd party cache
#
# Revision 1.6  2010-04-16 22:44:23  regan
# merged HEAD with changes for mmap and intrusive pointer
#
# Revision 1.4.2.2  2010-04-05 05:42:56  regan
# checkpoint mmaping input files
#
# Revision 1.4.2.1  2010-04-05 03:32:13  regan
# moved read file reader
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
 
