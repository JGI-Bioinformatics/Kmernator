#ifndef _KOMER_CONFIG_H
#define _KOMER_CONFIG_H

#include <vector>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

#ifdef ENABLE_OPENMP
#define _USE_OPENMP
//#define _USE_THREADSAFE_KMER
#include <omp.h>
const int OMP_NESTED_DEFAULT = omp_get_nested();
const int OMP_DYNAMIC_DEFAULT = omp_get_dynamic();
const int MAX_FILE_PARALLELISM = 4;
const int OMP_MAX_THREADS = omp_get_max_threads();
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_thread_num()  { return 0; }
const int MAX_FILE_PARALLELISM = 1;
const int OMP_MAX_THREADS = 1;
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <cassert>

namespace KoMer {
   typedef const char * RecordPtr;
   static const char REF_QUAL = 0xff;
   static const char FASTQ_START_CHAR = 64;
   static const char PRINT_REF_QUAL = 126;

   typedef boost::iostreams::mapped_file MmapFile;
   typedef boost::iostreams::mapped_file_source MmapSource;
   typedef boost::iostreams::stream< MmapSource > MmapIStream;
   typedef std::vector< MmapFile > MmapFileVector;


};

#endif

