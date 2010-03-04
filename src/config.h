#ifndef _KOMER_CONFIG_H
#define _KOMER_CONFIG_H

#ifdef ENABLE_OPENMP
#define _USE_OPENMP
//#define _USE_THREADSAFE_KMER
#include <omp.h>
const int OMP_NESTED_DEFAULT = omp_get_nested();
const int OMP_DYNAMIC_DEFAULT = omp_get_dynamic();
const int MAX_FILE_PARALLELISM = 4;
#endif

#endif

