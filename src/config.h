#ifndef _KOMER_CONFIG_H
#define _KOMER_CONFIG_H

#if(1)
#define _USE_OPENMP
//#define _USE_THREADSAFE_KMER
#include <omp.h>
const int OMP_NESTED_DEFAULT = omp_get_nested();
const int OMP_DYNAMIC_DEFAULT = omp_get_dynamic();
#endif

#endif

