//
// Kmernator/src/config.h
//
// Author: Rob Egan, Craig Furman
//
// Copyright 2010 The Regents of the University of California.
// All rights reserved.
//
// The United States Government has rights in this work pursuant
// to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
// W-7405-ENG-48 between the United States Department of Energy
// and the University of California.
//
// Redistribution and use in source and binary forms are permitted
// provided that: (1) source distributions retain this entire
// copyright notice and comment, and (2) distributions including
// binaries display the following acknowledgement:  "This product
// includes software developed by the University of California,
// JGI-PSF and its contributors" in the documentation or other
// materials provided with the distribution and in all advertising
// materials mentioning features or use of this software.  Neither the
// name of the University nor the names of its contributors may be
// used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE.
//


#ifndef _KMERNATOR_CONFIG_H
#define _KMERNATOR_CONFIG_H

#include "../version.h"

#include <vector>

#include <boost/cstdint.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#ifdef ENABLE_MPI

#define _USE_MPI
#include <boost/serialization/string.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include <mpi.h>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

#ifdef ENABLE_OPENMP

#define _USE_OPENMP
//#define _USE_THREADSAFE_KMER
#include <omp.h>

const int MAX_FILE_PARALLELISM = 4;

#else

inline int omp_get_max_threads() { return 1; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_thread_num()  { return 0; }
inline void omp_set_nested(int i) {}
inline int omp_get_nested() { return 0; }
inline void omp_set_dynamic(int i) {}
inline int omp_get_dynamic() { return 0; }
inline bool omp_in_parallel() { return false; }
inline void omp_set_num_threads(int t) { assert(t==1); }

const int MAX_FILE_PARALLELISM = 1;

#endif



const int OMP_NESTED_DEFAULT = omp_get_nested();
const int OMP_DYNAMIC_DEFAULT = omp_get_dynamic();
const int OMP_MAX_THREADS_DEFAULT = omp_get_max_threads();

// Some processes need memory to make a two bit sequence.
// if more than 128 bytes (512 sequence length) is needed, malloc will be called.
#define MAX_STACK_SIZE 1024
#define STACK_ALLOC(_TYPE, _VAR, _LENGTH) \
	bool _VAR_needMalloc = (_LENGTH * sizeof(_TYPE)) > MAX_STACK_SIZE; \
	_TYPE _VAR_buffer[_VAR_needMalloc ? 0 : _LENGTH]; \
	_TYPE *_VAR = _VAR_buffer; \
	if (_VAR_needMalloc) _VAR = new _TYPE[_LENGTH];
#define STACK_DEALLOC(_VAR) \
		if (_VAR_needMalloc) delete [] _VAR;

#ifndef MAX_KMER_MAP_BUCKETS
#define MAX_KMER_MAP_BUCKETS 33554432
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <cassert>

namespace Kmernator {
typedef const char * RecordPtr;
static const char REF_QUAL = 126;
static const char FASTQ_START_CHAR_ILLUMINA = 64;
static const char FASTQ_START_CHAR_STD = 33;
static const char PRINT_REF_QUAL = REF_QUAL;

typedef boost::uint8_t  UI8;
typedef boost::uint16_t UI16;
typedef boost::uint32_t UI32;
typedef boost::uint64_t UI64;

typedef boost::int8_t  I8;
typedef boost::int16_t I16;
typedef boost::int32_t I32;
typedef boost::int64_t I64;

#define MAX_I8   127
#define MAX_UI8  255u

#define MAX_I16  32767
#define MAX_UI16 65535u

#define MAX_I32  2147483647
#define MAX_UI32 4294967295u

#ifdef IS_64_BIT
#define MAX_I64  9223372036854775807
#define MAX_UI64 18446744073709551615u
#else
#error "only 64-bit environments are supported"
#define MAX_I64  MAX_I32
#define MAX_UI64 MAX_UI32
#endif


typedef UI32 SequenceLengthType;
#define MPISequenceLengthType MPI_UNSIGNED_LONG
#define MAX_SEQUENCE_LENGTH MAX_UI32
typedef UI8 SequenceLengthType1;
typedef UI16 SequenceLengthType2;

typedef UI64 ReadSetSizeType;
#define MPIReadSetSizeType MPI_UNSIGNED_LONG_LONG
#define MAX_READ_SET_SIZE   MAX_UI64

typedef UI64  KmerNumberType;
typedef SequenceLengthType  KmerIndexType;
#define MAX_KMER_INDEX MAX_SEQUENCE_LENGTH
typedef I64  KmerSizeType;


typedef boost::iostreams::mapped_file MmapFile;
typedef boost::iostreams::mapped_file_source MmapSource;
typedef boost::iostreams::stream< MmapSource > MmapIStream;
typedef boost::iostreams::filtering_istreambuf FilteredIStream;
typedef std::vector< MmapFile > MmapFileVector;
typedef std::vector< MmapSource > MmapSourceVector;


};

#endif // _KMERNATOR_CONFIG_H
