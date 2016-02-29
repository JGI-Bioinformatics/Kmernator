//
// Kmernator/src/config.h
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/


#ifndef _KMERNATOR_CONFIG_H
#define _KMERNATOR_CONFIG_H

#include "version.h"

#include <vector>

#include <boost/cstdint.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#ifdef ENABLE_MPI

#define _USE_MPI
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
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
inline int omp_get_num_procs() { return 1; }
inline int omp_get_thread_limit() { return 1024; }
inline int omp_get_thread_num()  { return 0; }
inline void omp_set_nested(int i) {}
inline int omp_get_nested() { return 0; }
inline void omp_set_dynamic(int i) {}
inline int omp_get_dynamic() { return 0; }
inline bool omp_in_parallel() { return false; }
inline void omp_set_num_threads(int t) { assert(t==1); }
inline int omp_get_level() { return 1; }

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
#define MAX_KMER_MAP_BUCKETS 67108864
#endif

#ifndef DEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

namespace Kmernator {
typedef const char * RecordPtr;
static const boost::uint8_t FASTQ_START_CHAR_ILLUMINA = 64;
static const boost::uint8_t FASTQ_START_CHAR_STD = 33;
static const boost::uint8_t FASTQ_START_CHAR_DEFAULT = FASTQ_START_CHAR_STD;
static const boost::uint8_t PRINT_REF_QUAL = FASTQ_START_CHAR_STD + 70; // 64 + 39 or 33 + 70 == 103 == 'j' - high by any base level
static const boost::uint8_t REF_QUAL = 127;

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

template<typename T>
class ScopedTempValue {
public:
	ScopedTempValue(T &variable, T value) : _var(variable), _oldValue(variable) {
		_var = value;
	}
	~ScopedTempValue() {
		_var = _oldValue;
	}
	operator T () {
		return _oldValue;
	}
private:
	T &_var;
	T _oldValue;
};


#endif // _KMERNATOR_CONFIG_H
