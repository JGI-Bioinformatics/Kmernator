//
// Kmernator/src/MmapTempFile.h
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

#ifndef _MMAP_TEMP_FILE_H
#define _MMAP_TEMP_FILE_H

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>

#include <boost/lexical_cast.hpp>
#include <boost/pool/simple_segregated_storage.hpp>
#include <boost/shared_ptr.hpp>

#include "config.h"
#include "Options.h"
#include "Log.h"
#include "Utils.h"

class MmapTempFile {
public:
	typedef Kmernator::MmapFile MmapFile;
	typedef unsigned long size_type;

	class FileHandle {
	public:
		typedef boost::shared_ptr< std::ofstream > OSPtr;
		std::string filename;
		OSPtr osPtr;
		MmapFile mmap;

		FileHandle(std::string _filename)
		: filename(_filename) {
			osPtr.reset(new std::ofstream(filename.c_str()));
		}
		FileHandle(const FileHandle &copy)
		: filename(copy.filename), osPtr(copy.osPtr), mmap(copy.mmap) {}
		FileHandle &operator=(const FileHandle &copy) {
			if (this == &copy)
				return *this;
			filename = copy.filename;
			osPtr = copy.osPtr;
			mmap = copy.mmap;
			return *this;
		}
		std::ofstream &getOS() {
			return *osPtr;
		}
	};

	static FileHandle buildNew(size_type size, std::string permanentFile) {
		std::string filename;
		if (permanentFile.empty())
			filename = Options::getOptions().getTmpDir() + UniqueName::generateUniqueName("/.tmp-Kmmap-");
		else
			filename = permanentFile;
		LOG_DEBUG_OPTIONAL(1, true, "Creating new file: " << filename << " " << size);
		FileHandle fh(filename);
		fh.getOS().seekp(size-1);
		fh.getOS() << '\0';
		return fh;
	}
	static MmapFile buildNewMmap(size_type size, std::string permanentFile = "") {
		FileHandle fh = buildNew(size, permanentFile);
		Kmernator::MmapFile mmap(fh.filename, std::ios_base::in | std::ios_base::out, size);
		LOG_DEBUG_OPTIONAL(1, true, "Created mmap with alignment " << mmap.alignment() << " at " << fh.filename);
		if (permanentFile.empty())
			unlink(fh.filename.c_str());
		return mmap;
	}

	static MmapFile openMmap(std::string filename) {
		MmapFile mmap;
		if (FileUtils::fileExists(filename)) {
			mmap = MmapFile(filename, std::ios_base::in | std::ios_base::out);
			assert(mmap.is_open());
			assert(mmap.data() != NULL);
			assert(mmap.size() > 0);
		}
		return mmap;
	}

	template<MmapTempFile::size_type blockSize, typename T = char>
	class MmapAllocator {
	private:

		typedef T Type;
		typedef Kmernator::MmapFile MmapFile;
		typedef Kmernator::MmapFileVector MmapFileVector;
		typedef char * CharPtr;
		typedef MmapTempFile::size_type size_type;
		typedef typename boost::simple_segregated_storage< size_type > SSS;

		static const size_type scale = 1024;
		static const size_type overage = 32;

		static MmapAllocator *getSingletons() {
			static MmapAllocator *singletons = NULL;
			if (singletons == NULL) {
				// dangling pointer!!!
				singletons = new MmapAllocator[omp_get_max_threads()];
			}
			return singletons;
		}
		static MmapAllocator &getSingleton() {
			return *(getSingletons() + omp_get_thread_num());
		}

		SSS sss;
		MmapFileVector mmaps;

		inline size_type _getBlocks(const size_type bytes) {
			return (bytes + (blockSize-1)) / blockSize;
		}

		void _addMap(const size_type bytes) {
			size_type blocks = scale;
			if (bytes * overage > blocks * blockSize)
				blocks = _getBlocks( bytes * overage );
			if (!mmaps.empty()) {
				if (blocks * blockSize < mmaps.back().size() * 1.3) {
					blocks = _getBlocks( mmaps.back().size() * 1.3);
				}
			}

			size_type headerSize = (blocks+1) * sizeof(size_type);
			MmapFile newMap = MmapTempFile::buildNewMmap( headerSize + (blocks * blockSize));

			// store headerSize in first size holder
			size_type &_headerSize = *_getVal( newMap.data() );
			_headerSize = headerSize;

			mmaps.push_back(newMap);
			sss.add_block(newMap.data() + headerSize, newMap.size() - headerSize, blockSize);
		}

		static inline size_type *_getVal(CharPtr ptr) {
			return (size_type*) ptr;
		}
		size_type *_getBlockAllocation(CharPtr allocation) {
			for(size_type mapIdx = 0; mapIdx < mmaps.size(); mapIdx++) {
				size_type *blockAllocation = _getBlockAllocation(mmaps[mapIdx], allocation);
				if (blockAllocation != NULL)
					return blockAllocation;
			}
			return NULL;
		}
		static size_type *_getBlockAllocation(MmapFile &mmap, CharPtr allocation) {
			if ( allocation >= mmap.data() && allocation < (mmap.data()+mmap.size()) ) {
				size_type *headerSize = _getVal(mmap.data());
				size_type blocks = (allocation - (mmap.data() + *headerSize)) / blockSize;
				return headerSize + blocks + 1;
			} else {
				return NULL;
			}
		}

		void _setSize(CharPtr allocation, const size_type blocks) {
			size_type *blockAllocation = _getBlockAllocation(allocation);
			*blockAllocation = blocks;
		}

		CharPtr _malloc(const size_type bytes) {

			size_type blocks = _getBlocks( bytes );

			CharPtr allocation = NULL;
			int iterations = 0;
			while (iterations++ < 5 && allocation == NULL) {
				if (iterations > 1)
					_addMap(bytes);
				allocation = (CharPtr) sss.malloc_n(blocks, blockSize);
			}
			_setSize(allocation, blocks);
			return allocation;
		}

		void _free(CharPtr allocation) {
			size_type *blockAllocation = _getBlockAllocation(allocation);
			if (blockAllocation == NULL) {
				// iterate through all thread allocations
				MmapAllocator *singletons = getSingletons();
				for(size_type thread = 0; thread < omp_get_max_threads(); thread++) {
					if (thread == omp_get_thread_num())
						continue;
					blockAllocation = (singletons + thread)->_getBlockAllocation(allocation);
					if (blockAllocation != NULL) {
						size_type blocks = *blockAllocation;
						(singletons + thread)->sss.free_n(allocation, blocks, blockSize);
						return;
					}
				}
				throw("Should never get here");
			} else {
				size_type blocks = *blockAllocation;
				sss.free_n(allocation, blocks, blockSize);
			}
		}
		MmapAllocator() {
		}

	public:
		~MmapAllocator() {
		}

		static char * malloc(const size_type bytes) {
			return getSingleton()._malloc(bytes);
		}
		static void free(char * const block) {
			return getSingleton()._free(block);
		}
	};

};

#endif

