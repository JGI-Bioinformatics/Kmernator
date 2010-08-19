//
// Kmernator/src/MmapTempFile.h
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

class MmapTempFile {
	static int getUnique() { static int id = 0; return id++; }
public:
	typedef KoMer::MmapFile MmapFile;
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

	static FileHandle buildNew(size_type size, std::string permamentFile) {
		std::string filename = Options::getTmpDir() + "/Kmmap-";
		filename += boost::lexical_cast<std::string>( getpid() );
		filename += "-" + boost::lexical_cast<std::string>( getUnique() );
		filename += getenv("HOST") == NULL ? "unknown" : getenv("HOST");
		std::cerr << "Creating new tmp file: " << filename << " " << size << std::endl;
		FileHandle fh(filename);
		fh.getOS().seekp(size-1);
		fh.getOS() << '\0';
		if ( !permamentFile.empty() ) {
			permamentFile = Options::getTmpDir() + "/" + permamentFile;
			link(filename.c_str(), permamentFile.c_str());
		}
		return fh;
	}
	static MmapFile buildNewMmap(size_type size, std::string permamentFile = "") {
		FileHandle fh = buildNew(size, permamentFile);
		KoMer::MmapFile mmap(fh.filename, std::ios_base::in | std::ios_base::out, size);
		unlink(fh.filename.c_str());
		return mmap;
	}

	template<MmapTempFile::size_type blockSize, typename T = char>
	class MmapAllocator {
	private:

		typedef T Type;
		typedef KoMer::MmapFile MmapFile;
		typedef KoMer::MmapFileVector MmapFileVector;
		typedef char * CharPtr;
		typedef MmapTempFile::size_type size_type;
		typedef typename boost::simple_segregated_storage< size_type > SSS;

		static const size_type scale = 1024;
		static const size_type overage = 32;

		static MmapAllocator *getSingletons() {
			static MmapAllocator *singletons = NULL;
			if (singletons == NULL) {
				// dangling pointer!!!
				singletons = new MmapAllocator[OMP_MAX_THREADS];
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
				for(size_type thread = 0; thread < OMP_MAX_THREADS; thread++) {
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

//
// $Log: MmapTempFile.h,v $
// Revision 1.2  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.1.2.1  2010-05-12 18:25:04  regan
// refactored
//
