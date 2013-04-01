/* Some code was incorporated from samtools 0.1.18 */
/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
 */
/*
 * BgzfStream.h
 *
 *  Created on: Sep 21, 2012
 *      Author: regan
 */
#ifndef BGZF_STREAM_H
#define BGZF_STREAM_H

#include <string.h>
#include <iostream>
#include <algorithm>
#include <vector>
extern "C" {
#include "bgzf.h"
}
#include "Log.h"
#include <zlib.h>

#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/categories.hpp>


class bgzf_detail  {
public:
	typedef char char_type;
	typedef uint8_t bgzf_byte_t;
	typedef std::vector< int64_t > FileOffsetVector;

	static const int BLOCK_HEADER_LENGTH = 18;
	static const int BLOCK_FOOTER_LENGTH = 8;
	static const int DEFAULT_BLOCK_SIZE = BGZF_MAX_BLOCK_SIZE;

	static const int GZIP_ID1 = 31;
	static const int GZIP_ID2 = 139;
	static const int CM_DEFLATE = 8;
	static const int FLG_FEXTRA = 4;
	static const int OS_UNKNOWN = 255;
	static const int BGZF_ID1 = 66; // 'B'
	static const int BGZF_ID2 = 67; // 'C'
	static const int BGZF_LEN = 2;
	static const int BGZF_XLEN = 6; // BGZF_LEN+4

	static const int GZIP_WINDOW_BITS = -15; // no zlib header
	static const int Z_DEFAULT_MEM_LEVEL = 8;

	class bgzf_compressor_impl {
	public:
		typedef char char_type;

		bgzf_compressor_impl(bool _addEOFBlock = true, bool compress = true) {
			fp = (BGZF*)(calloc(sizeof (BGZF), 1));
	        init(fp, compress);
	        compressed_block_length = 0;
	        compressed_block_written_offset = 0;
	        isEOF = false;
	        addEOFBlock = _addEOFBlock;
	   }

	    ~bgzf_compressor_impl()
	    {
	        reset(fp);
	        blockFileOffsets.clear();
	    }

	    void unsetEofBlock() {
	    	addEOFBlock = false;
	    }

	    void close()
	    {
	        assert(fp->block_offset == 0);
	        assert(compressed_block_length == 0);
	    }

	    // returns true if buffer is empty after operation
	    bool writeBuffer(char *& begin_out, char *end_out)
	    {
	        assert(compressed_block_written_offset <= compressed_block_length);
	        if(compressed_block_length == 0)
	            return true;

	        int bytesLefttoWrite = compressed_block_length - compressed_block_written_offset;
	        int bytesAvailableToWrite = end_out - begin_out;
	        int writeLen = std::min(bytesAvailableToWrite, bytesLefttoWrite);
	        char *buf = (char*)(fp->compressed_block);
	        memcpy(begin_out, buf + compressed_block_written_offset, writeLen);
	        compressed_block_written_offset += writeLen;
	        begin_out += writeLen;
	        if(compressed_block_written_offset == compressed_block_length){
	            compressed_block_written_offset = 0;
	            compressed_block_length = 0;
	        }
	        return compressed_block_length == 0;
	    }

	    // returns true when the input buffer is empty;
	    bool deflateBlock()
	    {
	        assert(compressed_block_length == 0);
	        assert(fp->block_offset <= BGZF_MAX_BLOCK_SIZE);
	        compressed_block_length = deflate_block(fp, fp->block_offset);
	        blockFileOffsets.push_back( fp->block_address );
	        return fp->block_offset == 0;
	    }

	    void readIntoBuffer(const char *end_in, const char *& begin_in)
	    {
	        int bytesInputAvailable = end_in - begin_in;
	        int bytesLeftInBuffer = BGZF_MAX_BLOCK_SIZE - fp->block_offset;
	        int copy_length = std::min(bytesInputAvailable, bytesLeftInBuffer);
	        if(copy_length > 0){
	            char *buffer = (char*)(fp->uncompressed_block);
	            memcpy(buffer + fp->block_offset, begin_in, copy_length);
	            fp->block_offset += copy_length;
	            begin_in += copy_length;
	        }
	    }

	    bool filter(const char *& begin_in, const char *end_in, char *& begin_out, char *end_out, bool flush)
	    {
	        // copy in, if possible
	        bool readIsEmpty = (end_in == begin_in);
	        int block_length = BGZF_MAX_BLOCK_SIZE;
	        readIntoBuffer(end_in, begin_in);
	        // copy out, if possible
			bool writeIsEmpty = writeBuffer(begin_out, end_out);

			if (flush || fp->block_offset == block_length) {
				if(writeIsEmpty) {
					// never write BGZF EOF marker here
					if (fp->block_offset > 0)
						deflateBlock();
					writeIsEmpty &= writeBuffer(begin_out, end_out);
				} else {
					return true;
				}
			}
			if (flush && fp->block_offset == 0 && readIsEmpty && writeIsEmpty && !isEOF) {
				isEOF = true;
				if (addEOFBlock) {
					deflateBlock(); // EOF marker, empty BGZF block
					writeIsEmpty &= writeBuffer(begin_out, end_out);
				}
			}

			if (flush) {
				if (writeIsEmpty & (fp->block_offset == 0)) {
					// all data has been forwarded
					return false;
				} else
					return true;
			} else {
				if (readIsEmpty & writeIsEmpty & (fp->block_offset == 0))
					// all data has been forwarded, no new data read
					return false;
				else
					return true;
			}

		}
	    const FileOffsetVector &getBlockFileOffsets() const {
	    	return blockFileOffsets;
	    }
	private:
		BGZF *fp;
		size_t compressed_block_length, compressed_block_written_offset;
		FileOffsetVector blockFileOffsets;
		bool isEOF, addEOFBlock;
	};

	class bgzf_decompressor_impl {
	public:
		typedef char char_type;
		bgzf_decompressor_impl() {
			fp = (BGZF*) calloc(sizeof(BGZF), 1);
			init(fp);
			uncompressed_block_length = 0;
			uncompressed_block_written_offset = 0;
			isEOF = false;
		}
		~bgzf_decompressor_impl() {
			reset(fp);
		}

		void close() {
			assert(fp->block_offset == 0);
			assert(uncompressed_block_length == 0);
		}
		// returns true if buffer is empty after operation
		bool writeBuffer(char*& begin_out, char *end_out) {
			assert(uncompressed_block_written_offset <= uncompressed_block_length);
			if (uncompressed_block_length == 0)
				return true;
			int bytesLefttoWrite = uncompressed_block_length - uncompressed_block_written_offset;
			int bytesAvailableToWrite = end_out - begin_out;
			int writeLen = std::min(bytesAvailableToWrite, bytesLefttoWrite);
			char *buf = (char*) fp->uncompressed_block;
			memcpy(begin_out, buf + uncompressed_block_written_offset, writeLen);
			uncompressed_block_written_offset += writeLen;
			begin_out += writeLen;
			if (uncompressed_block_written_offset == uncompressed_block_length) {
				uncompressed_block_written_offset = 0;
				uncompressed_block_length = 0;
			}
			return uncompressed_block_length == 0;
		}

		// returns true when the input buffer is empty;
		bool inflateBlock() {
			assert(uncompressed_block_length == 0);
			assert(fp->block_offset >= BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
			assert(fp->block_offset == getBlockLength());
			assert(fp->block_offset <= BGZF_MAX_BLOCK_SIZE);
			uncompressed_block_length = inflate_block(fp, fp->block_offset);
			fp->block_offset = 0;
			return true;
		}

		bool readHeader(const char*& begin_in, const char *end_in) {
			if (fp->block_offset >= BLOCK_HEADER_LENGTH)
				return true;
			int remainingHeader = BLOCK_HEADER_LENGTH - fp->block_offset;
			int availableToRead = end_in - begin_in;
			int readLength = std::min(remainingHeader, availableToRead);
			if (readLength > 0) {
				char *buffer = ((char*) (fp->compressed_block)) + fp->block_offset;
				memcpy(buffer, begin_in, readLength);
				fp->block_offset += readLength;
				begin_in += readLength;
			}
			return fp->block_offset >= BLOCK_HEADER_LENGTH;
		}
		int getBlockLength() {
			assert(fp->block_offset >= BLOCK_HEADER_LENGTH);
			return unpackInt16((uint8_t*) fp->compressed_block + 16) + 1;
		}

		// returns true when block is ready to inflate
		bool readIntoBuffer(const char*& begin_in, const char *end_in) {
			if (!readHeader(begin_in, end_in))
				return false;
			int block_length = getBlockLength();
			int bytesInputAvailable = end_in - begin_in;
			int bytesLeftInBuffer = block_length - fp->block_offset;
			assert(block_length <= BGZF_MAX_BLOCK_SIZE);

			int copy_length = std::min(bytesInputAvailable, bytesLeftInBuffer);
			if (copy_length > 0) {
				char* buffer = (char*) fp->compressed_block;
				memcpy(buffer + fp->block_offset, begin_in, copy_length);
				fp->block_offset += copy_length;
				begin_in += copy_length;
			}
			return fp->block_offset == block_length;
		}
		bool filter(const char*& begin_in, const char *end_in,
				char*& begin_out, char *end_out, bool flush) {
			// copy in, if possible
			bool readIsEmpty = (end_in == begin_in);

			// copy out, if possible
			bool cont = true;
			bool writeIsEmpty, haveFullBlock;
			while (cont) {
				writeIsEmpty = writeBuffer(begin_out, end_out);
				haveFullBlock = readIntoBuffer(begin_in, end_in);
				if (writeIsEmpty & haveFullBlock) {
					inflateBlock();
				} else
					cont = false;
			}

			if (flush) {
				if (writeIsEmpty & (fp->block_offset == 0)) {
					// all data has been forwarded
					return false;
				} else
					return true;
			} else {
				if (readIsEmpty & writeIsEmpty & (fp->block_offset == 0))
					// all data has been forwarded, no new data read
					return false;
				else
					return true;
			}

		}
	private:
		BGZF *fp;
		size_t uncompressed_block_length, uncompressed_block_written_offset;
		bool isEOF;
	};

protected:
	static void init(BGZF *fp, bool compress = true) {
		fp->fp = NULL;
		fp->compress_level = compress ? 7 : 0;
		fp->uncompressed_block = calloc(BGZF_MAX_BLOCK_SIZE, 1);
		fp->compressed_block = calloc(BGZF_MAX_BLOCK_SIZE, 1);
		fp->block_address = 0;
		fp->block_offset = 0;  // uncompressed offset
		fp->block_length = 0;
	}
	static void reset(BGZF *fp) {
		if (fp->uncompressed_block != NULL) free(fp->uncompressed_block);
		if (fp->compressed_block != NULL) free(fp->compressed_block);
		fp->uncompressed_block = NULL;
		fp->compressed_block = NULL;
		if (fp != NULL) free(fp);
		fp = NULL;
	}

public:
	// copied from bgzf.c
	/* The MIT License

	   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

	   Permission is hereby granted, free of charge, to any person obtaining a copy
	   of this software and associated documentation files (the "Software"), to deal
	   in the Software without restriction, including without limitation the rights
	   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	   copies of the Software, and to permit persons to whom the Software is
	   furnished to do so, subject to the following conditions:

	   The above copyright notice and this permission notice shall be included in
	   all copies or substantial portions of the Software.

	   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	   THE SOFTWARE.
	 */
	static const uint8_t *getGmagic() {
		static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";
		return g_magic;
	}

	static inline void packInt16(uint8_t *buffer, uint16_t value)
	{
		buffer[0] = value;
		buffer[1] = value >> 8;
	}

	static inline int unpackInt16(const uint8_t *buffer)
	{
		return buffer[0] | buffer[1] << 8;
	}

	static inline void packInt32(uint8_t *buffer, uint32_t value)
	{
		buffer[0] = value;
		buffer[1] = value >> 8;
		buffer[2] = value >> 16;
		buffer[3] = value >> 24;
	}

	static inline uint32_t unpackInt32(uint8_t * buffer) {
		return (buffer[0] | (buffer[1] << 8) | (buffer[2] << 16) | (buffer[3] << 24));
	}


	// Deflate the block in fp->uncompressed_block into fp->compressed_block. Also adds an extra field that stores the compressed block length.
	static int deflate_block(BGZF *fp, int block_length)
	{
		int comp_size = BGZF_MAX_BLOCK_SIZE;
		if (bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level) != 0) {
			fp->errcode |= BGZF_ERR_ZLIB;
			return -1;
		}
		fp->block_offset = 0;
		return comp_size;
	}

	// Inflate the block in fp->compressed_block into fp->uncompressed_block
	static int inflate_block(BGZF* fp, int block_length)
	{
		z_stream zs;
		zs.zalloc = NULL;
		zs.zfree = NULL;
		zs.next_in = ((Bytef*)fp->compressed_block) + 18;
		zs.avail_in = block_length - 16;
		zs.next_out = (Bytef*) fp->uncompressed_block;
		zs.avail_out = BGZF_MAX_BLOCK_SIZE;

		if (inflateInit2(&zs, -15) != Z_OK) {
			fp->errcode |= BGZF_ERR_ZLIB;
			return -1;
		}
		if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
			inflateEnd(&zs);
			fp->errcode |= BGZF_ERR_ZLIB;
			return -1;
		}
		if (inflateEnd(&zs) != Z_OK) {
			fp->errcode |= BGZF_ERR_ZLIB;
			return -1;
		}
		return zs.total_out;
	}

	static int bgzf_compress(void *_dst, int *dlen, void *src, int slen, int level)
	{
		uint32_t crc;
		z_stream zs;
		uint8_t *dst = (uint8_t*)_dst;

		// compress the body
		zs.zalloc = NULL; zs.zfree = NULL;
		zs.next_in  = (Bytef*)src;
		zs.avail_in = slen;
		zs.next_out = dst + BLOCK_HEADER_LENGTH;
		zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
		if (deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) return -1; // -15 to disable zlib header/footer
		if (deflate(&zs, Z_FINISH) != Z_STREAM_END) return -1;
		if (deflateEnd(&zs) != Z_OK) return -1;
		*dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
		// write the header
		memcpy(dst, getGmagic(), BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
		packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
		// write the footer
		crc = crc32(crc32(0L, NULL, 0L), (const Bytef*) src, slen);
		packInt32((uint8_t*)&dst[*dlen - 8], crc);
		packInt32((uint8_t*)&dst[*dlen - 4], slen);
		return 0;
	}


	static int check_header(const uint8_t *header)
	{
		return (header[0] == 31 && header[1] == 139 && header[2] == 8 && (header[3] & 4) != 0
				&& unpackInt16((uint8_t*)&header[10]) == 6
				&& header[12] == 'B' && header[13] == 'C'
				&& unpackInt16((uint8_t*)&header[14]) == 2);
	}

	//
	// end copied from bgzf.c
	//

	static bool checkInflate(char *compressedBuffer, char *uncompressedBuffer, int block_length, uint32_t crc, uint32_t uncompressedLength) {

		z_stream zs;
		int status;
		zs.msg = Z_NULL;
		zs.opaque = Z_NULL;
		zs.zalloc = Z_NULL;
		zs.zfree = Z_NULL;
		zs.next_in = (Bytef*) compressedBuffer + 18;
		zs.avail_in = block_length - 16;
		zs.next_out = (Bytef*) uncompressedBuffer;
		zs.avail_out = BGZF_MAX_BLOCK_SIZE;

		status = inflateInit2(&zs, GZIP_WINDOW_BITS);
		if (status != Z_OK) {
			fprintf(stderr, "inflate init failed");
			return -1;
		}
		status = inflate(&zs, Z_FINISH);
		if (status != Z_STREAM_END) {
			inflateEnd(&zs);
			fprintf(stderr, "inflate failed");
			return -1;
		}
		status = inflateEnd(&zs);
		if (status != Z_OK) {
			fprintf(stderr, "inflate failed");
			return -1;
		}
		uint32_t _crc = crc32(0L, Z_NULL, 0);
		_crc = crc32(_crc, (Bytef*) uncompressedBuffer, uncompressedLength);
		LOG_DEBUG_OPTIONAL(4, true, "checkInflate: crc: " << crc << " _crc: " << _crc << " length: " << uncompressedLength << " total_out: " << zs.total_out );
		return crc == _crc && uncompressedLength == zs.total_out;

	}

	// returns -1 if there is no block after offset and before eof
	static int64_t getNextBlockFileOffset(void *bgzf_fp, int64_t offset, uint16_t &compressedBlockLength, char *uncompblock, uint32_t &uncompressedBlockLength) {
		FILE *f = (FILE*) bgzf_fp;
		fseek(f, offset, SEEK_SET);
		size_t readSize = BGZF_MAX_BLOCK_SIZE*3;
		uint8_t *bgzfBuff = (uint8_t*) calloc(readSize, 1);
//		buffer[0] = GZIP_ID1;
//		buffer[1] = GZIP_ID2;
//		buffer[2] = CM_DEFLATE;
//		buffer[3] = FLG_FEXTRA;
//		buffer[4] = 0; // mtime
//		buffer[5] = 0;
//		buffer[6] = 0;
//		buffer[7] = 0;
//		buffer[8] = 0;
//		buffer[9] = OS_UNKNOWN;
//		buffer[10] = BGZF_XLEN;
//		buffer[11] = 0;
//		buffer[12] = BGZF_ID1;
//		buffer[13] = BGZF_ID2;
//		buffer[14] = BGZF_LEN;
//		buffer[15] = 0;
//		buffer[16] = 0; // placeholder for block length
//		buffer[17] = 0;
		size_t bytes = fread(bgzfBuff, 1, readSize, f);

		int64_t newOffset = -1;
		for(int i = 0; i < (int) bytes - 34; i++) {
			if (!bgzf_detail::check_header(bgzfBuff + i)) {
				LOG_DEBUG_OPTIONAL(5, true, "check_header failed for " << i);
				continue;
			}
			compressedBlockLength = unpackInt16((uint8_t*) bgzfBuff + i + 16);
			uint32_t crc = unpackInt32((uint8_t*) bgzfBuff + i + compressedBlockLength + 1 - 8);
			uncompressedBlockLength = unpackInt32((uint8_t*) bgzfBuff + i + compressedBlockLength + 1 - 4);
			if (!checkInflate((char*) bgzfBuff + i, uncompblock, compressedBlockLength, crc, uncompressedBlockLength)) {
				LOG_DEBUG_OPTIONAL(2, true, "checkInflate failed for " << i);
				continue;
			}
			newOffset = i + offset;
			if (i + compressedBlockLength < (int) bytes - 34 && ! bgzf_detail::check_header(bgzfBuff + i + compressedBlockLength + 1)) {
				LOG_DEBUG_OPTIONAL(2, true, "check_header for next block failed for " << i << " + " << compressedBlockLength);
				continue;
			}
			break;
		}
		free(bgzfBuff);

		if (newOffset >= 0) {
			LOG_DEBUG_OPTIONAL(2, true, "getNextBlockFileOffset(): Found " << newOffset << " (" << (newOffset - offset) << " tries)");
		} else {
			LOG_DEBUG_OPTIONAL(2, true, "getNextBlockFileOffset(): Did not find a new block after: " << offset);
		}

		return newOffset;
	};

};

template< typename Alloc = std::allocator<char> >
class basic_bgzf_compressor : public boost::iostreams::symmetric_filter< bgzf_detail::bgzf_compressor_impl, Alloc > {
private:
	typedef bgzf_detail::bgzf_compressor_impl                impl_type;
	typedef boost::iostreams::symmetric_filter<impl_type, Alloc>  base_type;
public:
	typedef typename base_type::char_type               char_type;
	typedef typename base_type::category                category;
	basic_bgzf_compressor(bool setEOFBlock = true) : base_type(BGZF_MAX_BLOCK_SIZE, setEOFBlock) {}
};
typedef basic_bgzf_compressor<> bgzf_compressor;

template< typename Alloc = std::allocator<char> >
class basic_bgzf_decompressor : public boost::iostreams::symmetric_filter< bgzf_detail::bgzf_decompressor_impl, Alloc > {
private:
	typedef bgzf_detail::bgzf_decompressor_impl                impl_type;
	typedef boost::iostreams::symmetric_filter<impl_type, Alloc>  base_type;
public:
	typedef typename base_type::char_type               char_type;
	typedef typename base_type::category                category;
	basic_bgzf_decompressor() : base_type(BGZF_MAX_BLOCK_SIZE) {}
};
typedef basic_bgzf_decompressor<> bgzf_decompressor;

class bgzf_ostream : public boost::iostreams::filtering_ostream
{
public:
	template< typename OUT >
	bgzf_ostream(OUT &_os, bool addEOFBlock = true) : comp(addEOFBlock) {
		this->push(comp);
		this->push(_os);
        }
	virtual ~bgzf_ostream() {}
private:
	bgzf_compressor comp;
};

class bgzf_istream : public boost::iostreams::filtering_istream
{
public:
	template< typename IN >
	bgzf_istream(IN &_is) : decomp() {
		this->push(decomp);
		this->push(_is);
        }
	virtual ~bgzf_istream() {}
private:
	bgzf_decompressor decomp;
};


#endif
