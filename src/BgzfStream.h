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
#include "bgzf.h"
#include <zlib.h>

#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/filtering_stream.hpp>



class bgzf_detail  {
public:
	typedef char char_type;
	typedef int8_t bgzf_byte_t;
	static const int DEFAULT_BLOCK_SIZE = 64 * 1024;
	static const int MAX_BLOCK_SIZE = 64 * 1024;
	static const int BLOCK_HEADER_LENGTH = 18;
	static const int BLOCK_FOOTER_LENGTH = 8;

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
		static bool &getAddEOFBlock() {
			static bool _x = true;
			return _x;
		}
		static void setAddEOFBlock(bool v) {
			getAddEOFBlock() = v;
		}
		bgzf_compressor_impl(bool compress = true) {
			fp = (BGZF*)(calloc(sizeof (BGZF), 1));
	        init(fp, compress);
	        compressed_block_length = 0;
	        compressed_block_written_offset = 0;
	        isEOF = false;
	   }

	    ~bgzf_compressor_impl()
	    {
	        reset(fp);
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
	        int block_length = fp->uncompressed_block_size;
	        assert(fp->block_offset <= block_length);
	        compressed_block_length = deflate_block(fp, fp->block_offset);
	        return fp->block_offset == 0;
	    }

	    void readIntoBuffer(const char *end_in, const char *& begin_in)
	    {
	        int bytesInputAvailable = end_in - begin_in;
	        int bytesLeftInBuffer = fp->uncompressed_block_size - fp->block_offset;
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
	        int block_length = fp->uncompressed_block_size;
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
				if (getAddEOFBlock()) {
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
	private:
		BGZF *fp;
		size_t compressed_block_length, compressed_block_written_offset;
		bool isEOF;
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
			assert(fp->block_offset <= fp->compressed_block_size);
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
			int max_block_length = fp->compressed_block_size;
			int block_length = getBlockLength();
			int bytesInputAvailable = end_in - begin_in;
			int bytesLeftInBuffer = block_length - fp->block_offset;
			assert(block_length <= max_block_length);

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
		fp->file_descriptor = -1;
		fp->open_mode = 'w';
		fp->owned_file = 0;
		fp->compress_level = compress ? 7 : 0;
		fp->uncompressed_block_size = DEFAULT_BLOCK_SIZE;
		fp->uncompressed_block = calloc(MAX_BLOCK_SIZE, 1);
		fp->compressed_block_size = MAX_BLOCK_SIZE;
		fp->compressed_block = calloc(MAX_BLOCK_SIZE, 1);
		fp->block_address = 0;
		fp->block_offset = 0;  // uncompressed offset
		fp->block_length = 0;
		fp->error = NULL;
	}
	static void reset(BGZF *fp) {
		if (fp->uncompressed_block != NULL) free(fp->uncompressed_block);
		if (fp->compressed_block != NULL) free(fp->compressed_block);
		if (fp != NULL) free(fp);
		fp->uncompressed_block = NULL;
		fp->compressed_block = NULL;
		fp = NULL;
	}

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

	static inline
	void
	packInt16(uint8_t* buffer, uint16_t value)
	{
		buffer[0] = value;
		buffer[1] = value >> 8;
	}

	static inline
	int
	unpackInt16(const uint8_t* buffer)
	{
		return (buffer[0] | (buffer[1] << 8));
	}

	static inline
	void
	packInt32(uint8_t* buffer, uint32_t value)
	{
		buffer[0] = value;
		buffer[1] = value >> 8;
		buffer[2] = value >> 16;
		buffer[3] = value >> 24;
	}

	static inline
	int
	bgzf_min(int x, int y)
	{
		return (x < y) ? x : y;
	}

	static
	void
	report_error(BGZF* fp, const char* message) {
		fp->error = message;
	}

	static
	int
	deflate_block(BGZF* fp, int block_length)
	{
		// Deflate the block in fp->uncompressed_block into fp->compressed_block.
		// Also adds an extra field that stores the compressed block length.

		bgzf_byte_t* buffer = (bgzf_byte_t*) fp->compressed_block;
		int buffer_size = fp->compressed_block_size;

		// Init gzip header
		buffer[0] = GZIP_ID1;
		buffer[1] = GZIP_ID2;
		buffer[2] = CM_DEFLATE;
		buffer[3] = FLG_FEXTRA;
		buffer[4] = 0; // mtime
		buffer[5] = 0;
		buffer[6] = 0;
		buffer[7] = 0;
		buffer[8] = 0;
		buffer[9] = OS_UNKNOWN;
		buffer[10] = BGZF_XLEN;
		buffer[11] = 0;
		buffer[12] = BGZF_ID1;
		buffer[13] = BGZF_ID2;
		buffer[14] = BGZF_LEN;
		buffer[15] = 0;
		buffer[16] = 0; // placeholder for block length
		buffer[17] = 0;

		// loop to retry for blocks that do not compress enough
		int input_length = block_length;
		int compressed_length = 0;
		while (1) {
			z_stream zs;
			zs.zalloc = NULL;
			zs.zfree = NULL;
			zs.next_in = (Bytef*) fp->uncompressed_block;
			zs.avail_in = input_length;
			zs.next_out = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
			zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

			int status = deflateInit2(&zs, fp->compress_level, Z_DEFLATED,
					GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
			if (status != Z_OK) {
				report_error(fp, "deflate init failed");
				return -1;
			}
			status = deflate(&zs, Z_FINISH);
			if (status != Z_STREAM_END) {
				deflateEnd(&zs);
				if (status == Z_OK) {
					// Not enough space in buffer.
					// Can happen in the rare case the input doesn't compress enough.
					// Reduce the amount of input until it fits.
					input_length -= 1024;
					if (input_length <= 0) {
						// should never happen
						report_error(fp, "input reduction failed");
						return -1;
					}
					continue;
				}
				report_error(fp, "deflate failed");
				return -1;
			}
			status = deflateEnd(&zs);
			if (status != Z_OK) {
				report_error(fp, "deflate end failed");
				return -1;
			}
			compressed_length = zs.total_out;
			compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
			//std::cout << " total_out: " << zs.total_out << " ";
			if (compressed_length > MAX_BLOCK_SIZE) {
				// should never happen
				report_error(fp, "deflate overflow");
				return -1;
			}
			break;
		}

		packInt16((uint8_t*)&buffer[16], compressed_length-1);
		uint32_t crc = crc32(0L, NULL, 0L);
		crc = crc32(crc, (const Bytef*)fp->uncompressed_block, input_length);
		packInt32((uint8_t*)&buffer[compressed_length-8], crc);
		packInt32((uint8_t*)&buffer[compressed_length-4], input_length);


		int remaining = block_length - input_length;
		if (remaining > 0) {
			if (remaining > input_length) {
				// should never happen (check so we can use memcpy)
				report_error(fp, "remainder too large");
				return -1;
			}
			memcpy(fp->uncompressed_block,
					((Bytef*)fp->uncompressed_block) + input_length,
					remaining);
		}
		fp->block_offset = remaining;
		//std::cout << "compressed length " << compressed_length << " crc: " << crc << " input_length: " << input_length << " remaining: " << remaining << std::endl;
		return compressed_length;
	}

	static
	int
	inflate_block(BGZF* fp, int block_length)
	{
		// Inflate the block in fp->compressed_block into fp->uncompressed_block

		z_stream zs;
		int status;
		zs.zalloc = NULL;
		zs.zfree = NULL;
		zs.next_in = (Bytef*) fp->compressed_block + 18;
		zs.avail_in = block_length - 16;
		zs.next_out = (Bytef*) fp->uncompressed_block;
		zs.avail_out = fp->uncompressed_block_size;

		status = inflateInit2(&zs, GZIP_WINDOW_BITS);
		if (status != Z_OK) {
			report_error(fp, "inflate init failed");
			return -1;
		}
		status = inflate(&zs, Z_FINISH);
		if (status != Z_STREAM_END) {
			inflateEnd(&zs);
			report_error(fp, "inflate failed");
			return -1;
		}
		status = inflateEnd(&zs);
		if (status != Z_OK) {
			report_error(fp, "inflate failed");
			return -1;
		}
		return zs.total_out;
	}

	static
	int
	check_header(const bgzf_byte_t* header)
	{
		return (header[0] == GZIP_ID1 &&
				header[1] == (bgzf_byte_t) GZIP_ID2 &&
				header[2] == Z_DEFLATED &&
				(header[3] & FLG_FEXTRA) != 0 &&
				unpackInt16((uint8_t*)&header[10]) == BGZF_XLEN &&
				header[12] == BGZF_ID1 &&
				header[13] == BGZF_ID2 &&
				unpackInt16((uint8_t*)&header[14]) == BGZF_LEN);
	}
	// end copied from bgzf.c


};

template< typename Alloc = std::allocator<char> >
class basic_bgzf_compressor : public boost::iostreams::symmetric_filter< bgzf_detail::bgzf_compressor_impl, Alloc > {
private:
	typedef bgzf_detail::bgzf_compressor_impl                impl_type;
	typedef boost::iostreams::symmetric_filter<impl_type, Alloc>  base_type;
public:
	typedef typename base_type::char_type               char_type;
	typedef typename base_type::category                category;
	basic_bgzf_compressor() : base_type(bgzf_detail::DEFAULT_BLOCK_SIZE) {}
	void setAddEOFBlock(bool v) { impl_type::setAddEOFBlock(v); }
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
	basic_bgzf_decompressor() : base_type(bgzf_detail::DEFAULT_BLOCK_SIZE) {}
};
typedef basic_bgzf_decompressor<> bgzf_decompressor;

class bgzf_ostream : public boost::iostreams::filtering_ostream
{
public:
	bgzf_ostream(std::ostream &_os, bool addEOFBlock = true) : dest(_os), comp() {
		comp.setAddEOFBlock(addEOFBlock);
		this->push(comp);
		this->push(dest);
        }
	virtual ~bgzf_ostream() {}
private:
	std::ostream &dest;
	bgzf_compressor comp;
};

class bgzf_istream : public boost::iostreams::filtering_istream
{
public:
	bgzf_istream(std::istream &_is) : dest(_is), decomp() {
		this->push(decomp);
		this->push(dest);
        }
	virtual ~bgzf_istream() {}
private:
	std::istream &dest;
	bgzf_decompressor decomp;
};


#endif
