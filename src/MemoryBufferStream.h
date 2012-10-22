/*
 * MemoryBufferStream.h
 *
 *  Created on: Oct 22, 2012
 *      Author: regan
 */

#ifndef MEMORY_BUFFER_STREAM_H_
#define MEMORY_BUFFER_STREAM_H_

#include <string>
#include <iostream>
#include <deque>
#include <algorithm>

#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/shared_ptr.hpp>

template< std::streamsize BUFFER_SIZE = 262144 >
class memory_buffer_detail {
public:
	typedef boost::shared_ptr< char > BufferType;

	memory_buffer_detail() : writeOffset(BUFFER_SIZE), readOffset(0), writeCount(0), readCount(0) {
		addBuffer();
	}
	virtual ~memory_buffer_detail() {
		clear();
	}
	std::streamsize write(const char* s, std::streamsize n) {
		std::streamsize wrote = 0, remaining = n;
		while (remaining > 0) {
			wrote = writeSome(s, remaining);
			s+=wrote;
			remaining -= wrote;
		}
		return n;
	}
	std::streamsize writeSome(const char* s, std::streamsize n) {
		std::streamsize wrote = std::min(BUFFER_SIZE - writeOffset, n);
		memcpy(buffers.back().get() + writeOffset, s, wrote);
		writeOffset += wrote;
		writeCount += wrote;
		if (writeOffset == BUFFER_SIZE)
			addBuffer();
		return wrote;
	}
	std::streamsize read(char *s, std::streamsize n) {
		std::streamsize totalReadSize = 0, remaining = n;
		while (remaining > 0) {
			std::streamsize readSize = readSome(s, remaining);
			s+=readSize;
			remaining -= readSize;
			totalReadSize += readSize;
			if (readSize == 0)
				break;
		}
		return totalReadSize;
	}
	std::streamsize readSome(char *s, std::streamsize n) {
		int readSize;
		if (buffers.size() == 0)
			return 0;
		else if (buffers.size() == 1) {
			// front buffer may not be full
			readSize = std::min(n, (writeOffset - readOffset));
			if (readSize == 0)
				return 0;
		} else {
			// front buffer is full
			readSize = std::min(n, BUFFER_SIZE - readOffset);
		}
		memcpy(s, buffers.front().get() + readOffset, readSize);
		readOffset += readSize;
		readCount += readSize;
		if (readOffset == BUFFER_SIZE)
			removeFirstBuffer();
		return readSize;
	}
	std::streamsize tellp() const {
		return writeCount;
	}
	std::streamsize tellg() const {
		return readCount;
	}

protected:
	void clear() {
		writeOffset = readOffset = writeCount = readCount = 0;
		buffers.clear();
	}
	void addBuffer() {
		assert(writeOffset == BUFFER_SIZE);
		BufferType p( new char[BUFFER_SIZE] );
		buffers.push_back( p );
		writeOffset = 0;
	}
	void removeFirstBuffer() {
		assert(readOffset == BUFFER_SIZE);
		assert(!buffers.empty());
		buffers.pop_front();
		readOffset = 0;
	}
private:
	std::deque< BufferType > buffers;
	std::streamsize writeOffset, readOffset;
	std::streamsize writeCount, readCount;
};

template< std::streamsize BUFFER_SIZE = 262144 >
class memory_buffer {
public:
	typedef char char_type;
	typedef boost::iostreams::bidirectional_device_tag category;
	typedef memory_buffer_detail< BUFFER_SIZE > Impl;
	typedef boost::shared_ptr< Impl > ImplPtr;
	memory_buffer() : _impl( new Impl() ) {}
	memory_buffer(const memory_buffer &copy) {
		_impl = copy._impl;
	}
	~memory_buffer() {}
	memory_buffer &operator=(const memory_buffer &copy) {
		_impl = copy._impl;
		return *this;
	}
	std::streamsize write(const char* s, std::streamsize n) {
		return _impl->write(s, n);
	}
	std::streamsize read(char *s, std::streamsize n) {
		return _impl->read(s, n);
	}
	std::streamsize tellp() const {
		return _impl->tellp();
	}
	std::streamsize tellg() const {
		return _impl->tellg();
	}

	class ostream : public boost::iostreams::filtering_ostream
	{
	public:
		ostream(memory_buffer &_os) {
			this->push(_os);
		}
		virtual ~ostream() {}
	};

	class istream : public boost::iostreams::filtering_istream
	{
	public:
		istream(memory_buffer &_is) {
			this->push(_is);
		}
		virtual ~istream() {}
	};

private:
	ImplPtr _impl;
};

typedef memory_buffer<> MemoryBuffer;

#endif // MEMORY_BUFFER_STREAM_H_
