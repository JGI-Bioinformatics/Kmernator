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
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/lexical_cast.hpp>

template< std::streamsize BUFFER_SIZE = 131072 >
class memory_buffer_detail {
public:
	class Buffer {
	public:
		Buffer() : readOffset(0), writeOffset(0), buf(new char[BUFFER_SIZE]) {}
		~Buffer() { delete [] buf; }
		volatile std::streamsize readOffset, writeOffset;
		char * buf;
	};
	typedef boost::shared_ptr< Buffer > BufferType;

	memory_buffer_detail(bool _setExplicitWriteClose = false)
	: writeCount(0), readCount(0),
	  setExplicitWriteClose(_setExplicitWriteClose), writeClosed(false), blocked(false) {
		addBuffer();
		LOG_DEBUG_OPTIONAL(2, true, "memory_buffer_detail(): " << this);
	}
	virtual ~memory_buffer_detail() {
		LOG_DEBUG_OPTIONAL(2, true, "~memory_buffer_detail(): " << this << " readCount: " << readCount << " writeCount: " << writeCount << " buffers:" << buffers.size());
		clear();
	}
	// returns true if reading while write is still open and buffer is empty
	bool isBlocked() {
		return blocked;
	}

	// signals write is complete so last reads can complete too
	void writeClose() {
		writeClosed = true;
		if (setExplicitWriteClose == true) {
			boost::mutex::scoped_lock mylock(writeCloseLock);
			readReady.notify_one();
		}
		LOG_DEBUG_OPTIONAL(2, true, "memory_buffer_detail::writeClose(): " << this << " at " << writeCount);
	}

	// never returns less than n
	std::streamsize write(const char* s, std::streamsize n) {
		assert(!writeClosed);
		std::streamsize wrote = 0, remaining = n;
		while (remaining > 0) {
			wrote = writeSome(s, remaining);
			s+=wrote;
			remaining -= wrote;
		}
		return n;
	}
	// returning less than n signifies EOF
	// so only do so if writing has stopped
	std::streamsize read(char *s, std::streamsize n) {
		std::streamsize totalReadSize = 0, remaining = n;
		while (remaining > 0) {
			std::streamsize readSize = readSome(s, remaining);
			s+=readSize;
			remaining -= readSize;
			totalReadSize += readSize;
			if (readSize == 0) {
				if (setExplicitWriteClose && !writeClosed) {
					blocked = true;
					boost::mutex::scoped_lock mylock(writeCloseLock);
					readReady.timed_wait(mylock, boost::get_system_time() + boost::posix_time::milliseconds(500));
				} else {
					break;
				}
			}
		}
		blocked = false;
		return totalReadSize;
	}
	std::streamsize tellp() const {
		return writeCount;
	}
	std::streamsize tellg() const {
		return readCount;
	}

	std::streamsize writeSome(const char* s, std::streamsize n) {
		BufferType b = back;
		volatile std::streamsize &writeOffset = b->writeOffset;
		std::streamsize wrote = std::min(BUFFER_SIZE - writeOffset, n);

		memcpy(b->buf + writeOffset, s, wrote);
		{
			// make copies to prevent race between read & write threads
			std::streamsize tmp = writeOffset + wrote;
			writeOffset = tmp;
		}
		writeCount += wrote;
		if (writeOffset == BUFFER_SIZE)
			addBuffer();
		if (setExplicitWriteClose == true)
			readReady.notify_one();
		return wrote;
	}

	std::streamsize readSome(char *s, std::streamsize n) {
		int readSize;
		BufferType b = front;
		if (b.get() == NULL)
			return 0;
		else {
			std::streamsize writeOffset = b->writeOffset; // make a copy
			volatile std::streamsize &readOffset = b->readOffset;  // direct reference
			readSize = std::min(n, writeOffset - readOffset);
			memcpy(s, b->buf + readOffset, readSize);
			readOffset += readSize;
			readCount += readSize;
			if (readOffset == BUFFER_SIZE || (readOffset == b->writeOffset && writeClosed))
				removeFirstBuffer();
		}
		return readSize;
	}

	void concat(memory_buffer_detail &src) {
		assert(writeClosed && src.writeClosed);
		readCount += src.readCount;
		writeCount += src.writeCount;
		buffers.insert(buffers.end(), src.buffers.begin(), src.buffers.end());
		front = buffers.front();
		back = buffers.back();
		src.clear();
	}

protected:
	void clear() {
		{
			writeCount = readCount = 0;
			buffers.clear();
			front.reset();
			back.reset();
		}
	}

	void addBuffer() {
		BufferType p( new Buffer() );
		{
			buffers.push_back( p );
			if (front.get() != buffers.front().get())
				front = buffers.front();
			back = p;
		}
	}
	void removeFirstBuffer() {
		assert(!buffers.empty());
		{
			buffers.pop_front();
			if (!buffers.empty()) {
				front = buffers.front();
				if (back.get() != buffers.back().get())
					back = buffers.back();
			} else {
				front.reset();
				back.reset();
			}
		}
	}

private:
	std::streamsize writeCount, readCount;
	std::deque< BufferType > buffers;
	bool setExplicitWriteClose, writeClosed, blocked;
	BufferType front, back;
	boost::mutex writeCloseLock;
	boost::condition_variable readReady;
};


template< typename Impl >
class stream_impl_template {
public:
	typedef char char_type;
	typedef boost::iostreams::bidirectional_device_tag category;
	typedef boost::shared_ptr< Impl > ImplPtr;

	stream_impl_template() : _impl( new Impl() ) {}

	template<typename U1>
	stream_impl_template(U1 &u1) : _impl( new Impl(u1) ) {};

	template<typename U1, typename U2>
	stream_impl_template(U1 &u1, U2 &u2) : _impl( new Impl(u1, u2) ) {};

	template<typename U1, typename U2, typename U3>
	stream_impl_template(U1 &u1, U2 &u2, U3 &u3) : _impl( new Impl(u1, u2, u3) ) {};

	~stream_impl_template() {}
	stream_impl_template &operator=(const stream_impl_template &copy) {
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

	template<typename U>
	stream_impl_template &operator<<(U u) {
		std::string s = boost::lexical_cast<std::string>(u);
		write(s.c_str(), s.length());
		return *this;
	}

	// specialized for memory_buffer_detail
	bool isBlocked() {
		return _impl->isBlocked();
	}
	void writeClose() {
		return _impl->writeClose();
	}
	void concat(stream_impl_template &src) {
		_impl->concat(*src._impl);
	}

	class ostream : public boost::iostreams::filtering_ostream
	{
	public:
		ostream(stream_impl_template &_os) {
			this->push(_os);
		}
/*
		template<typename U>
		ostream(U &u, stream_impl_template &_os) {
			this->push(u);
			this->push(_os);
		};
*/
		virtual ~ostream() {}
	};

	class istream : public boost::iostreams::filtering_istream
	{
	public:
		istream(stream_impl_template &_is) {
			this->push(_is);
		}
/*
		template<typename U>
		istream(U &u, stream_impl_template &_is) {
			this->push(u);
			this->push(_is);
		};
*/
		virtual ~istream() {}
	};

private:
	ImplPtr _impl;
};

typedef stream_impl_template< memory_buffer_detail< > > MemoryBuffer;
typedef boost::shared_ptr< MemoryBuffer > MemoryBufferPtr;
typedef std::vector< MemoryBufferPtr > MemoryBufferVector;

#endif // MEMORY_BUFFER_STREAM_H_
