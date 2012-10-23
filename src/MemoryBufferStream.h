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

template< std::streamsize BUFFER_SIZE = 262144 >
class memory_buffer_detail {
public:
	typedef boost::shared_ptr< char > BufferType;

	memory_buffer_detail(bool _setExplicitWriteClose = false)
	: writeOffset(BUFFER_SIZE), readOffset(0), writeCount(0), readCount(0),
	  writerTid(-1), readerTid(-1), setExplicitWriteClose(_setExplicitWriteClose), writeClosed(false) {
		addBuffer();
	}
	virtual ~memory_buffer_detail() {
		clear();
	}
	void writeClose() {
		writeClosed = true;
		if (setExplicitWriteClose == true)
			readReady.notify_one();
	}

	// never returns less than n
	std::streamsize write(const char* s, std::streamsize n) {
		if (writerTid < 0)
			writerTid = omp_get_thread_num();
		assert(writerTid == omp_get_thread_num());
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
	std::streamsize read(char *s, std::streamsize n) {
		if (readerTid < 0)
			readerTid = omp_get_thread_num();
		assert(readerTid == omp_get_thread_num());
		std::streamsize totalReadSize = 0, remaining = n;
		while (remaining > 0) {
			std::streamsize readSize = readSome(s, remaining);
			s+=readSize;
			remaining -= readSize;
			totalReadSize += readSize;
			if (readSize == 0) {
				if (setExplicitWriteClose && !writeClosed) {
					boost::mutex::scoped_lock mylock(writeCloseLock);
					readReady.wait(mylock);
				} else {
					break;
				}
			}
		}
		return totalReadSize;
	}
	std::streamsize tellp() const {
		return writeCount;
	}
	std::streamsize tellg() const {
		return readCount;
	}

protected:
	void clear() {
		{
			bool locked = false;
			if (needsLock()) { lock.lock(); locked = true;}
			writeOffset = readOffset = writeCount = readCount = 0;
			buffers.clear();
			front.reset();
			back.reset();
			if (locked) lock.unlock();
		}
		writerTid = readerTid = -1;
	}

	std::streamsize writeSome(const char* s, std::streamsize n) {
		std::streamsize wrote = std::min(BUFFER_SIZE - writeOffset, n);
		BufferType b = back;
		memcpy(b.get() + writeOffset, s, wrote);
		#pragma omp atomic
		writeOffset += wrote;
		writeCount += wrote;
		if (writeOffset == BUFFER_SIZE)
			addBuffer();
		if (setExplicitWriteClose == true)
			readReady.notify_one();
		return wrote;
	}

	std::streamsize readSome(char *s, std::streamsize n) {
		int readSize;
		if (front.get() == NULL)
			return 0;
		else {
			BufferType b = front;

			if (b.get() == back.get()) {
				// front buffer may not be full
				bool locked = false;
				if (needsLock()) { lock.lock(); locked = true; }
				if (b.get() == back.get()) {
					readSize = std::min(n, (writeOffset - readOffset));
					if (locked) lock.unlock();
					if (readSize == 0)
						return 0;
				} else {
					// front buffer is full
					if (locked) lock.unlock();
					readSize = std::min(n, BUFFER_SIZE - readOffset);
				}
			} else {
				// front buffer is full
				readSize = std::min(n, BUFFER_SIZE - readOffset);
			}

			memcpy(s, b.get() + readOffset, readSize);
			#pragma omp atomic
			readOffset += readSize;
			readCount += readSize;
			if (readOffset == BUFFER_SIZE)
				removeFirstBuffer();
		}
		return readSize;
	}
	void addBuffer() {
		assert(writeOffset == BUFFER_SIZE);
		BufferType p( new char[BUFFER_SIZE] );
		{
			bool locked = false;
			if (needsLock()) { lock.lock(); locked = true;}
			buffers.push_back( p );
			if (front.get() != buffers.front().get())
				front = buffers.front();
			back = p;
			writeOffset = 0;
			if (locked) lock.unlock();
		}
	}
	void removeFirstBuffer() {
		assert(readOffset == BUFFER_SIZE);
		assert(!buffers.empty());
		{
			bool locked = false;
			if (needsLock()) { lock.lock(); locked = true; }
			buffers.pop_front();
			if (!buffers.empty()) {
				front = buffers.front();
				if (back.get() != buffers.back().get())
					back = buffers.back();
			} else {
				front.reset();
				back.reset();
			}
			readOffset = 0;
			if (locked) lock.unlock();
		}
	}

	bool needsLock() {
		return writerTid != readerTid && writerTid >=0 && readerTid >= 0;
	}

private:
	std::deque< BufferType > buffers;
	std::streamsize writeOffset, readOffset;
	std::streamsize writeCount, readCount;
	int writerTid, readerTid;
	bool setExplicitWriteClose, writeClosed;
	BufferType front, back;
	boost::mutex lock, writeCloseLock;
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

	void writeClose() {
		return _impl->writeClose();
	}

	class ostream : public boost::iostreams::filtering_ostream
	{
	public:
		ostream(stream_impl_template &_os) {
			this->push(_os);
		}
		virtual ~ostream() {}
	};

	class istream : public boost::iostreams::filtering_istream
	{
	public:
		istream(stream_impl_template &_is) {
			this->push(_is);
		}
		virtual ~istream() {}
	};

private:
	ImplPtr _impl;
};

typedef stream_impl_template< memory_buffer_detail< > > MemoryBuffer;

#endif // MEMORY_BUFFER_STREAM_H_
