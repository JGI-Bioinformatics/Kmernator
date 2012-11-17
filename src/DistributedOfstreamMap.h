/*
 * DistributedOfstreamMap.h
 *
 *  Created on: Oct 24, 2012
 *      Author: regan
 */

#ifndef DISTRIBUTEDOFSTREAMMAP_H_
#define DISTRIBUTEDOFSTREAMMAP_H_

#include "mpi.h"

#include "config.h"
#include "Options.h"
#include "Log.h"
#include "ReadSet.h"
#include "Utils.h"

/*
 * DistributedOfstreamMap
 */
class DistributedOfstreamMap : public OfstreamMap
{
public:
	static const int64_t WRITE_BLOCK_SIZE = (8 * 1024 * 1024); // Write in 8MB chunks if possible

private:
	mpi::communicator _world;
	std::string _tempPrefix;
	std::string _realOutputPrefix;

	static string getSharedGlobalUnique(mpi::communicator &world) {
		unsigned long unique = 0;
		if (world.rank() == 0) {
			unsigned int seed = (unsigned int) (((unsigned long) (MPI_Wtime()*1000)) & 0xffffffff);
			unique = LongRand::rand(seed);
			LOG_DEBUG_OPTIONAL(1, true, "getSharedGlobalUnique(): " << unique);
		}
		MPI_Bcast(&unique, 1, MPI_LONG_LONG_INT, 0, world);
		return boost::lexical_cast<string>( unique );
	}
	static string &getGlobalTempDir() {
		static string tempDir;
		return tempDir;
	}
	static string setGlobalTempDir(mpi::communicator &world, std::string tempPath) {
		getGlobalTempDir() = tempPath + "/" + getSharedGlobalUnique(world) + "/";
		mkdir(getGlobalTempDir().c_str(), 0777); // all ranks need to make in case FS is local
		LOG_DEBUG_OPTIONAL(1, world.rank()==0, "setGlobalTempDir(): " << getGlobalTempDir());
		return getGlobalTempDir();
	}
	static void rmGlobalTempDir() {
		rmdir(getGlobalTempDir().c_str());
		getGlobalTempDir().clear();
	}
	static string &getLocalTempDir() {
		static string tempDir;
		return tempDir;
	}
	static string setLocalTempDir(mpi::communicator &world, std::string tempPath) {
		getLocalTempDir() = DistributedDirectoryManagement::makeRankSubDir(world, setGlobalTempDir(world, tempPath));
		mkdir(getLocalTempDir().c_str(), 0777);
		return getLocalTempDir();
	}
protected:
	virtual void close() {
		LOG_VERBOSE_OPTIONAL(2, _world.rank() == 0, "Concatenating all MPI rank files");
		KeySet keys = getGlobalKeySet();
		if (isBuildInMemory())
			writeGlobalFiles(keys);
		OfstreamMap::close();
		if (!isBuildInMemory())
			concatenateMPI(keys);
	}

public:
	DistributedOfstreamMap(mpi::communicator &world, std::string outputFilePathPrefix = Options::getOptions().getOutputFile(), std::string suffix = FormatOutput::getDefaultSuffix(), std::string tempPath = Options::getOptions().getTmpDir())
	:  OfstreamMap(setLocalTempDir(world, tempPath)+"tmp", suffix), _world(world), _tempPrefix(), _realOutputPrefix(outputFilePathPrefix) {

		_tempPrefix = OfstreamMap::getOutputPrefix();
		LOG_DEBUG(3, "DistributedOfstreamMap(world, " << outputFilePathPrefix << ", " << suffix << ", " << getLocalTempDir() << " (from " << tempPath <<") )");
		setBuildInMemory(Options::getOptions().getBuildOutputInMemory());
	}

	~DistributedOfstreamMap() {
		LOG_DEBUG_OPTIONAL(2, _world.rank() == 0, "~DistributedOfstreamMap()");
		this->clear();
		DistributedDirectoryManagement::rmRankSubDir(_world, getGlobalTempDir());
		rmGlobalTempDir(); // all ranks need to rmdir in case it is local
	}

	virtual std::string getRank() const {
		return std::string("--MPIRANK-") + boost::lexical_cast<std::string>(_world.rank());
	}
	virtual void clear() {
		LOG_DEBUG_OPTIONAL(2, true, "DistributedOfstreamMap::clear()");
		this->close();
		_clear();
	}
	virtual std::string getRealFilePath(std::string key) const {
		return _realOutputPrefix + key + getSuffix();
	}
	// gets global keys to rank0.  All other ranks may have partial set...
	KeySet getGlobalKeySet() {
		LOG_DEBUG_OPTIONAL(2, _world.rank()==0, "Calling DistributedOfstreamMap::getGlobalKeySet()");

		// Send all filenames (minus Rank) to master
		KeySet keys = getKeySet();

		if (_world.rank() != 0) {
			_world.send(0, 0, keys);
		} else {
			for(int i = 1 ; i < _world.size() ; i++) {
				KeySet newKeys;
				_world.recv(i, 0, newKeys);
				keys.insert(newKeys.begin(), newKeys.end());
			}
			LOG_DEBUG_OPTIONAL(2, true, "getGlobalKeySet(): Collectively writing " << keys.size() << " files");
			for(KeySet::iterator it = keys.begin(); it != keys.end(); it++)
				LOG_DEBUG(3, "File key: " << *it);
		}

		return keys;
	}
	void writeGlobalFiles(KeySet &keys) {
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Calling DistributedOfstreamMap::writeGlobalFiles()");
		assert(isBuildInMemory());

		int size = _world.size();
		MPI_Comm world = _world;

		// synchronize all files
		int numFiles = keys.size();
		mpi::broadcast(_world, numFiles, 0);

		KeySet::iterator itF = keys.begin();
		for(int fileNum = 0; fileNum < numFiles; fileNum++) {
			std::string key;
			if (_world.rank() == 0) {
				assert(itF != keys.end());
				key = *(itF++);
			}
			mpi::broadcast(_world, key, 0);
			std::string fullPath = getRealFilePath(key);
			LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0, "writeGlobalFiles(): Collectively writing: " << fullPath);

			std::string contents;

			Iterator it = this->_map->find(key);

			if (it == this->_map->end()) {
				LOG_DEBUG_OPTIONAL(1, true, "Could not find " << key << " in DistributedOfstreamMap, contributing 0 bytes");
			} else {
				assert(it->second.isStringStream());
				contents = it->second.getFinalString();
				assert(it->second.empty()); // File must be closed / in a string already
			}

			int64_t mySize = contents.length();
			int64_t sendPos[size], recvPos[size], totalSize = 0, myStart = 0;
			for(int i = 0; i < size; i++)
				sendPos[i] = _world.rank() == i ? mySize : 0;
			MPI_Allreduce(&sendPos, &recvPos, size, MPI_LONG_LONG_INT, MPI_SUM, _world);
			for(int i = 0; i < size; i++) {
				if (_world.rank() == i)
					myStart = totalSize;
				totalSize += recvPos[i];
			}

			LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Opening " << fullPath);
			int err;
			MPI_File ourFile;
			err = MPI_File_open(world, const_cast<char*>(fullPath.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ourFile);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not open " << fullPath << " collectively");
			}
			err = MPI_File_set_size(ourFile, totalSize);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not set the size for " << fullPath << " to " << totalSize);
			}
			LOG_DEBUG(2, "Writing " << mySize << " at " << myStart << " to " << fullPath);
			char *data = const_cast<char*>(contents.data());
			int64_t offset = 0;
			int64_t maxwrite = 0xf000000; // keep writes to less than max int size at a time to avoid MPI overflows
			int64_t bytesToEndofBlock = WRITE_BLOCK_SIZE - (myStart % WRITE_BLOCK_SIZE);
			bytesToEndofBlock = std::min(bytesToEndofBlock, mySize - offset);
			bytesToEndofBlock = bytesToEndofBlock == 0 ? maxwrite : bytesToEndofBlock;
			while (offset < mySize) {
				int64_t thisWriteSize = std::min(maxwrite, mySize - offset);
				thisWriteSize = std::min(bytesToEndofBlock, thisWriteSize);
				MPI_File_write_at(ourFile, myStart+offset, data+offset, thisWriteSize, MPI_BYTE, MPI_STATUS_IGNORE);
				offset += thisWriteSize;
				bytesToEndofBlock = maxwrite;
			}
			LOG_DEBUG_OPTIONAL(1, _world.rank()==0, "Closing " << fullPath);
			MPI_File_close(&ourFile);
		}

	}
	void concatenateMPI(KeySet &keys) {
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Calling DistributedOfstreamMap::concatenateMPI()");

		// synchronize all files
		int numFiles = keys.size();
		mpi::broadcast(_world, numFiles, 0);

		KeySet::iterator itF = keys.begin();
		for(int fileNum = 0; fileNum < numFiles; fileNum++) {
			std::string key;
			if (_world.rank() == 0) {
				assert(itF != keys.end());
				key = *(itF++);
			}
			mpi::broadcast(_world, key, 0);
			std::string fullPath = getRealFilePath(key);
			LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0, "concatenateMPI(): Collectively writing: " << fullPath);

			Iterator it = this->_map->find(key);
			std::string myFilePath;
			if (it == this->_map->end()) {
				LOG_DEBUG_OPTIONAL(1, true, "Could not find myFilePath " << myFilePath << " (" << key << ") in DistributedOfstreamMap");
			} else {
				myFilePath = getFilePath(key);
				assert(it->second.empty()); // File must be closed already
			}
			mergeFiles(_world, myFilePath, fullPath, true);
		}
	}

	static void mergeFiles(mpi::communicator &world, std::string rankFile, std::string globalFile, bool unlinkAfter = false) {
		MPI_Offset mySize = 0;
		char *buf[2];
		int bufSize = WRITE_BLOCK_SIZE;
		buf[0] = new char[bufSize];
		buf[1] = new char[bufSize];
		int bufId = 0;
		MPI_Info info(MPI_INFO_NULL);

		int rank = world.rank();
		int size = world.size();
		int err;
		bool isOpen = false;
		MPI_File myFile;
		if (!rankFile.empty()) {
			err = MPI_File_open(MPI_COMM_SELF, const_cast<char*>(rankFile.c_str()), MPI_MODE_RDONLY, info, &myFile);
			if (err == MPI_SUCCESS) {
				isOpen = true;
				err = MPI_File_get_size(myFile, &mySize);
				if (err != MPI_SUCCESS) {
					LOG_WARN(1, "Could not get the size of " << rankFile << " setting to 0 bytes");
					mySize = 0;
				}
			} else {
				err = MPI_File_close(&myFile);
				LOG_DEBUG_OPTIONAL(1, true, "Could not open " << rankFile << " (for " << globalFile << ") myself.  Merging 0 bytes");
				isOpen = false;
				rankFile.clear();
			}
		} else {
			LOG_DEBUG_OPTIONAL(1, true, "No myFile to merge");
		}

		int64_t sendPos[size], recvPos[size], totalSize = 0, myStart = 0, myPos = 0;
		for(int i = 0; i < size; i++)
			sendPos[i] = (rank == i) ? mySize : 0;
		MPI_Allreduce(&sendPos, &recvPos, size, MPI_LONG_LONG_INT, MPI_SUM, world);
		for(int i = 0; i < size; i++) {
			if (rank == i)
				myStart = totalSize;
			totalSize += recvPos[i];
		}

		LOG_DEBUG_OPTIONAL(1, true, "Writing to '" << globalFile << "' at " << myStart << " for " << mySize << " total: "<< totalSize << " bytes");

		MPI_File ourFile;
		err = MPI_File_open(world, const_cast<char*>(globalFile.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &ourFile);
		if (err != MPI_SUCCESS) {
			LOG_THROW("Could not open " << globalFile << " collectively");
		}
		err = MPI_File_set_size(ourFile, totalSize);
		if (err != MPI_SUCCESS) {
			LOG_THROW("Could not set the size for " << globalFile << " to " << totalSize);
		}

		MPI_Status status, writeStatus;
		MPI_Request writeRequest = MPI_REQUEST_NULL;
		myPos = myStart;
		MPI_Offset bytesToEndofBlock = WRITE_BLOCK_SIZE - (myStart % WRITE_BLOCK_SIZE);
		int lastWrite = 0;
		bytesToEndofBlock = std::min(bytesToEndofBlock, mySize);
		while (isOpen && myPos < myStart + mySize) {
			err = MPI_File_read(myFile, buf[bufId % 2], bytesToEndofBlock == 0 ? bufSize : bytesToEndofBlock, MPI_BYTE, &status);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not read from " << rankFile);
			}
			int sendBytes;
			err = MPI_Get_count(&status, MPI_BYTE, &sendBytes);
			LOG_DEBUG_OPTIONAL(3, true, "Read " << sendBytes << " from " << rankFile);
			if (sendBytes == 0 || err != MPI_SUCCESS)
				break;

			err = MPI_Wait(&writeRequest, &writeStatus);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not wait for write of " << globalFile);
			}
			if (lastWrite > 0) {
				int x;
				err = MPI_Get_count(&writeStatus, MPI_BYTE, &x);
				if (err != MPI_SUCCESS || x != lastWrite) {
					LOG_THROW("Last write was for a different size than requested! " << lastWrite << " vs " << x);
				}
			}
			err = MPI_File_iwrite_at(ourFile, myPos, buf[bufId % 2], sendBytes, MPI_BYTE, &writeRequest);
			lastWrite = sendBytes;
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not write to " << globalFile);
			}
			LOG_DEBUG_OPTIONAL(3, true, "Writing at " << myPos << " for " << sendBytes << " to " << globalFile);
			bufId++;
			myPos += sendBytes;
			bytesToEndofBlock = 0;
		}
		assert(myPos == myStart + mySize);
		err = MPI_Wait(&writeRequest, &writeStatus);
		if (err != MPI_SUCCESS)
			LOG_THROW("Error waiting for write request " << rankFile << " to " << globalFile);
		if (lastWrite > 0) {
			int x;
			err = MPI_Get_count(&writeStatus, MPI_BYTE, &x);
			if (err != MPI_SUCCESS || x != lastWrite) {
				LOG_THROW("Last write was for a different size than requested! " << lastWrite << " vs " << x);
			}
		}

		err = MPI_File_close(&ourFile);
		if (err != MPI_SUCCESS)
			LOG_THROW("Error closing for global file: " << globalFile);

		if (isOpen)
			err = MPI_File_close(&myFile);
		if (isOpen && err != MPI_SUCCESS)
			LOG_THROW("Error closing for rankfile file: " << rankFile);

		delete [] buf[0];
		delete [] buf[1];

		if (isOpen && unlinkAfter)
			unlink(rankFile.c_str());
		LOG_DEBUG_OPTIONAL(1, true, "Finished with temporary: " << rankFile);
	}

	static std::string writeGlobalReadSet(mpi::communicator &world, const ReadSet &readSet, std::string outputFile = Options::getOptions().getOutputFile(), std::string suffix = "", FormatOutput format = FormatOutput::Fasta(), bool trimmed = true)
	{
		DistributedOfstreamMap om(world, outputFile, suffix);
		om.setBuildInMemory();
		std::string fileKey = "";
		readSet.writeAll(om.getOfstream(fileKey), format, trimmed);
		om.clear();
		std::string filename = om.getRealFilePath(fileKey);

		return filename;
	}

};


#endif /* DISTRIBUTEDOFSTREAMMAP_H_ */
