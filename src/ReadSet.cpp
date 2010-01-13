// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.cpp,v 1.20 2010-01-13 23:47:44 regan Exp $
//

#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>

#include <locale>
#include <algorithm>

#include "ReadSet.h"

using  namespace std;

class ReadFileReader
{
private:
    class SequenceStreamParser;
    
    SequenceStreamParser  *_parser;
    string _path; 
    ifstream _ifs;
    ifstream _qs;
    istringstream _iss;

public:

   ReadFileReader(string fastaFilePath,string qualFilePath):
      _parser(NULL),
      _path(fastaFilePath)
   {
      _ifs.open(fastaFilePath.c_str());

      if (_ifs.fail())
        throw  runtime_error("Could not open : " + fastaFilePath) ;

      if (!qualFilePath.empty())
      {
        _qs.open (qualFilePath.c_str());
        if (_qs.fail())
            throw  runtime_error("Could not open : " + qualFilePath) ;
            
        _parser = new FastaQualStreamParser(_ifs,_qs);
      } else {
        // test for an implicit qual file
        _qs.open( (fastaFilePath + ".qual").c_str() );
        if (! _qs.fail()) {
          _parser = new FastaQualStreamParser(_ifs,_qs);
	  } else if ( _ifs.peek() == '@')
          _parser = new FastqStreamParser(_ifs);
        else
          _parser = new FastaStreamParser(_ifs);
      }
      
    }
    
    ReadFileReader(string &fasta) : _iss(fasta) {
    	_parser = new FastaStreamParser(_iss);
    }
    
    ~ReadFileReader()
    {
      _ifs.close();
      _qs.close();
      if (_parser)
         delete _parser;
    }

    bool nextRead(string &name,string &bases, string &quals)
    {
       try {
         name = _parser->getName();
         if (name.empty())
            return false;

        bases = _parser->getBases();
        std::transform(bases.begin(),bases.end(),bases.begin(),::toupper);

        quals = _parser->getQuals();

        if (quals.length() != bases.length())
          if ( !(quals.length() == 1 && quals[0] == Read::REF_QUAL) )
            throw  runtime_error("Number of bases and quals not equal");
        return true;
       }
       
       catch (runtime_error &e)
       {
          stringstream error;
          error << e.what() << " in file '" << _path << "' at line " << _parser->lineNumber();
          throw runtime_error(error.str());
       } 
    }
    
    unsigned long getFileSize()
    {
    	ifstream::streampos current = _ifs.tellg();
    	_ifs.seekg(0, ios_base::end);
    	unsigned long size = _ifs.tellg();
    	_ifs.seekg(current);
    	return size;
    }
 
    unsigned long getBlockSize(unsigned int numThreads)
    {
        return getFileSize() / numThreads;
    }
    
    inline unsigned long getPos()
    {
    	return _parser->getPos();
    }
    
    void seekToNextRecord(unsigned long minimumPos) {
    	_parser->seekToNextRecord(minimumPos);
    }
    int getType() {
    	return _parser->getType();
    }

private:

  class SequenceStreamParser
  {
    istream  *    _stream;
    unsigned long   _line;    
    char          _marker;
    unsigned long   _pos;

  protected:
    string _nameBuffer;
    string _basesBuffer;
    string _qualsBuffer;
    string _lineBuffer;
    
  public:
    inline string &nextLine(string &buffer) {
        _line++;
        buffer.clear();
    	getline(*_stream, buffer);
    	_pos += buffer.length()+1;
    	return buffer;
    }
    inline string &nextLine () {
         nextLine(_lineBuffer);
        return _lineBuffer;
    }
    
    SequenceStreamParser(istream &stream,char marker) : _stream(&stream), _line(0),_marker(marker),_pos(0) {  }

    virtual ~SequenceStreamParser() {  }
    
    unsigned long lineNumber() { return _line;}
    bool endOfStream()         { return _stream->eof();}
    int  peek()                { return _stream->peek();}

    virtual string &getName()
    {
        nextLine(_nameBuffer);
      
        while (_nameBuffer.length() == 0) // skip empty lines at end of stream
        {
          if (endOfStream()) {
          	_nameBuffer.clear();
          	return _nameBuffer;
          }
          nextLine(_nameBuffer);
        }

        if (_nameBuffer[0] != _marker)
          throw  runtime_error((string("Missing '") + _marker + "'").c_str());

//         char *space = strchr(name,' ');
//         if (space != NULL)
//             space = '\0';
//         char *tab = strchr(name,'\t');
//         if (tab !=NULL)
//             tab = '\0';

      return _nameBuffer.erase(0,1);
    }
    
    virtual string &getBases() = 0;
    virtual string &getQuals() = 0;
    virtual int getType() = 0;
    
    virtual unsigned long getPos() {
    	return _pos;
    }
    virtual void seekToNextRecord(unsigned long minimumPos)
    {
    	// get to the first line after the pos
    	if (minimumPos > 0) {
    	  _pos = minimumPos -1;
    	  _stream->seekg(_pos);
    	  if (endOfStream())
    	    return;
    	  if (peek() == '\n') {
    	    _pos = minimumPos;
    	    _stream->seekg(_pos);
    	  }
    	  else
    	    nextLine();
    	}
    	else {
    	  _pos = 0;
    	  _stream->seekg(0);
    	}
    	while(_stream->peek() != _marker && ! endOfStream())
    	  nextLine();
    	
    }
    
  };


  class FastqStreamParser : public SequenceStreamParser
  {

  public:
    FastqStreamParser(istream &s) : SequenceStreamParser(s,'@') { }

    string &getBases()
    {
      nextLine(_basesBuffer);

      if (_basesBuffer.empty()  || _basesBuffer.length() == _basesBuffer.max_size())
          throw  runtime_error("Missing or too many bases");

      string &qualName = nextLine();
      if (qualName.empty() || qualName[0] != '+')
          throw  runtime_error("Missing '+'");

      return _basesBuffer;
    }

    string &getQuals() { return nextLine(_qualsBuffer); }
    int getType() { return 0; }
  };


  class FastaStreamParser : public SequenceStreamParser
  {
    static const char fasta_marker = '>';
    
  public:

    FastaStreamParser(istream &s) : SequenceStreamParser(s,fasta_marker) {    }
  
    virtual string &getBases()
    {
      string &bases = getBasesOrQuals();
      _qualsBuffer.assign(bases.length(), Read::REF_QUAL); 
      return bases;
    }
    
    virtual string &getQuals() {  return _qualsBuffer; }
    
    virtual string &getBasesOrQuals()
    {
      _basesBuffer.clear();
      while (!endOfStream())
      {
         nextLine();
          
          if (_lineBuffer.empty())
            break;
            
         _basesBuffer += _lineBuffer;
         if ( _basesBuffer.size() == _basesBuffer.max_size())
            throw  runtime_error("Sequence/Qual too large to read");
          
          if (peek() == fasta_marker)
            break;
      }
      return _basesBuffer;
    }
    int getType() { return 1; }
  };


  class FastaQualStreamParser : public FastaStreamParser
  {
    FastaStreamParser _qualParser;  // odd, but it works

  public:
    FastaQualStreamParser(istream &fastaStream,istream &qualStream) : FastaStreamParser(fastaStream) , _qualParser(qualStream)
    {   }

    string &getName() {
      string &qualName = _qualParser.getName();
      string &name     =  FastaStreamParser::getName();
      if (name != qualName)
        throw runtime_error("fasta and quals have different names");
      return name;
    }

    string &getBases() {
      return getBasesOrQuals( );
    }

    string &getQuals() {
      // odd, but it works
      string &qualValues = _qualParser.getBasesOrQuals(); 
      istringstream ss(qualValues);

      _qualsBuffer.clear();
      while (!ss.eof()) {
        int qVal;
        ss >> qVal;
        if (ss.fail())
          break;
        qVal += 64;
        _qualsBuffer.push_back(qVal);
      }
      return _qualsBuffer;
    }
    int getType() { return 1; }
  };
};



/*
std::ifstream::pos_type _fileSize(ifstream &f)
{
  std::ifstream::pos_type current = f.tellg();
  f.seekg(0, std::ios_base::beg);
  std::ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, std::ios_base::end);
  std::ifstream::pos_type size = f.tellg() - begin_pos;
  f.seekg(current,std::ios_base::beg);

  return size;
}*/


ReadSet::ReadSet():
_baseCount(0)
{

}

ReadSet::~ReadSet()
{
   
}

void ReadSet::addRead(Read &read) {
	_reads.push_back( read );
	_baseCount += read.getLength();
}

void ReadSet::appendAnyFile(string filePath, string filePath2)
{
	ReadFileReader reader(filePath, filePath2);
	switch (reader.getType()) {
		case 0: appendFastq(filePath); break;
		case 1: appendFasta(reader);
	}
}

void ReadSet::appendAllFiles(Options::FileListType &files)
{
	
    #ifdef _USE_OPENMP
	ReadSet myReads[ files.size() ];
	omp_set_nested(1);
	#pragma omp parallel for schedule(dynamic)
	#endif
	for(long i = 0; i < (long) files.size(); i++) {
		#ifdef _USE_OPENMP
		#pragma omp critical
		#endif
		{ std::cerr << "reading " << files[i] << std::endl; }
		
		#ifdef _USE_OPENMP
		// append int this thread's ReadSet buffer (note: line continues)
		myReads[ i ].		
		#endif
		appendAnyFile(files[i]);
		
		#ifdef _USE_OPENMP
		#pragma omp critical
		#endif
		{ std::cerr << "finished reading " << files[i] << std::endl; }
	}
	#ifdef _USE_OPENMP
	omp_set_nested(OMP_NESTED_DEFAULT);
	{ std::cerr << "concatenating ReadSet buffers" << std::endl; }
	for(int i = 0; i< (long) files.size(); i++)
	  append(myReads[i]);
	#endif
	
}

void ReadSet::append(ReadSet &reads)
{
    unsigned long oldSize = _reads.size();
    unsigned long newSize = oldSize + reads._reads.size();
    _reads.resize(newSize);
    #ifdef _USE_OPENMP
	#pragma omp parallel for
	#endif
	for(long i = 0; i < (long) reads._reads.size(); i++)
	  _reads[ oldSize + i ] = reads._reads[i];
	
	unsigned long oldPairSize = _pairs.size();
	unsigned long newPairSize = oldPairSize + reads._pairs.size();
	_pairs.resize(newPairSize);
	
    #ifdef _USE_OPENMP
	#pragma omp parallel for
	#endif
	for(long i = 0; i < (long) reads._pairs.size(); i++) {
	  Pair &tmp =  reads._pairs[i];
	  _pairs[ oldPairSize + i ] = Pair( (tmp.read1 == MAX_READ_IDX ? MAX_READ_IDX : tmp.read1+oldSize), (tmp.read2 == MAX_READ_IDX ? MAX_READ_IDX : tmp.read2+oldSize) );
	}
	_baseCount += reads._baseCount;
		  
}

void ReadSet::appendFasta(string fastaFilePath,string qualFilePath)
{
    ReadFileReader reader(fastaFilePath,qualFilePath);
    appendFasta(reader);
}
void ReadSet::appendFasta(ReadFileReader &reader)
{
    string name,bases,quals;
    while (reader.nextRead(name,bases,quals))
    {
       Read read(name, bases, quals);
       addRead( read );
    }
}

void ReadSet::appendFastaFile(string &str)
{
	ReadFileReader reader(str);
	appendFasta(reader);
}


#ifdef _USE_OPENMP

void ReadSet::appendFastqBlockedOMP(string fastaFilePath,string qualFilePath)
{
  
  unsigned long startIdx = _reads.size();
  unsigned long blockSize = 0;
  unsigned long tmpBaseCount = 0;
  
  // set OMP variables
  omp_set_nested(1);
  int numThreads = omp_get_max_threads();
  if (numThreads > 8)
    numThreads = 8;
  unsigned long numReads[ numThreads ];
  unsigned long seekPos[ numThreads ];
  for(int i = 0; i < numThreads ; i++)
    numReads[i] = 0;
  
  #pragma omp parallel reduction(+:tmpBaseCount) num_threads(numThreads)
  {
  	
  	if (omp_get_num_threads() != numThreads)
  	  throw "OMP thread count discrepancy!";
  	  
  	ReadSet myReads;
    ReadFileReader reader(fastaFilePath,qualFilePath);
    
    #pragma omp single
	{
		blockSize = reader.getBlockSize(numThreads);
		if (blockSize < 100)
		  blockSize = 100;
		std::cerr << "Reading " << fastaFilePath << " with " << numThreads << " threads" << std::endl;
	}
	reader.seekToNextRecord( blockSize * omp_get_thread_num() );
	seekPos[ omp_get_thread_num() ] = reader.getPos();

    string name,bases,quals;
    bool hasNext = omp_get_thread_num() < (long) numThreads;
    if (hasNext)
      hasNext = (reader.getPos() < blockSize * (omp_get_thread_num() +1));
    
    if (!hasNext)
      seekPos[ omp_get_thread_num() ] = 0;
    
    #pragma omp barrier
	if (hasNext && omp_get_thread_num() != 0 && seekPos[omp_get_thread_num()] == seekPos[omp_get_thread_num()-1])
	  hasNext = false;

	//#pragma omp critical
	//{ std::cerr << omp_get_thread_num() << " " << blockSize << " seeked to " << reader.getPos() << " will read: " << hasNext << std::endl; }
    
    while (hasNext && reader.nextRead(name,bases,quals))
    {
       Read read(name, bases, quals);
       myReads.addRead( read );
       tmpBaseCount += read.getLength();
       // todo look into getPos() to optimize
       hasNext = (reader.getPos() < blockSize * (omp_get_thread_num() +1));
    }
    numReads[ omp_get_thread_num() ] = myReads.getSize();

	//#pragma omp critical
	//{ std::cerr << omp_get_thread_num() << " " << blockSize << " finished at " << reader.getPos() << " with " << myReads.getSize() << " " << numReads[ omp_get_thread_num() ]<< std::endl; }
    
    #pragma omp barrier
	
	#pragma omp single
    {
    	unsigned long newReads = 0; 
    	for(int i=0; i < numThreads; i++) {
    	  unsigned long tmp = newReads;
    	  newReads += numReads[i];
    	  numReads[i] = tmp;
    	}
        
        _reads.resize(startIdx + newReads);
    }
    
  	for(unsigned long j = 0; j <  myReads.getSize(); j++)
  	  _reads[startIdx + numReads[ omp_get_thread_num() ] + j] = myReads.getRead(j);
  }
  _baseCount += tmpBaseCount;
  
  // reset omp variables
  omp_set_nested(OMP_NESTED_DEFAULT);
  
}

// not used as it speeds up only marginally
void ReadSet::appendFastqBatchedOMP(string fastaFilePath,string qualFilePath)
{
    ReadFileReader reader(fastaFilePath,qualFilePath);
    
    long batchSize = 100000;
    long batchIdx = 0;
    std::string names[batchSize];
    std::string bases[batchSize];
    std::string quals[batchSize];
    unsigned long tmpBaseCount = 0;
    bool hasNext = true;
    while (hasNext)
    {
    	// single threaded reading...
        for(long idx = 0 ; idx < batchSize; idx++)
          if( ! reader.nextRead(names[idx],bases[idx],quals[idx]) ) {
          	hasNext = false;
          	batchSize = idx;
          	break;
          }
        
        // allocate space
        _reads.resize( batchIdx + batchSize );
        
        
        #pragma omp parallel for reduction(+:tmpBaseCount)
		for(long idx = 0; idx < batchSize ; idx++) {
	      Read read(names[idx], bases[idx], quals[idx]);
		  _reads[batchIdx+idx] = read;
		  tmpBaseCount += read.getLength();
		}
		  
		batchIdx += batchSize;
    }
    _baseCount += tmpBaseCount;
    
}

void ReadSet::appendFastq(string fastaFilePath)
{
   appendFastqBlockedOMP(fastaFilePath);
}

#else // _USE_OPENMP

void ReadSet::appendFastq(string fastaFilePath)
{
   appendFasta(fastaFilePath);
}

#endif // _USE_OPENMP

std::string _commonName(const std::string &readName)
{
	return readName.substr(0, readName.length() - 1);
}
int _readNum(const std::string &readName)
{
	int retVal = 0;
	int len = readName.length();
	char c = readName[ len - 1];
	switch (c) {
		case '1': if (readName[len-2] != '/')  break;
		case 'A':
		case 'F': retVal = 1; break;
		case '2': if (readName[len-2] != '/')  break;
		case 'B':
		case 'R': retVal = 2; break;
	}
	return retVal;
}

ReadSet::ReadSetSizeType ReadSet::identifyPairs()
{
  ReadSetSizeType size = getSize();
  
  // find the last Pair to be identified
  ReadSetSizeType readIdx = 0; 
  if (_pairs.size() > 0) {
    Pair &lastPair = _pairs[_pairs.size() - 1];
    if (isValidRead(lastPair.read1))
      readIdx = lastPair.read1;
    if (isValidRead(lastPair.read2) && readIdx < lastPair.read2)
      readIdx = lastPair.read2;
    readIdx++; // start with the next read
  }
  if (size <= readIdx)
    return _pairs.size();
  
  boost::unordered_map< std::string, ReadSetSizeType > unmatchedNames;
  boost::unordered_map< std::string, ReadSetSizeType >::iterator unmatchedIt;
  
  for(ReadSetSizeType i = 0 ; i < _pairs.size() ; i++) {
  	Pair &pair = _pairs[i];
  	if ( (pair.read1 == MAX_READ_IDX || pair.read2 == MAX_READ_IDX) && pair.read1 != pair.read2) {
  	  ReadSetSizeType idx = pair.lesser();
  	  unmatchedNames[ _commonName( _reads[idx].getName()) ] = idx;	
  	}
  }
  
  std::string lastName, common;
  int readNum = 0;
  
  while (readIdx < size) {
  
    string name = getRead(readIdx).getName();
    readNum = _readNum(name); 
    common = _commonName(name);
    
    unmatchedIt = unmatchedNames.find( common );
    if (unmatchedIt != unmatchedNames.end()) {
      Pair &test = _pairs[ unmatchedIt->second ];
      if (readNum == 2) {
      	if (test.read2 != MAX_READ_IDX)
      	  throw;
      	_pairs[ unmatchedIt->second ].read2 = readIdx;
      } else {
      	if (test.read1 != MAX_READ_IDX)
      	  throw;
        _pairs[ unmatchedIt->second ].read1 = readIdx;
      }
      unmatchedNames.erase(unmatchedIt);
    } else {
      Pair pair;
      if (readNum == 2)
        pair.read2 = readIdx;
      else
        pair.read1 = readIdx;
      if (readNum > 0)
        unmatchedNames[ common ] = _pairs.size();
      _pairs.push_back(pair);
    }
    	
    readIdx++;
  }
  
  return _pairs.size();
}

//
// $Log: ReadSet.cpp,v $
// Revision 1.20  2010-01-13 23:47:44  regan
// made const class modifications
// fixed identify pairs
//
// Revision 1.19  2010-01-13 00:26:49  regan
// fixed some parallelism
// started pair indentification
//
// Revision 1.18  2010-01-08 06:23:07  regan
// fixed to support nested parallel reading of files
//
// Revision 1.17  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.16  2010-01-05 06:44:39  regan
// fixed warnings
//
// Revision 1.15  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.14  2009-12-23 07:16:52  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.13  2009-12-22 18:31:41  regan
// parallelized reading fastq if openmp is enabled
//
// Revision 1.12  2009-12-21 22:04:38  regan
// minor optimization
//
// Revision 1.11  2009-12-21 07:54:18  regan
// minor parallelization of reading files step
//
// Revision 1.10  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.9  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.7  2009-11-04 20:14:46  cfurman
// added conversion to uppercase
//
// Revision 1.6  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.5  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.4  2009-10-23 07:06:59  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.3  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//

