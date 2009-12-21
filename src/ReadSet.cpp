// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.cpp,v 1.12 2009-12-21 22:04:38 regan Exp $
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
 

private:

  class SequenceStreamParser
  {
    istream  *    _stream;
    unsigned long   _line;    
    char          _marker;

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
    	return buffer;
    }
    inline string &nextLine () {
         nextLine(_lineBuffer);
        return _lineBuffer;
    }
    
    SequenceStreamParser(istream &stream,char marker) : _stream(&stream), _line(0),_marker(marker) {  }

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
  };


  class FastaStreamParser : public SequenceStreamParser
  {
    static const char fasta_marker = '>';
    
  public:

    FastaStreamParser(istream &s) : SequenceStreamParser(s,fasta_marker) {    }
  
    virtual string &getBases()
    {
      string &bases = getBasesOrQuals();
      _qualsBuffer.assign(bases.length(),255); // TODO make 255 a const ??
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

void ReadSet::appendFasta(string fastaFilePath,string qualFilePath)
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
   appendFasta(fastaFilePath);
}



ReadSetSizeType ReadSet::getSize()
{
  return _reads.size();
}

Read &ReadSet::getRead(ReadSetSizeType index)
{
  return _reads[index];
}

//
// $Log: ReadSet.cpp,v $
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

