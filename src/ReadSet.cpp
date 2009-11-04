// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.cpp,v 1.7 2009-11-04 20:14:46 cfurman Exp $
//

#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <locale>
#include <algorithm>

#include "ReadSet.h"

#define MAX_LINE_LENGTH 1024*1024

using  namespace std;


std::ifstream::pos_type _fileSize(ifstream &f)
{
  std::ifstream::pos_type current = f.tellg();
  f.seekg(0, std::ios_base::beg);
  std::ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, std::ios_base::end);
  std::ifstream::pos_type size = f.tellg() - begin_pos;
  f.seekg(current,std::ios_base::beg);

  return size;
}


ReadSet::ReadSet():
_baseCount(0)
{

}

ReadSet::~ReadSet()
{
   
}

void ReadSet::push_back(Read &read) {
	_reads.push_back( read );
	_baseCount += read.getLength();
}


 
class fast_stream_base : public std::ifstream
{
private: 
   std::string _path;
   unsigned long _line;
      
public:
  
  fast_stream_base(std::string path) : _path(path), _line(0)
  { 
     open(_path.c_str());

    if (fail()) {
       throw  std::invalid_argument("Could not open : " + _path) ;
    }     
  }

  virtual ~fast_stream_base()
  {
  }
  
  istream& getline (char* s, streamsize n )
  {
      _line++;
      return std::istream::getline(s,n);
  }


   std::string get_name(char marker)
   {
        char name[1024];
        getline(name,sizeof (name));
        if (strlen(name) == 0 || name[0] != marker)
        {
          if (eof())
             return string();
          else if (strlen(name) == 0) {
              // throw away last empty line
              getline(name, sizeof(name));
              if (eof())
                return string();
          }
          throw  this->exception(string("Missing ") + marker);
        }
        char *space = strchr(name,' ');
        if (space != NULL)
            space = '\0';
        char *tab = strchr(name,'\t');
        if (tab !=NULL)
            tab = '\0';

    return string(name+1);
  }

  runtime_error exception(const std::string &msg) 
  {
    char buffer[1024];
    sprintf(buffer, "%s in file '%s' at line %ld", msg.c_str(),_path.c_str(), _line);
    return runtime_error(buffer);
  }

  virtual string getName() = 0;
  virtual string getBases() = 0;
  virtual string getQuals() = 0;
};


class fastq_stream : public fast_stream_base
{

public:
  fastq_stream(const string &path) : fast_stream_base(path)
  {
  }

  std::string getName() { return get_name('@'); }

  std::string getBases()
  {
    char bases[MAX_LINE_LENGTH+1];
    char qualname[1024];

    getline(bases,sizeof (bases));
    unsigned long basesSize = strlen(bases);

    if (basesSize == 0 || basesSize > MAX_LINE_LENGTH-1)
        throw  exception("Missing or too many bases");

    getline(qualname,sizeof (qualname));
    if (strlen(qualname) == 0 || qualname[0] != '+')
        throw  exception("Missing '+'");

    return string(bases);
  }

  std::string getQuals()
  {
    char quals[MAX_LINE_LENGTH+1];
    getline(quals,sizeof (quals));
    return string(quals);
  }
};


class fasta_stream : public fast_stream_base
{
  static const char fasta_marker = '>';
  string _qualString;
public:
  fasta_stream(const string &path) : fast_stream_base(path)
  {
    
  }

  virtual std::string getName() {  return get_name(fasta_marker); }

 
  virtual std::string getBases()
  {
    char bases[MAX_LINE_LENGTH+1];

    string baseString;
    while (good())
    {
       getline(bases,sizeof (bases));
       unsigned long basesSize = strlen(bases);

      if ( basesSize > MAX_LINE_LENGTH-1)
         throw  exception("Too many bases");
      
      if (basesSize == 0)
        break;

      if (peek() == fasta_marker)
         break;

      baseString += bases;
    }
  
    _qualString.assign(baseString.length(),255);
    return baseString;
  }
  

  virtual std::string getQuals()
  {
    return _qualString;
  }
};


class fasta_qual_stream : public fasta_stream
{
  fasta_stream _qual_stream;
public:
  fasta_qual_stream(const string &path,const string &qualPath) : fasta_stream(path) , _qual_stream(qualPath)
  {
      
  }

  std::string getName()
  {
     string qualName = _qual_stream.getName();
     string name     =  fasta_stream::getName();
     if (name != qualName)
       throw exception("fasta and quals have different names");
     return name;
  }


  std::string getQuals()
  {
    return _qual_stream.getBases(); // odd!!!
  }
};




void ReadSet::appendFasta(string fastaFilePath,string qualFilePath)
{
    fast_stream_base  *reader;    

    if (qualFilePath.empty())
    {        
       // determine file format
       ifstream ifs(fastaFilePath.c_str());
    
      if (ifs.fail())
        throw  runtime_error("Could not open : " + fastaFilePath) ;
    
      bool isFastq = ifs.peek() == '@';
      ifs.close();

      if (isFastq)
        reader = new fastq_stream(fastaFilePath);
      else
        reader = new fasta_stream(fastaFilePath);
    }
    else
       reader = new fasta_qual_stream(fastaFilePath,qualFilePath);

try {
    while (!reader->eof()) {

         string name = reader->getName();
         if (name.empty())
            break;

        string bases = reader->getBases();
        std::transform(bases.begin(),bases.end(),bases.begin(),::toupper);

        string quals = reader->getQuals();

        if (quals.length() != bases.length())
            throw  reader->exception("Number of bases and quals not equal");

        Read read(name, bases, quals);
        push_back( read );
    }
}

catch(...)
{
 delete reader;
 throw;
}
  delete reader;
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

