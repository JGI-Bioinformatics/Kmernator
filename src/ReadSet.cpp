// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.cpp,v 1.4 2009-10-23 07:06:59 regan Exp $
//

#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>

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


class ReadSet::body
{
public:

    std::vector<Read>  seqs;
};

ReadSet::ReadSet():
_my(*new(body))
{
}

ReadSet::~ReadSet()
{
  delete &_my;
}


string fileErrorMsg(string msg, string name, unsigned long lineCount)
{

  char buffer[1024];

  sprintf(buffer, "%s in %s at %ld", msg.c_str(),name.c_str(), lineCount);
  return string(buffer);
}

void ReadSet::appendFastq(string fastqFilePath)
{
 try {
    ifstream ifs(fastqFilePath.c_str());

    if (ifs.fail()) {
       throw  std::invalid_argument("Could not open : " + fastqFilePath) ;
    }
    unsigned long fileSize = _fileSize(ifs);
    unsigned long lineCount = 0;

    char name[1024];
    char bases[MAX_LINE_LENGTH+1];
    char qualname[1024];
    char quals[MAX_LINE_LENGTH+1];
    while (!ifs.eof()) {

        ifs.getline(name,sizeof (name));
        if (strlen(name) == 0 || name[0] != '@')
        {
          if (ifs.eof())
             break;
          else if (strlen(name) == 0) {
          	 // throw away last empty line
          	 ifs.getline(name, sizeof(name));
          	 if (ifs.eof())
          	   break;
          }   
          throw   std::invalid_argument(fileErrorMsg("Missing @", fastqFilePath, lineCount));
        }
        char *space = strchr(name+1,' ');
        if (space != NULL)
           space = '\0';
        char *tab = strchr(name+1,'\t');
        if (tab !=NULL)
            tab = '\0';

        ifs.getline(bases,sizeof (bases));
        unsigned long basesSize = strlen(bases);
        if (basesSize == 0 || basesSize > MAX_LINE_LENGTH-1)
              throw  std::invalid_argument(fileErrorMsg("Missing bases or too many", fastqFilePath, lineCount));

        ifs.getline(qualname,sizeof (qualname));
        if (strlen(qualname) == 0 || qualname[0] != '+')
              throw  std::invalid_argument(fileErrorMsg("Missing '+'", fastqFilePath, lineCount));

        ifs.getline(quals,sizeof (quals));
        if (strlen(quals) != basesSize)
              throw  std::invalid_argument(fileErrorMsg("Wrong number of quals", fastqFilePath, lineCount));

         if (lineCount == 0)  // Estimate set size
         {
             //_my.seqs.reserve(1+ _my.seqs.size() + fileSize/((unsigned long)ifs.tellg() - 10UL));
             cerr << "Capacity : " << _my.seqs.capacity()<< endl;
        }

         lineCount += 4;

        _my.seqs.push_back( Read(name+1 , bases, quals));


    }
    ifs.close();
  } catch (...) {
    throw;
  }
}

ReadSetSizeType ReadSet::getSize()
{
  return _my.seqs.size();
}

Read &ReadSet::getRead(ReadSetSizeType index)
{
  return _my.seqs[index];
}

//
// $Log: ReadSet.cpp,v $
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

