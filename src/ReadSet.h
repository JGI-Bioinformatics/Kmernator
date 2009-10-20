#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>
#include "Sequence.h"


typedef unsigned int ReadSetSizeType;

class ReadSet {
    class body;
    body &_my;
public:
    ReadSet();
    ~ReadSet();

    void appendFastq(std::string fastqFilePath);

    ReadSetSizeType getSize();
     
    Read &getRead(ReadSetSizeType index);

};
#endif
