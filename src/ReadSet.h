// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.2 2009-10-20 17:25:50 regan Exp $
//

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

//
// $Log: ReadSet.h,v $
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//
