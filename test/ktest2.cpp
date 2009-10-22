// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.7 2009-10-22 07:04:03 regan Exp $
//

#include <iostream>

#include <cstdlib>
#include <cstring>

#include <tr1/unordered_map>

#include "ReadSet.h"
#include "Kmer.h"

using namespace std;

int main(int argc, char *argv[]) {
    
    cerr << "Hello, world! " << endl;

    ReadSet store;

    store.appendFastq(argv[1]);

    cerr << "loaded " << store.getSize() << " Reads" << endl;
    for (int i=0 ; i < store.getSize(); i++)
    {
       cout << store.getRead(i).toFastq();
    }


    Read s;

    cerr << "And even 0 length Reads work!: " << s.getQuals() << s.getFasta() << s.getName() << s.getLength()<< endl;

    KmerSizer::set(44);
    for (int i=0 ; i < store.getSize(); i++)
    {
       KmerArray kmers(store.getRead(i).getTwoBitSequence(),store.getRead(i).getLength());
       for (int j=0; j < kmers.size(); j++)
       {
          cout << "Read " << i << " kmer "<< j << ' '
               << kmers[j].toFasta() << endl;
       }
       cout << "----------------" << endl;
    }
}


//
// $Log: ktest2.cpp,v $
// Revision 1.7  2009-10-22 07:04:03  regan
// added a few unit tests
// minor refactor
//
// Revision 1.6  2009-10-22 01:39:46  cfurman
// bug fix in kmer.h
//
// Revision 1.5  2009-10-21 06:51:37  regan
// bug fixes
// build lookup tables for twobitsequence
//
// Revision 1.4  2009-10-21 00:02:02  cfurman
// working on kmers....
//
// Revision 1.3  2009-10-20 20:56:29  cfurman
// Got it to compile!
//
// Revision 1.2  2009-10-20 17:25:53  regan
// added CVS tags
//
//
