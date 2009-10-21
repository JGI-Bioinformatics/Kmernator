// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.4 2009-10-21 00:02:02 cfurman Exp $
 
#include <iostream>

#include <cstdlib>
#include <cstring>

#include <tr1/unordered_map>

#include "../src/ReadSet.h"
#include "../src/Kmer.h"


using namespace std;




int main(int argc, char *argv[]) {


    KmerReadSetStatsMap kmap;

    Kmer my_mer("ACTGGTTT");

    kmap[my_mer] = KmerReadSetStats(10.1,sdasdasdasdas);
    
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
}

//
// $Log: ktest2.cpp,v $
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
