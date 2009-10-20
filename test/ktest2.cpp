/***************************************************************************
 *   Copyright (C) 2009 by Craig Furman   *
 *   cfurman@phoenix   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
 
#include <cstdlib>
#include <cstring>
 
#include "ReadSet.h"

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
}
