/*
 * Log.cpp
 *
 *  Created on: Oct 13, 2010
 *      Author: regan
 */

#include "Log.h"

Logger Log::verboseOstream = Logger( std::cerr, "INFO", 1 );
Logger Log::debugOstream = Logger( std::cerr, "DEBUG", 0 );
Logger Log::warningOstream = Logger( std::cerr, "WARNING", 1 );
Logger Log::errorOstream = Logger( std::cerr, "ERROR", 1 );
