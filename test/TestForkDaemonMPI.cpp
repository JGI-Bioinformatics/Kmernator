/*
 * TestForkDaemonMPI.cpp
 *
 *  Created on: Oct 29, 2013
 *      Author: regan
 */


#include "Log.h"
#include "Options.h"
#include "Utils.h"
#include "MPIUtils.h"


class _MPITestForkDaemonOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		GeneralOptions::getOptions().getDebug() = 2;
		GeneralOptions::getOptions().getVerbose() = 1;
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		return ret;
	}
};

typedef OptionsBaseTemplate< _MPITestForkDaemonOptions > MPITestForkDaemonOptions;

int main(int argc, char **argv) {

	{
		ForkDaemon::initialize();
		ForkDaemon &forkDaemon = *ForkDaemon::getInstance();

		{
			ScopedMPIComm<MPITestForkDaemonOptions> world(argc, argv);

			pid_t pid1 = forkDaemon.startNewChild("/bin/uname -a");
			pid_t pid2 = forkDaemon.startNewChild("/bin/sleep 2");
			forkDaemon.startNewChild("/bin/sleep 7");

			int status = 0;
			status = forkDaemon.waitChild(pid1);
			LOG_DEBUG(1, "uname returned: " << status << " (" << pid1 << ")");
			status = forkDaemon.waitChild(pid2);
			LOG_DEBUG(1, "sleep 2 returned: " << status << " (" << pid2 << ")");

		}

		LOG_DEBUG(1, "MPI now out of scope");

		pid_t pid1 = forkDaemon.startNewChild("/bin/uname -a");
		pid_t pid2 = forkDaemon.startNewChild("/bin/sleep 1");
		forkDaemon.startNewChild("/bin/sleep 3");

		int status = 0;
		status = forkDaemon.waitChild(pid1);
		LOG_DEBUG(1, "uname returned: " << status << " (" << pid1 << ")");
		status = forkDaemon.waitChild(pid2);
		LOG_DEBUG(1, "sleep 1 returned: " << status << " (" << pid2 << ")");

	}

	ForkDaemon::finalize();
	LOG_DEBUG(1, "forkDaemon now out of scope");
}
