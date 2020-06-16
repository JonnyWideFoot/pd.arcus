#include "global.h"
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "forcefields/forcefield.h"
#include "workspace/snapshot.h"
#include "protocols/md.h"
#include "energy.h"

using namespace Physics;
using namespace Maths;

namespace Protocol
{
	int Energy::runcore()
	{
		if(ensureFFSetup()!=0) return -1;
		refreshNeighborList();
		ff->calcEnergies();
		return 1;
	}
/*
	int ReadPickles::runcore(){

		size_t is=0;
		
		FILE *file;
		file = fopen(filename.c_str(),"r");
		if(file==NULL){ printf("Cannot find file \n"); return -1; }

		printf("Loading all pickles and saving in trajectory \n");

		while(!feof(file)){
			if((is>0)&&((is%100)==0)) printf("%d...",is);
			SnapShot psp;
			psp.readMIME(file);
			getWSpace().load(psp);
			getWSpace().outtra.append();
			is++;
		}
		printf("Read %d snapshots \n",is);
		fclose(file);
		return 0;
	}
*/

} // namespace 'Protocol'


