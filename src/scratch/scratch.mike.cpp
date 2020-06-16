#include "global.h"

#include "../pd/accessory.h"
#include "tools/statclock.h"

#include "system/fundamentals.h"

#include "system/genpolymer.h"
#include "workspace/workspace.h"

#include "sequence/sequence.h"
#include "sequence/alignment.h"

#include "forcefields/ffbonded.h"
#include "forcefields/ffnonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"

#include "fileio/tra.h"
#include "fileio/pdb.h"
#include "tools/io.h"

#include "protocols/energy.h"
#include "protocols/minimise.h"
#include "protocols/torsionalminimisation.h"
#include "protocols/md.h"

#include "workspace/space.h"

using namespace std;
using namespace Physics;
using namespace Protocol;
using namespace Tra;
using namespace PDB;
using namespace Sequence;

#include <omp.h>

int main()
{
  omp_set_dynamic(0);
  omp_set_num_threads(4);

	double fpub = 1;
		int its=0;
	#pragma omp parallel
	{	 
		int i;

		int j;
		#pragma omp for //schedule(dynamic, 10)
		for(i=0;i<100;i++){
			//fpub = 0;
			//fpub += (double)i;
			//printf("[%2d]:   i = %d  fpub: %f \n", omp_get_thread_num(), i, fpub);
			//its[omp_get_thread_num()] += 1; // remember howmany were executed
			for(j=0;j<100;j++){
				its++;
			}
		}
	}

	printf("result:  %d \n", its);


}



