#ifndef __ARCUS_STITCH_H
#define __ARCUS_STITCH_H

#include "arcus_base.h"

namespace Protocol
{
	class PD_API ArcusStitch : public ArcusBase
	{
	public:
		ArcusStitch( Physics::Forcefield& _ffs, Physics::Forcefield& _ff, Library::AngleSet& _as );
		virtual ArcusStitch* clone() const 
		{ 
			return new ArcusStitch(*this); 
		}
		virtual int runcore();
	};
}

#endif

