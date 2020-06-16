#ifndef __LOOP_CA_DIST_FILTER_H
#define __LOOP_CA_DIST_FILTER_H

#include "tools/cloneholder.h"
#include "filters/filterbase.h"

class PickResidueRange;
class ProximityGrid;

class PD_API LoopCADistFilter : public FilterBase
{
public:
	LoopCADistFilter();
	virtual ~LoopCADistFilter(){};
	virtual LoopCADistFilter* clone() const { return new LoopCADistFilter(*this); }

	void setTo( const ProximityGrid& _StaticCAGrid, const MoleculeBase& mol, const PickResidueRange& picker, size_t _NOffset, size_t _COffset  );
	virtual bool passes();
	virtual std::string reason();

	double SqrCAMinReach;

protected:
	const MoleculeBase* m_Mol;
	const ProximityGrid* m_StaticCAGrid;
	CloneHolder< PickResidueRange > m_Picker;
	size_t m_NOffset; 
	size_t m_COffset;
};

#endif

