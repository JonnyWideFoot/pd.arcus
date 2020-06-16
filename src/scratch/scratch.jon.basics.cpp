#include "global.h"

#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "forcefields/nonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"
#include "forcefields/ffsoftvdw.h"
#include "workspace/workspace.h"

#include "scratch.jon.basics.h"

using namespace Physics;

Forcefield createffs( WorkSpace& wspace, bool useBreakableFF, bool summary )
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	SoftVDWForcefield *sff = new SoftVDWForcefield(wspace);
	ff.addWithOwnership(sff);

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createffts( WorkSpace& wspace, bool useBreakableFF, bool summary)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		bonds->DoBonds = false;
		bonds->DoAngles = false;
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		bonds->DoBonds = false;
		bonds->DoAngles = false;
		ff.addWithOwnership( bonds ) ;
	}

	SoftVDWForcefield *sff = new SoftVDWForcefield(wspace);
	ff.addWithOwnership(sff);

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createffVac(WorkSpace& wspace, bool useBreakableFF, bool summary)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	FF_NonBonded* nb = new FF_NonBonded(wspace);
	nb->Cutoff = 12.0;
	nb->InnerCutoff = 6.0;
	ff.addWithOwnership( nb );

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createff(WorkSpace& wspace, bool useBreakableFF, double dielec, bool summary )
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	GB_Still* gbsa = new GB_Still( wspace ); // used to take nb
	gbsa->FastMode = true;
	gbsa->DielectricSolute = dielec;
	ff.addWithOwnership( gbsa );

	SASA_LCPO* sasa = new SASA_LCPO(wspace);
	sasa->GlobalASP = 0.009;
	ff.addWithOwnership( sasa );

	if( summary ) ff.printEnergySummary();

	return ff;
}

double getMeRMS( const std::vector<Maths::dvector>& native, const std::vector<Maths::dvector>& conformer )
{
	double sumSqr = 0.0;
	for( size_t i = 0; i < native.size(); i++ )
	{
		sumSqr += native[i].sqrdist( conformer[i] );
	}
	return std::sqrt( sumSqr / (double)native.size() );
}

