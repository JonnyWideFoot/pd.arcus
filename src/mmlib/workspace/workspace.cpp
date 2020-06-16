#include "global.h"

#include "forcefields/forcefield.h"
#include "system/workspacecreator.h"
#include "system/system.h"
#include "fileio/tra.h"
#include "exception.h"

// Used workspace operators
#include "neighbourlist.h"
#include "bondorder.h"
#include "rotbond.h"
#include "space.h"

#include "workspace.h"

// namespace includes
using namespace Maths;
using namespace Physics;
using namespace Library;

/// Generic constructor - workspace is created using a workspace
/// creator passed by argument
WorkSpace::WorkSpace(const WorkspaceCreatorBase &wc) : 
	MoleculeBase(wc.ffps()),
	outtra( *this )
{
	CommonInit();
	wc.create(*this);
	parentspec = wc.system;
	//reinitAll(); // Should already be done by 'create()' above
}

/// This is a 'shortcut' version which directly takes a sysspec
/// and automatically creates the WorkSpace
WorkSpace::WorkSpace(System &sysspec) : 
	MoleculeBase(sysspec.ffps()),
	outtra( *this )
{
	CommonInit();
	WSfull(sysspec).create(*this);
	parentspec = &sysspec;
	//reinitAll(); // Should already be done by 'create()' above
}

WorkSpace::WorkSpace(const WorkSpace &copywspace):
	outtra( *this )
{
	CommonInit();
	(*this) = copywspace;
	reinitAll();
}

void WorkSpace::CommonInit()
{
	Step = 0;
	default_rotbond   = new RotBond_Dummy(); 
	default_nlist     = new NeighbourList();
	default_bondorder = new BondOrder();
	default_boundary  = new InfiniteSpace();

	ptr_rotbondlist = default_rotbond;
	ptr_nlist       = default_nlist;
	ptr_bondorder   = default_bondorder;
	ptr_boundary    = default_boundary;

	m_CheckSum = 0;

	//outtra.linkWSpace(*this); // link the outtra to its parent wspace to trigger ene calcs
}

void WorkSpace::Allocate(int _natoms)
{
	cur = SnapShot(_natoms);
	old = SnapShot(_natoms);
	m_CheckSum ++;
}

void WorkSpace::reinitAll()
{
	ptr_bondorder->reinit_base(this);
	ptr_nlist->reinit_base(this);
	ptr_rotbondlist->reinit_base(this);
	m_CheckSum ++;
	parseGroups();
}


void WorkSpace::parseGroups()
{
	int igroup=0;
	size_t incTo = atom.size() - 1;
	for(size_t i = 0; i < incTo; i++)
	{
		atom[i].igroup = igroup;
		
		if( mol[atom[i].imol].natoms() <= 3)
		{
			if(atom[i].imol != atom[i+1].imol) 
			{
				igroup++; // only advance on molecule boundaries
			}
			continue;
		}
		
		igroup++;
	}
	atom[incTo].igroup = igroup;

	group_index.clear();
	group_index.push_back(0);
	for(size_t i = 1; i < atom.size(); i++)
	{
		if( atom[i].igroup != atom[i-1].igroup )
		{
			group_index.push_back(i);
		}
	}
}

WorkSpace::~WorkSpace() 
{
	delete default_rotbond;
	delete default_nlist;	
	delete default_bondorder;
	delete default_boundary;
}

char WorkSpace::getChainID( size_t iAtom ) const
{
	if( iAtom >= nAtoms())
	{
		THROW(OutOfRangeException,"getChainID() request out of bounds");
	}
	return mol[ atom[iAtom].imol ].chainID;
}

void WorkSpace::setRotatableBondList(RotBond *new_rotbondlist)
{
	if(new_rotbondlist == NULL) throw(ArgumentException("Cannot set RotatableBondList to NULL pointer\n"));
	ptr_rotbondlist = new_rotbondlist;
	ptr_rotbondlist->reinit_base(this);
	m_CheckSum++;
}


// setting operators
void WorkSpace::setNeighbourList(NeighbourListBase *new_nlist)
{
	if(new_nlist == NULL) throw(ArgumentException("Cannot set DeprecatedNeighbourList to NULL pointer\n"));
	ptr_nlist = new_nlist;
	ptr_nlist->reinit_base(this);
	m_CheckSum ++;
}

// setting operators
void WorkSpace::setSpace(Space *new_boundary)
{
	if(new_boundary == NULL) throw(ArgumentException("Cannot set Space to NULL pointer\n"));
	ptr_boundary = new_boundary;
	m_CheckSum ++;
}

// returns a copy of cur + the current energy (i.e. ehatever happens to be in ene.epot)
SnapShot WorkSpace::save() const 
{ 
	// yuk - that means copying it twice!(i think) but
	// cant think of other way of doing it without breaking 
	// const requirement of this function
	SnapShot psp(cur); 
	psp.epot = ene.epot;
	return psp;
}

void WorkSpace::load(const SnapShot &psp)
{ 
	if(nAtoms() != psp.natoms){	
		throw(ArgumentException(std::string("ERROR: SnapShot and WorkSpace are	incompatible\n - number of atoms does not match:\nWorkspace: ")
			+ int2str(nAtoms()) + std::string("  SnapShot:") + int2str(psp.natoms)));
	}
	cur = psp;
}


void WorkSpace::load_forced(const SnapShot &psp)
{ 
	if(nAtoms() != psp.natoms){	
		printf("WARNING: SnapShot and WorkSpace are incompatible: ");
		printf("number of atoms does not match:\n");
		printf("  WorkSpace: %d  SnapShot: %d",nAtoms(),psp.natoms);
		printf("  Loading forced: %d atoms loaded\n",min((int)nAtoms(),(int)psp.natoms));
	}
	cur.setToPartial(psp);
}

int WorkSpace::updateSystemPositions(){
	// Should check if sysspec has not changed
	// check if all atoms have backlinks
	for(size_t i = 0; i < nAtoms() ; i++)
	{
		if((isysmol[i] < 0)||(isysatom[i] < 0))
		{
		  THROW(ProcedureException,"Cannot update SystemSpec because the binding is broken.");
		}
	}

	// update the parent sysspec
	for(size_t i=0;i<nAtoms();i++)
	{
		parentspec->getMolecule(isysmol[i]).getAtom(isysatom[i]).pos() = cur.atom[i].p;
	}
	return 0;
}

void WorkSpace::addTraWithOwnership(IO::OutputTrajectory *newtra)
{
	outtra.addWithOwnership(newtra);
}

void WorkSpace::addTra(IO::OutputTrajectory &newtra)
{
	outtra.add(newtra);
}

IO::OutTra_BTF& WorkSpace::addStdTra(const std::string &_filestem)
{
	IO::OutTra_BTF* tra = new IO::OutTra_BTF(_filestem,*this);
	outtra.addWithOwnership(tra);
	return *tra;
}

IO::OutTra_BTF& WorkSpace::addExtdTra(const std::string &_filestem, bool _MoreModes )
{
	// Create BristolTrajectoryFormat With Extra Extensions ...
	IO::OutTra_BTF* tra = new IO::OutTra_BTF(_filestem,*this);
	
	IO::BTF_Block_Vector* vect = new IO::BTF_Block_Vector(); // due to a bug in dave, BTF_Block_Vector must be added 1st...
	tra->addOwnedBlock(vect);
	IO::BTF_Block_Comment* cmnt = new IO::BTF_Block_Comment();
	tra->addOwnedBlock(cmnt);
	
	if( _MoreModes )
	{
		tra->create(IO::AllIncludes);
	}
	else
	{
		tra->create(IO::DefaultIncludes);
	}

	outtra.addWithOwnership(tra);

	return *tra;
}

double WorkSpace::calcCRMS_AllAtom( bool useOnlyKnownAtom ) const
{
	double cRMS = MoleculeBase::calcCRMS_AllAtom( useOnlyKnownAtom );
	ene.cRMS = cRMS;
	return cRMS;
}

double WorkSpace::calcCRMS_HeavyAtom( bool useOnlyKnownAtom ) const
{
	double cRMS = MoleculeBase::calcCRMS_HeavyAtom( useOnlyKnownAtom );
	ene.cRMS = cRMS;
	return cRMS;
}

double WorkSpace::calcCRMS_CA( bool useOnlyKnownAtom ) const
{
	double cRMS = MoleculeBase::calcCRMS_CA( useOnlyKnownAtom );
	ene.cRMS = cRMS;
	return cRMS;
}

void WorkSpace::save( IO::OutputFile &_output )
{
	_output.save( *this );
}


void WorkSpace::zeroForces()
{
	ene.zero();
	for(size_t i = 0; i < nAtoms(); i++)
	{
		cur.atom[i].f.setTo(0, 0, 0);
		getAtom(i).epot = 0;
	}

}

void WorkSpace::scaleVelocities(double factor)
{
	for(size_t i=0; i < nAtoms(); i++)
	{
		cur.atom[i].v.mul( factor );
	}
}

double WorkSpace::getVolume() const
{
	ClosedSpace *cboundary = dynamic_cast<ClosedSpace*> (ptr_boundary);
	if(cboundary == NULL){
			return -1.0; // no volume for open spaces
	}
	return cboundary->volume();        // return volume otherwise 
}

double WorkSpace::getDensity() const
{
	double V = getVolume();
	return getTotalMass()*1000/(V*1E-24);  // return density in g/ml
}

void   WorkSpace::scaleSystemMolecular(double factor){
	for(size_t imol=0; imol < mol.size(); imol++)
	{
		dvector cog = getCentreOfGeometry(imol);
		dvector newcog = cog;
		newcog.mul(factor);
		newcog.sub(cog);
		moveMolecule(imol,newcog);
	}
	boundary().scale(factor);
}

void   WorkSpace::scaleSystem(double factor){
	for(size_t i=0; i < nAtoms(); i++)
	{
		cur.atom[i].p.mul(factor);
	}
	boundary().scale(factor);
}


// Move stray molecules back into simulation box
void WorkSpace::cleanSpace()
{
	for(size_t imol=0; imol < nMolecules(); imol++){
		dvector cog(getCentreOfGeometry(imol));
		dvector cog_inbox(cog);
		ptr_boundary->moveIntoBox(cog_inbox);
		cog_inbox.sub(cog);
		moveMolecule(imol,cog_inbox);
	}
}

// -------------------------------------------------------------------------------
// Class: WorkSpace
// Function: printForces(); <debug>
// -------------------------------------------------------------------------------
// Parameters: none
//
// prints forces for debugging purposes

void WorkSpace::printForces(){
	int i;
	printf(" nr atname restype resnum atommass force speed \n");

	for(i = 0; i < nAtoms(); i++) 
	{
		printf(" %4i %-4s%3s %3i % 6.3e % 6.3e % 6.3e % 6.3e % 6.3e % 6.3e \n",
			i,
			getAtom(i).pdbname.c_str(),
			getAtom(i).parentl3name.c_str(),
			getAtom(i).ir, getAtom(i).mass, cur.atom[i].f.x, cur.atom[i].f.y, cur.atom[i].f.z, cur.atom[i].f.mag(), cur.atom[i].v.mag()
			);
	}
}

size_t WorkSpace::memuse(int level)
{
	level--;
	size_t s_self=sizeof(MoleculeBase);
	size_t s_atom = 0;
	for(size_t i = 0; i < nAtoms(); i++) s_atom += getAtom(i).memuse(level);
	size_t s_res = 0;
	for(size_t i = 0; i < nResidues(); i++) s_res += getRes(i).memuse(level);
	size_t s_mol = 0;
	for(size_t i = 0; i < nMolecules(); i++) s_mol += getMol(i).memuse(level);
	size_t s_isysxxx = 2 * sizeof(int)*nAtoms();
	size_t s_cur = cur.memuse(level);
	size_t s_old = old.memuse(level);
	size_t s_nlist = ptr_nlist->memuse(level);   
	size_t s_rotlist = ptr_rotbondlist->memuse(level);
	size_t s_bondorder = ptr_bondorder->memuse(level);

	if(level >= 0)
	{
		printf("WorkSpace (total)  %d : \n",(int)(s_self+ s_atom+ s_res+ s_mol+
							s_isysxxx+s_cur+s_old+s_nlist+s_rotlist+s_bondorder));
		printf("  self: %d\n"		 ,(int)s_self);
		printf("  atom: %d\n"		 ,(int)s_atom);
		printf("  res: %d\n"		 ,(int)s_res);
		printf("  mol: %d\n"		 ,(int)s_mol);			
		printf("  isysxx: %d\n"  ,(int)s_isysxxx); 
		printf("  cur: %d\n"		 ,(int)s_cur);
		printf("  old: %d\n"		 ,(int)s_old);
		printf("  ptr_nlist: %d\n"   ,(int)s_nlist);   
		printf("  rotlist: %d\n" ,(int)s_rotlist);
		printf("  ptr_bondorder: %d\n",(int)s_bondorder);
	}
	return s_self+ s_atom+ s_res+ s_mol+
				 s_isysxxx+s_cur+s_old+
				 s_nlist+s_rotlist+s_bondorder;
}

