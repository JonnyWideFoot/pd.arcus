#ifndef __FILTER_BASE_H
#define __FILTER_BASE_H

#include <string>
#include <vector>

#include "object.h"

class MoleculeBase;

//-------------------------------------------------
//
/// \brief  Abstract base class used to define filters based via an arbritary internal function.
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API FilterBase : public Object
{
public:
	/// Public constructor logic
	FilterBase();
	virtual ~FilterBase();
	virtual FilterBase* clone() const = 0;

	virtual bool passes() = 0; ///< Returns true if filter is valid
	virtual std::string reason() { return "Undefined"; } ///< Returns the reason for the pass or fail

protected:
	virtual void setInternalName();
};


//-------------------------------------------------
//
/// \brief  Abstract base class used to define structural filters based via an arbritary internal function.
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API MolFilterBase : public FilterBase
{
public:
	/// Public constructor logic
	MolFilterBase();
	MolFilterBase( const MoleculeBase& _mol );
	virtual ~MolFilterBase(){}
	virtual MolFilterBase* clone() const = 0;

protected:
	virtual void setInternalName();
	const MoleculeBase* m_Mol; ///< A constant pointer to the underlying molecule. We dont want to change the molecule, only observe its properties.	
};

#ifdef SWIG
%template(ObjectContainer_Filter) ObjectContainer<FilterBase>;
#endif

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API FilterContainer: public FilterBase, public ObjectContainer<FilterBase> 
{
public:
	FilterContainer();
	virtual ~FilterContainer(){}
	virtual FilterContainer* clone() const { return new FilterContainer(*this); }

	virtual bool passes(); ///< Returns true if all child filters pass
	virtual std::string reason(); ///< Returns the reason for each child if they have failed

protected:
	virtual void setInternalName();
};

#endif

