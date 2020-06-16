// Example Module
//
// Shows how functionality can be added to mmlib through a module. Use this as a template
// for creating real modules.
//
#ifndef __EXAMPLE_H
#define __EXAMPLE_H

#include "forcefields/forcefield.h"

namespace Physics
{
	void exampleStandAloneFunction();






	//-------------------------------------------------
	//
	/// \brief  This is an example class which derives from a base class defined in PD 
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author You! 
	///
	/// \todo STATE OF DEVELOPMENT?
	///
	/// \bug BUGS?
	///
	class ExampleModuleForcefield: public Physics::ForcefieldBase
	{
	public:
		ExampleModuleForcefield( WorkSpace &newwspace ): 
			ForcefieldBase( newwspace ) 
		{

		}

		virtual ~ExampleModuleForcefield()
		{
		}

		virtual ExampleModuleForcefield *clone() const 
		{
			return new ExampleModuleForcefield(*this);
		}

		virtual int setup(){ return 0; }
		virtual void calcForces(){}

	protected:
	};

}

#endif

