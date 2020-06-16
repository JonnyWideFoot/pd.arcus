#include "global.h"
#include "manipulators/movebase.h"
#include "workspace/workspace.h"

namespace Manipulator
{
	MoveBase::MoveBase(WorkSpace& _wspace) 
		: wspace(&_wspace)
	{
	}

	MoveBase::~MoveBase()
	{
	}

	// Applies the move 'rounds' times for testing purposes
	void MoveBase::test(int rounds)
	{
		for(int i = 0; i < rounds; i++) 
		{
			apply();
			wspace->outtra.append();
		}
	}

	MoveSet::MoveSet(WorkSpace& newwspace)
		: MoveBase(newwspace) 
	{
		name = "MoveSet";
	}
		
	MoveSet::~MoveSet()
	{
	}

	MoveSet* MoveSet::clone() const 
	{ 
		return new MoveSet(*this); 
	}

	int MoveSet::apply()
	{
		for(size_t i = 0; i < size(); i++)
		{
			element(i).apply();
		}
		return 0;
	}
}

