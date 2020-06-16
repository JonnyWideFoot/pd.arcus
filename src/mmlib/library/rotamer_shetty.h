#ifndef __ROTLIBCONVERT_SHETTY_H
#define __ROTLIBCONVERT_SHETTY_H

#include "library/rotamerlib.h" // Base class

namespace Library
{
	class PD_API RotLibConvert_Shetty : public RotLibConvertBase
	{
		public: // Very limited publc interface - you can just make one...
			RotLibConvert_Shetty();
		protected:
			/// The core of this class is used to call addRotamer() on the library
			virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib );
	};
}

#endif

