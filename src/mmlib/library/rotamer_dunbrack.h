#ifndef __ROTLIBCONVERT_DUNBRACK_H
#define __ROTLIBCONVERT_DUNBRACK_H

#include "library/rotamerlib.h" // Base class

namespace Library
{
	class PD_API RotLibConvert_Dunbrack_BBInd : public RotLibConvertBase
	{
		public: // Very limited publc interface - you can just make one...
			RotLibConvert_Dunbrack_BBInd();
		protected:
			/// The core of this class is used to call addRotamer() on the library
			virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib );
		private:
			/// Read a single torsional definition line, derived from the import file, into the library
			void readDefinition( StringBuilder& sb, RotamerLibrary& _RotLib ) const;
	};	

	class PD_API RotLibConvert_Dunbrack_BBDep : public RotLibConvertBase
	{
		public: // Very limited publc interface - you can just make one...
			RotLibConvert_Dunbrack_BBDep();

			bool switchImportMode;

		protected:
			/// The core of this class is used to call addRotamer() on the library
			virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib );
		private:
			/// Read a single torsional definition line, derived from the import file, into the library
			void readDefinition( StringBuilder& sb, RotamerLibrary& _RotLib ) const;
	};	
}

#endif

