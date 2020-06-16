#ifndef __CONFORMER_BUILDER_H
#define __CONFORMER_BUILDER_H

#include "conformerbuilderbase.h"

namespace Manipulator
{
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
	class ConfBuild_RandomWhole: public ConfBuilderBase_Random
	{
	public:
		ConfBuild_RandomWhole( 
			Library::AngleSet& angleSet, 
			SegmentDef& segdef
			);

		virtual void info() const;
		virtual bool next();
	};


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
	class ConfBuild_RandomSingle: public ConfBuilderBase_Random
	{
	public:
		ConfBuild_RandomSingle( 
			const Library::AngleSet& angleSet, 
			SegmentDef& segdef
			);

		virtual void info() const;
		virtual bool next();
		bool next( size_t _Index );
	};


	//-------------------------------------------------
	//
	/// \brief Conformer builder that will enumerate all the permutations covered by a given descriptor
	/// \details 
	/// Example descriptors: (Examples assume conformer length is 10)
	/// '*' - Exhaustive descriptor, encompases all permutations. (Default)
	/// 'NNNNNNNNNN' or DescriptorType::NativeBin. Encompases only the angles in the same torsional bin as the native state.
	/// 'CCCCCCCCCC' or DescriptorType::NativeClosest. Encompases only the single angle deemed the native state.
	/// '555??AA555' An arbritary conformer: 5,5,5,any,any,AlphaOnly,AlphaOnly,etc...
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class ConfBuild_Enum: public ConfBuilderBase_Enum
	{
	public:
		ConfBuild_Enum( 
			Library::AngleSet& angleSet, 
			SegmentDef& segdef
			);
		virtual void info() const;
	};


	//-------------------------------------------------
	//
	/// \brief Conformer builder that can store conformers from another conformer builder as it is being enumerated
	/// \details  Usage:
	/// 1) Invoke any conformer builder (builder)
	/// 2) Invoke an instance of ConfBuild_BuilderBuffer (store)
	/// 3) Call builder.next() until a desirable conformer is found
	/// 4) Call store.push_back() to cache that conformer
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class ConfBuild_BuilderBuffer: public ConfBuilderBase_Buffer
	{
	public:
		ConfBuild_BuilderBuffer( ConfBuilderBase& _Builder );
		bool push_back(); ///< Append the current ConformerBuilder conformer to the buffer. Returns false if we are over capacity.
		virtual void info() const;
	private:
		void init() const;
		ConfBuilderBase* m_Builder;
	};


	//-------------------------------------------------
	//
	/// \brief  Conformer builder that load a conformer set from a given file
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
	class ConfBuild_FromFileBuffer: public ConfBuilderBase_Buffer
	{
	public:
		ConfBuild_FromFileBuffer( 
			Library::AngleSet& angleSet, 
			SegmentDef& segdef
			);
		void load( const std::string& fileName );
		virtual void info() const;
	};
}
#endif

