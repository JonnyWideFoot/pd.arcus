#include "global.h"
#include <sstream>
#include "library/angleset.h"
#include "library/backbonetorsions.h"
#include "conformerbuilder.h"

// namespace includes
using namespace Library;

namespace Manipulator
{
	// --------------------------------------------------------
	// ConfBuild_RandomWhole
	// --------------------------------------------------------
	ConfBuild_RandomWhole::ConfBuild_RandomWhole( 
		Library::AngleSet& angleSet, 
		SegmentDef& segdef
		) : 
		// Initialise base classes
		ConfBuilderBase_Random(angleSet, segdef),
		SegmentDefUser(segdef)
	{
		reset();
	}

	void ConfBuild_RandomWhole::info() const
	{
		std::cout << "Random Whole Conformer Change Builder" << std::endl;
		ConfBuilderBase::info(); // explicitly call the base class info() function

		printf("\t");

		for( size_t i = 0; i < segLength(); i++ )
		{
			printf("%d", m_CurrentConformer[i] );
		}
		printf(": Is the current random conformer\n");
		if( getRandCount() == 0 )
		{
			printf("\tMaximim random changes is defined as infinity!\n");
		}
		else
		{
			printf("\tA maximum of %ld random whole-conformer changes will be performed in total.\n", getRandCount());
		}
	}

	bool ConfBuild_RandomWhole::next()
	{
		if( getRandCount() != 0 && m_CurrentIndex >= getRandCount() ) return false; ///< we have reached the end
		m_CurrentIndex++; ///< increment
		RandomiseWholeConformer();
		return true;
	}

	// --------------------------------------------------------
	// ConfBuild_RandomSingle
	// --------------------------------------------------------

	ConfBuild_RandomSingle::ConfBuild_RandomSingle( 
		const Library::AngleSet& angleSet, 
		SegmentDef& segdef
		): 
		// Initialise base classes
		ConfBuilderBase_Random(angleSet, segdef),
		SegmentDefUser(segdef)
	{
		reset();
	}

	void ConfBuild_RandomSingle::info() const
	{
		std::cout << "Random Single Residue Change Builder" << std::endl;
		ConfBuilderBase::info(); // explicitly call the base class info() function
		ConfBuilderBase_Descriptor::infoDescriptor();
		printf("\t");

		for( size_t i = 0; i < segLength(); i++ )
		{
			printf("%d", m_CurrentConformer[i] );
		}

		printf(": Is the current random conformer\n");
		if( getRandCount() == 0 )
		{
			printf("\tMaximim random changes is defined as infinity!\n");
		}
		else
		{
			printf("\tA maximum of %ld random single-residue conformer changes will be performed in total.\n", getRandCount());
		}
	}

	bool ConfBuild_RandomSingle::next()
	{		
		if( getRandCount() != 0 && m_CurrentIndex >= getRandCount() ) 
			return false; // we have reached the end
		m_CurrentIndex++; // increment
		RandomiseSingleConformerPos();
		return true;
	}

	bool ConfBuild_RandomSingle::next( size_t _Index )
	{
		if( getRandCount() != 0 && m_CurrentIndex >= getRandCount() ) 
			return false; // we have reached the end
		m_CurrentIndex++; // increment
		RandomiseSingleConformerPos( _Index );
		return true;
	}

	// --------------------------------------------------------
	// ConfBuild_Enum
	// --------------------------------------------------------

	ConfBuild_Enum::ConfBuild_Enum(
		Library::AngleSet& angleSet, 
		SegmentDef& segdef
		):
		// Initialise base classes
		ConfBuilderBase_Enum(angleSet, segdef),
		SegmentDefUser(segdef)
	{
		addEntireAnglesetLibrary();
		reset();
	}

	void ConfBuild_Enum::info() const
	{
		std::cout << "Enumerating Conformer Builder" << std::endl;
		ConfBuilderBase::info(); // explicitly call the base class info() function
		ConfBuilderBase_Descriptor::infoDescriptor();
	}

	// --------------------------------------------------------
	// ConfBuild_FromFileBuffer
	// --------------------------------------------------------

	ConfBuild_FromFileBuffer::ConfBuild_FromFileBuffer( Library::AngleSet& angleSet, SegmentDef& segdef ):
		// Initialise base classes
		ConfBuilderBase_Buffer(angleSet, segdef),
		SegmentDefUser(segdef)
	{
		// angleset
		// NOTE: The conformer file should also contain the filename of the .angleset file to which it refers!
		// and that should be checked here!
		addEntireAnglesetLibrary(); ///< import all the angles from the anglset.
	}

	void ConfBuild_FromFileBuffer::info() const
	{
		std::cout << "ConformerBuffer ('FromFile' submode)...." << std::endl;
		ConfBuilderBase::info();
	}

	void ConfBuild_FromFileBuffer::load( const std::string& fileName )
	{
		printf("ConformerBuilder file load completed beginning:" );

		// setup the file we are extracting from and create our buffer
		const int BUFFER_LENGTH = 512;
		char buffer[BUFFER_LENGTH];

		// Proxies
		size_t loopLength = getSegDef().getNRes();
		std::string seq = getSegDef().getSequence().printToStringSingle();

		FILE* confFile = fopen(fileName.c_str(),"r");
		if(confFile == NULL)
		{
			std::stringstream stream;
			stream << "Could not find conformer file: '" << fileName << "'";
			THROW(ArgumentException,stream.str());
		}

		// Do sanity checks:
		// Our loop sequence should matche the sequence header in the file we are importing. If not then file is invalid.
		// The sequence length also dictates the length that each of the imported conformers needs to be.
		fgets(&buffer[0], BUFFER_LENGTH, confFile);
		truncLF(&buffer[0]); ///< remove the windows line feeds and replace with a linux-style '/0' string termination
		if( 0 != strcmp(&buffer[0],seq.c_str()) )
		{
			std::stringstream stream;
			stream << "The seqenence imported from the conformer file: '" << fileName << "' does not match that of the loop!";
			THROW(ArgumentException,stream.str());
		}

		// the length of the conformers should be that of the sequence if all is well
		Conformer confBuffer( loopLength ); ///< make our temporary collection buffer.
		size_t buffer_overflow = 0;
		do
		{
			fgets(&buffer[0], BUFFER_LENGTH, confFile);
			truncLF(&buffer[0]); ///< again, remove the windows line feeds
			if( loopLength != strlen(&buffer[0]) )
			{
				std::stringstream stream;
				stream << "Conformer length does not match that of the sequence! ConformerFile: '"
					<< fileName << "', ConformerLine: '" << size() << "'";
				THROW(ProcedureException,stream.str());
			}

			if( (capacity() <= 0) || (size() == capacity()) )
			{
				buffer_overflow++;
			}
			else
			{
				for( size_t i = 0; i < loopLength; i++ )
				{
					confBuffer[i] = (conformer_type)(buffer[i] - 48); ///< convert the ASCII character to a small positive number and store in the conformer_type
				}

				if( !push_back(confBuffer))
				{
					std::stringstream stream;
					stream << "Invalid conformer found! ConformerFile: '"
						<< fileName << "', ConformerLine: '" << m_Conformers.size() << "'";
					THROW(ProcedureException,stream.str());
				}
			}
		}
		while(!feof(confFile)); ///< keep reading the conformer strings while they are present in the file.

		fclose(confFile); ///< close our file handle, we are done with it.

		rotateFromBuffer(0); ///< apply the 1st conformation (important)
		std::cout << "ConformerBuffer now contains a total of '" << m_Conformers.size() << "' conformers" << std::endl;

		if( buffer_overflow != 0 )
		{
			std::cout << "WARNING: File contents exceeded buffer capacity. The final '" <<  buffer_overflow << "' conformers were completely ignored!" << std::endl;
			std::cout << "ConformerBuffer file load completed with warnings!" << std::endl;
		}
		else
		{
			std::cout << "ConformerBuffer file load completed successfully!" << std::endl;
		}
	}


	// --------------------------------------------------------
	// ConfBuild_BuilderBuffer
	// --------------------------------------------------------

	ConfBuild_BuilderBuffer::ConfBuild_BuilderBuffer( ConfBuilderBase& _Builder ):
		// Initialise base classes
		ConfBuilderBase_Buffer(_Builder.getAngleSet(), _Builder.getSegDef() ),
		SegmentDefUser(_Builder.getSegDef()),
		// Internal members
		m_Builder(&_Builder)
	{
		THROW(NotImplementedException,"'ConfBuild_BuilderBuffer' is totally untested, you need to test it!");

		// Copy the whole internal anglset
		cloneBackboneTorsionSubSet( _Builder ); // import all the angles from the anglset.

		// **** This line looks dubious ****... surely we need to get the current builder state 1st ...
		rotateFromBuffer(0); // apply the 1st conformation (important)
	}

	void ConfBuild_BuilderBuffer::info() const
	{
		std::cout << "ConformerBuffer ('FromBuilder' submode)..." << std::endl;
		ConfBuilderBase::info();
		printf("Buffer contains a total of %ld conformers\n", m_Conformers.size());
	}

	bool ConfBuild_BuilderBuffer::push_back()// add the current ConformerBuilder conformer to the buffer
	{
		if( capacity() > 0 && size() == capacity() ) // m_InternalCapacity can be '0' to turn this off too.
		{
			return false; ///< internal capacity overflow. (We want this in case we want to limit the number of conformers passed to the next stage.
		}
		else
		{
			// add the builders current conformer to the underlying buffer
			m_Conformers.push_back( m_Builder->getConformer() ); 
			return true;
		}
	}
}

