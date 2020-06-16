#ifndef __STREAM_WRITER
#define __STREAM_WRITER

#include <string>
#include <fstream>






//-------------------------------------------------
//
/// \brief  Defines a class that will write to an arbritrary stream.
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
class StreamWriter
{
public:
	StreamWriter();
	StreamWriter( const std::string& _filename, bool _Append = false );
	StreamWriter( std::ostream& _stream, bool _Append = false );

	void streamToScreen();
	void streamTo( const std::string& _filename, bool _Append = false, bool _ClearFile = true );
	void streamTo( std::ostream& _stream );

	void assertStream() const ; ///< Throws a code exception if the underlying stream is not open.
	bool isStreaming() const { return m_Streaming; } ///< Resturns wether the current internal stream is open and ready for business.

protected:
	void beginStream();
	void endStream();
	std::ostream& stream() const { return *m_Stream; }

private:
	enum InternalStream
	{
		Screen,
		File,
		Stream
	};

	InternalStream m_Mode; /// Keep tags on what we are actually writing to
	bool m_Streaming; ///< Is a streaming operation currently underway? Controlled by the beginStream() and endStream() calls.
	bool m_Appending; ///< Should we append or open fresh? - applies to the file stream only.
	std::string m_FileName; ///< If we are writing to a file, whats it called?
	std::ofstream* m_FileStream; ///< File stream is opened once per transaction - Appending mode is required if previous transactions are not to be over-written
	std::ostream* m_Stream; ///< Pointer to the current stream.
};

#endif

