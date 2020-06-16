#include "global.h"

#include "exception.h"

using namespace std;

ExceptionBase::ExceptionBase(const std::string &_message) :
	std::exception(),
	message(_message)
{
	Details();
}

ExceptionBase::~ExceptionBase() throw()
{
}

void ExceptionBase::Details()
{
	if( message.size() == 0 ) message = "No detail available";
	  std::cerr << std::endl << "\n-----! PD Exception !----------------------------" << std::endl;
	  std::cerr << "Description: " << message << std::endl;
}

NotImplementedException::NotImplementedException(const std::string &_message) :
	ExceptionBase( _message )
{
}

NotImplementedException::~NotImplementedException() throw()
{
}


OutOfMemoryException::OutOfMemoryException(const std::string &_message) :
	ExceptionBase( _message )
{
}

OutOfMemoryException::~OutOfMemoryException() throw()
{
}


CodeException::CodeException(const std::string &_message ) :
	ExceptionBase( _message )
{
}

CodeException::~CodeException() throw()
{
}



ArgumentException::ArgumentException(const std::string &_message) :
	ExceptionBase( _message )
{
}

ArgumentException::~ArgumentException() throw()
{
}


NullArgumentException::NullArgumentException(const std::string &_message) :
	ArgumentException( _message )
{
}

NullArgumentException::~NullArgumentException() throw()
{
}

NullInternalException::NullInternalException(const std::string &_message) :
	ArgumentException( _message )
{
}

NullInternalException::~NullInternalException() throw()
{
}


ProcedureException::ProcedureException(const std::string &_message) :
	ExceptionBase( _message )
{
}

ProcedureException::~ProcedureException() throw()
{
}



IOException::IOException(const std::string &_message) :
	ExceptionBase( _message )
{
}

IOException::~IOException() throw()
{
}




ParseException::ParseException(const std::string &_message) :
	ExceptionBase( _message )
{
}

ParseException::ParseException(
	const std::string &_message, 
	const std::string &_filename,
	const int          _linenumber):
	ExceptionBase( _filename +   ( _linenumber < 0 ? " " : ( std::string( ":" + int2str(_linenumber) ) ) ) + "   Error: " + _message ) 
{
}

void	PrintFileLineError(const std::string &_message, 
								 				 const std::string &_filename,
								 				 const int  _linenumber )
{
	
	std::string msg( _filename +   ( _linenumber < 0 ? " " : ( std::string( ":" + int2str(_linenumber) ) ) ) + "   Error: " + _message ); 
	printf("%s\n",msg.c_str() );
}


void	PrintFileLineWarning(const std::string &_message, 
								 				 const std::string &_filename,
								 				 const int  _linenumber )
{
	
	std::string msg( _filename +   ( _linenumber < 0 ? " " : ( std::string( ":" + int2str(_linenumber) ) ) ) + "   Warning: " + _message ); 
	printf("%s\n",msg.c_str() );
}



ParseException::~ParseException() throw()
{
}



UninitialisedException::UninitialisedException(const std::string &_message) :
	ExceptionBase( _message )
{
}

UninitialisedException::~UninitialisedException() throw()
{
}



OutOfRangeException::OutOfRangeException(const std::string &_message) :
	ExceptionBase( _message )
{
}

OutOfRangeException::~OutOfRangeException() throw()
{
}

void testpdException(){
	throw(CodeException("TEST"));
}


