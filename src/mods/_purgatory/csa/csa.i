%include "std_string.i"
%include "exception.i"
#define pd_API
// Propagate Exceptions thrown in C++ to Python
%exception{
		try {
		$action
		} catch ( ExceptionBase ){
			SWIG_exception(SWIG_RuntimeError, "pd Exception occured" );
		}
}


%module csa
%{
#include "global.h"
#include "csa.h"
%}
%include csa.h


