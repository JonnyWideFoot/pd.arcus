%include "std_string.i"
%include "exception.i"

#define PD_API

// Propagate Exceptions thrown in C++ to Python
%exception{
		try {
		$action
		} catch ( ExceptionBase ){
			SWIG_exception(SWIG_RuntimeError, "PD Exception occured" );
		}
}


%module arcus
%{
// Include the PD/MMLIB headers
#include "mmlib.h"

// Include the module header files
#include "arcus.h"

%}

// This import allows type information to be shared between the module and 
// the pd code. 
%import ../../pd/pd.i

%include "workspace/segdef.h"
%include "system/segbuilder.h"
%include "forcefields/segrejoin.h"
%include "filters/segdistfan.h"
%include "filters/surface.h"
%include "manipulators/segtorsions.h"
%include "manipulators/segmove.h"
%include "protocols/arcus_base.h"
%include "protocols/arcus_refiner.h"
%include "protocols/arcus_stitch.h"
%include "protocols/arcus_main.h"
%include "protocols/arcus_replay.h"

