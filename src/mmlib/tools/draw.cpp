#include "global.h"
#include "draw.h"

void drawOrigin( IDrawProvider& _Vect, int colourCode )
{
	Maths::dvector origin(0.0,0.0,0.0);
	drawPoint( _Vect, origin, colourCode );
}

void drawPoint( IDrawProvider& _Vect, const Maths::dvector& pos, int colourCode )
{
	const float ptSize = 0.2f;
	DrawingVector* dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.x -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.x += ptSize;
		dv->colourCode = colourCode;
	}
	dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.y -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.y += ptSize;
		dv->colourCode = colourCode;
	}
	dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.z -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.z += ptSize;
		dv->colourCode = colourCode;
	}
}

