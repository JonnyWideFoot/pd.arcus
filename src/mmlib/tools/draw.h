#ifndef __DRAW_TOOLS_H
#define __DRAW_TOOLS_H

/// \brief Draw the origin (0,0,0) as a Point via the IDrawProvider
void drawOrigin( IDrawProvider& _Vect, int colourCode );

/// \brief Draw the given cartesian coordinate as a Point via the IDrawProvider
void drawPoint( IDrawProvider& _Vect, const Maths::dvector& pos, int colourCode );

#endif

