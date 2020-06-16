#ifndef __INTERFACE_H
#define __INTERFACE_H

// Defines a list of public ABC interfaces 
// that can be used throughout PD

// Forward declarations
struct DrawingVector; // defined in primitives.h

/// \brief  IDrawProvider
/// \details Defiens a class interface for those which can draw lines for user feedback
/// \author  Jon Rea 
class IDrawProvider
{
public:
	virtual DrawingVector* request() = 0;
};

/// \brief ICloneable - any clonable class should define itself as ICloneable
/// \details ICloneable - any clonable class should define itself as ICloneable
/// \author  Jon Rea 
class ICloneable
{
public:
	virtual ICloneable* clone() const = 0;
};

#endif

