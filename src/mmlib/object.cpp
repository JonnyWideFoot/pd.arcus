#include "global.h"

long Object::object_nextid = 0;

Object::Object()
{
	object_id = object_nextid;
	object_nextid++;
	setInternalName();
}

Object::Object( const Object & _clone )
{
	(*this) = _clone;
}

Object & Object::operator= ( const Object & _clone )
{
	if( &_clone == this ) return *this; // Handle self-assignement
	name = _clone.name ;							 
	object_id = object_nextid;
	object_nextid++;
	return (*this);
}

bool operator== (const Object &a, const Object &b)
{
	return a.object_id == b.object_id;
}

void Object::setInternalName()
{
	name = "Object";
}

