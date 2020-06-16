#ifndef __CONSTANTS_H
#define __CONSTANTS_H

// these can be directly casted to integers for interpretation... i.e. red is 0 in DAVE
class Colour
{
	private:
		Colour();
	public:
		static const int Red;
		static const int Blue;
		static const int Yellow;
		static const int Green;
		static const int Cyan;
		static const int Magenta;
		static const int Orange;
		static const int Turquoise;
};

#endif

