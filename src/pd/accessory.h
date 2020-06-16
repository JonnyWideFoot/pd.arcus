#ifndef __ACCESSORY_H
#define __ACCESSORY_H

#include <stdlib.h>
#include <string>

/// Print the introductory application header
void PrintFullpdHeader();
void Printpd();
void PrintQuote();

/// Python wrapper for the Random Number Seed used by the C Standard Library
/// sets a specific seed
void cseed(unsigned int newseed);

/// Python wrapper for the Random Number Seed used by the C Standard Library
/// sets the seed to the system time
void ctimeseed();

void cprint(std::string text);

/// Full header plus info()
void detail(); 

/// prints Run Time & Compile Time Info 
void info(); 

void timer(bool reset=false);		

// OS dependent functions
void printTimeStamp();
void printHostname();
void printWorkingDirectory();
void printTitle(const std::string &_Title);

class  classA{
	public:

	 int test(){ printf("A\n"); };
};


class classB{
	public:

	 int test(){ printf("B\n"); };
#ifndef SWIG	
	operator classA& (){
		return my_a;
	}
	operator int (){
		return 4;
	}
#endif
	classA my_a;
};

#ifdef SWIG
 %implicitconv classA;
 %implicitconv classB;
#endif


class classC{
	public:
		classC( const classA & _a ): 
		a(_a){  };

	classA a;
	classB b;
};

class printer{
public:
void printint( int myint ){
	printf("--> %d \n", myint );
}

};

/*
void test(){
	A a;
	B b;

	C c1(a);
	C c2(b);

}
*/

#endif
