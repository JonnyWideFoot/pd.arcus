#include "global.h"

#include "elements.h"

namespace Library
{
	// ----------------------------------------------------------------------------------
	// --- The Elements
	// ----------------------------------------------------------------------------------

	// incomplete, stops before Lanthanoids, but hey .. :)
	const char *element[] = { "-",
		"H", "He",
		"Li","Be", "B" ,"C" , "N", "O", "F","Ne",
		"Na","Mg", "Al","Si", "P", "S","Cl","Ar",
		"K" ,"Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn", "Ga","Ge","As","Se","Br","Kr",
		"Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", "In","Sn","Sb","Te","I", "Xe",
		"Cs","Ba" };

	int ElementSymbol2Z(const std::string &Symbol){
		char buffer1[15];
		char buffer2[15];
		strncpy(&buffer1[0],Symbol.c_str(),15);
		strlc(&buffer1[0]);
		for(int i=0;i<56;i++){
			strcpy(&buffer2[0],element[i]);
			strlc(&buffer2[0]);
			if(strcmp(&buffer1[0],&buffer2[0])==0){
				return i;
			}
		}
		return -1;
	}
}

