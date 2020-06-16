#include "global.h"

#include "tools/rdstdout.h"

char *rdstdout_file = NULL;

void set_rdstdout_file(char *filename){
	rdstdout_file = filename;
}

void clearrdstdout(){
	if(rdstdout_file != NULL) {
		FILE *rdfile = fopen(rdstdout_file, "w");
		if(rdfile == NULL)
			return;
		fclose(rdfile);
	}
}

void rdprintf(const char *a){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
