#include "Utilities.h"

// Source: https://stackoverflow.com/a/2736841
char *remove_ext(char* mystr, char dot, char sep) {
	char *retstr, *lastdot, *lastsep;

	// Error checks and allocate string.

	if (mystr == NULL)
		return NULL;
	if ((retstr = (char*)malloc(strlen(mystr) + 1)) == NULL)
		return NULL;

	// Make a copy and find the relevant characters.

	strcpy(retstr, mystr);
	lastdot = strrchr(retstr, dot);
	lastsep = (sep == 0) ? NULL : strrchr(retstr, sep);

	// If it has an extension separator.

	if (lastdot != NULL) {
		// and it's before the extenstion separator.

		if (lastsep != NULL) {
			if (lastsep < lastdot) {
				// then remove it.

				*lastdot = '\0';
			}
		}
		else {
			// Has extension separator with no path separator.

			*lastdot = '\0';
		}
	}

	// Return the modified string.

	return retstr;
}

// Source: https://stackoverflow.com/a/12911465
const char* get_field(char* line, int num) {
    const char* tok;
    for(tok = strtok(line, ";"); tok && *tok; tok =strtok(NULL, ";\n")) {
        if(!--num) {
            return tok;
        }
    }

    return NULL;
}
