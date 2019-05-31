#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *remove_ext(char* mystr, char dot, char sep);
const char* get_field(char* line, int num);
