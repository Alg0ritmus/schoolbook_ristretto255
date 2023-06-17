#ifndef _UTILS_H
#define _UTILS_H

#include "ristretto255.h"
#include "helpers.h"
/*		UTILS		*/
void print(field_elem o);

void print_32(const u8* o);

void pack_and_print_32(field_elem o);
#endif //_UTILS_H

