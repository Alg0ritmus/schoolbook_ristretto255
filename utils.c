#include "ristretto255.h"
#include "utils.h"
void print(field_elem o){

	for (int i=0;i<16;i++){
		printf("%llx ", o[i]);
		
	}
	printf("\n");
}


void print_32(const u8* o){

	for (int i=0;i<32;i++){
		printf("%02hx ", o[i]);
		
	}
	printf("\n");
}

void pack_and_print_32(field_elem o){
	u8 temp[32];
	pack25519(temp,o);
	print_32(temp);
}