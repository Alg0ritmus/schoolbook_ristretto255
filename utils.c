#include "ristretto255.h"
#include "utils.h"
#include "helpers.h"

void print(field_elem o){

    for (int i=0;i<FIELED_ELEM_SIZE;i++){
        printf("%llx ", o[i]);
        
    }
    printf("\n");
}


void print_32(const u8* o){

    for (int i=0;i<BYTES_ELEM_SIZE;i++){
        printf("%02hx ", o[i]);
        
    }
    printf("\n");
}

void pack_and_print_32(field_elem o){
    u8 temp[BYTES_ELEM_SIZE];
    pack25519(temp,o);
    print_32(temp);
}