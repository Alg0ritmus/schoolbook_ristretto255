#ifndef _RISTRETTO255_H
#define _RISTRETTO255_H

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#define FIELED_ELEM_SIZE 16 // i64 field_elem[16];
#define BYTES_ELEM_SIZE 32 // byte representation of field elements u8[32]
#define HASH_BYTES_SIZE 2*BYTES_ELEM_SIZE // u8[64]

typedef uint8_t  u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

typedef i64 field_elem[FIELED_ELEM_SIZE];

typedef struct ge_point25519{
    field_elem x,y,z,t;
}ristretto255_point;

int ristretto255_decode(ristretto255_point *ristretto_out, const unsigned char bytes_in[BYTES_ELEM_SIZE]);
int ristretto255_encode(unsigned char bytes_out[BYTES_ELEM_SIZE], const ristretto255_point *ristretto_in);
int hash_to_group(u8 bytes_out[BYTES_ELEM_SIZE], const u8 bytes_in[HASH_BYTES_SIZE]);
void ristretto255_scalarmult(ristretto255_point* p, ristretto255_point* q,const u8 *s);

void inverse_mod_l(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]);



#endif //_RISTRETTO255_H