#ifndef _RISTRETTO255_H
#define _RISTRETTO255_H

#include <stdio.h>

typedef unsigned char u8;
typedef long long i64;
typedef i64 field_elem[16];

typedef struct ge_point25519{
	field_elem x,y,z,t;
}ristretto255_point;



// CONSTANTS
extern const field_elem F_ZERO;
extern const field_elem F_ONE;
extern const field_elem F_TWO;
extern const field_elem F_BIGGEST;
extern const field_elem _121665;
extern const field_elem F_MODULUS;
extern const field_elem SQRT_M1;
extern const field_elem EDWARDS_D;


void print(field_elem o);
void print_32(const u8* o);
void pack_and_print_32(field_elem o);

void unpack25519(field_elem out, const u8 *in);
void carry25519(field_elem elem);
void fadd(field_elem out, const field_elem a, const field_elem b);
void fsub(field_elem out, const field_elem a, const field_elem b);
void fmul(field_elem out, const field_elem a, const field_elem b);/* out = a * b */
void finverse(field_elem out, const field_elem in);
void swap25519(field_elem p, field_elem q, int bit);
void pack25519(u8 *out, const field_elem in);
void scalarmult(u8 *out, const u8 *scalar, const u8 *point);

// mine funkcie
void fcopy(field_elem out, const field_elem in);
void fneg(field_elem out, const field_elem in);
void pow2(field_elem out, const field_elem a);
void pow3(field_elem out, const field_elem a);
void pow7(field_elem out, const field_elem a);
void pow_xtimes(field_elem out, const field_elem a, int n);
void curve25519_pow_two5mtwo0_two250mtwo0(field_elem b);
// je mozne nahradit tymto kodom (sv pow2523)?
// https://github.com/sbp/tweetnacl-tools/blob/master/tweetnacl.c#L382 
void curve25519_pow_two252m3(field_elem two252m3, const field_elem z);
void inv_sqrt(field_elem out,const field_elem u, const field_elem v);
int feq( const field_elem a,  const field_elem b); // return 1 if two are equal, otherwise 0
int bytes_eq_32( const unsigned char a[32],  const unsigned char b[32]); // return 1 if two are equal, otherwise 0

// funguje???
void fneg(field_elem out, const field_elem in);
int is_neg(const field_elem in); // return 1 if it's negative

// is negative ?
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_51.h#L243

// ristretto functions

int ristretto255_decode(ristretto255_point *ristretto_out, const unsigned char bytes_in[32]);
int ristretto255_encode(unsigned char bytes_out[32], const ristretto255_point *ristretto_in);

#endif //_RISTRETTO255_H