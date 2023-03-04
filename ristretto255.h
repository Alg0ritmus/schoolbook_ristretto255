#ifndef _RISTRETTO255_H
#define _RISTRETTO255_H

#include <stdio.h>

typedef unsigned char u8;
typedef long long i64;
typedef i64 field_elem[16];



// CONSTANTS
extern const field_elem F_ZERO;
extern const field_elem F_ONE;
extern const field_elem F_BIGGEST;
extern const field_elem _121665;


void print(field_elem o);
void print_32(const u8* o);

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
void curve25519_pow_two252m3(field_elem two252m3, const field_elem z);
void inv_sqrt(field_elem out,const field_elem u, const field_elem v);
int feq( const field_elem a,  const field_elem b); // return 1 if two are equal, otherwise 0


#endif //_RISTRETTO255_H