#ifndef _RISTRETTO255_H
#define _RISTRETTO255_H

#include <stdio.h>

typedef unsigned char u8;
typedef long long i64;
typedef i64 field_elem[16];

typedef struct ge_point25519{
	field_elem x,y,z,t;
}ristretto255_point;



extern const field_elem F_ZERO;
extern const field_elem F_ONE;
extern const field_elem F_TWO;
extern const field_elem F_BIGGEST;
extern const field_elem _121665;
extern const field_elem F_MODULUS;
extern const field_elem SQRT_M1;
extern const field_elem EDWARDS_D;


// helper functions (mostly tweeked tweetNaCl)
void unpack25519(field_elem out, const u8 *in);
void pack25519(u8 *out, const field_elem in);

// mine functions
void fcopy(field_elem out, const field_elem in);
void b_copy(u8 out[32], const u8 in[32]);
int feq( const field_elem a,  const field_elem b); // return 1 if two are equal, otherwise 0
int bytes_eq_32( const unsigned char a[32],  const unsigned char b[32]); // return 1 if two are equal, otherwise 0
void fabsolute(field_elem out, const field_elem in);


void fneg(field_elem out, const field_elem in);
int is_neg(const field_elem in); // return 1 if it's negative
int is_neg_bytes(const u8 in[32]);

int ristretto255_decode(ristretto255_point *ristretto_out, const unsigned char bytes_in[32]);
int ristretto255_encode(unsigned char bytes_out[32], const ristretto255_point *ristretto_in);
int hash_to_group(u8 bytes_out[32], const u8 bytes_in[64]);
void ristretto255_scalarmult(ristretto255_point* p, ristretto255_point* q,const u8 *s);



#endif //_HELPERS_H