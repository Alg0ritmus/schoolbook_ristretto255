#ifndef _HELPERS_H
#define _HELPERS_H

#include "ristretto255.h"
// helper functions (mostly tweeked tweetNaCl)
void unpack25519(field_elem out, const u8 *in);
void pack25519(u8 *out, const field_elem in);

// mine functions
void fcopy(field_elem out, const field_elem in);
void b_copy(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]);
int feq( const field_elem a,  const field_elem b); // return 1 if two are equal, otherwise 0
int bytes_eq_32( const u8 a[BYTES_ELEM_SIZE],  const u8 b[BYTES_ELEM_SIZE]); // return 1 if two are equal, otherwise 0
void fabsolute(field_elem out, const field_elem in);


void fneg(field_elem out, const field_elem in);
int is_neg(const field_elem in); // return 1 if it's negative
int is_neg_bytes(const u8 in[BYTES_ELEM_SIZE]);


#endif //_HELPERS_H