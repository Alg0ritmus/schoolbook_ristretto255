#ifndef _TWEETNACLLIB_H
#define _TWEETNACLLIB_H



/*------------ CONSTANTS ------------*/

// little-endian order --> a = a0*2^0 + a1*2^16 + a2*2^32 + ... + a15*2^240 (vzdy o ax*2^ o15 vyssie)
// conversion 2^255 - 19 is done throgh PY script Python_impl/convertLib.py
// CONSTANTS taken from ristretto draft: https://datatracker.ietf.org/doc/draft-hdevalence-cfrg-ristretto/
static const field_elem F_ZERO = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    F_ONE = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    F_MODULUS = {0xffed, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0x7fff},
    SQRT_M1 = {0xa0b0, 0x4a0e, 0x1b27, 0xc4ee,0xe478, 0xad2f, 0x1806, 0x2f43,0xd7a7, 0x3dfb, 0x99,  0x2b4d,0xdf0b, 0x4fc1, 0x2480, 0x2b83},
    EDWARDS_D = {0x78a3, 0x1359, 0x4dca, 0x75eb, 0xd8ab, 0x4141, 0xa4d, 0x70, 0xe898, 0x7779, 0x4079, 0x8cc7, 0xfe73, 0x2b6f, 0x6cee, 0x5203},
    EDWARDS_D2 = {0xf159, 0x26b2, 0x9b94, 0xebd6, 0xb156, 0x8283, 0x149a, 0xe0, 0xd130, 0xeef3, 0x80f2, 0x198e, 0xfce7, 0x56df, 0xd9dc, 0x2406},
    INVSQRT_A_MINUS_D  = {0x40ea, 0x805d, 0xfdaa, 0x99c8, 0x72be, 0x5a41, 0x1617, 0x9d2f, 0xd840, 0xfe01, 0x7b91, 0x16c2, 0xfca2, 0xcfaf, 0x8905, 0x786c},
    ONE_MINUS_D_SQ = {0xc176, 0x945f, 0x9c1, 0xe27c, 0x350f, 0xcd5e, 0xa138, 0x2c81, 0xdfe4, 0xbe70, 0xabdd, 0x9994, 0xe0d7, 0xb2b3, 0x72a8, 0x290},
    D_MINUS_ONE_SQ = {0x4d20, 0x44ed, 0x5aaa, 0x31ad, 0x1999, 0xb01e, 0x4a2c, 0xd29e, 0x4eeb, 0x529b, 0xd32f, 0x4cdc, 0x2241, 0xf66c, 0xb37a, 0x5968},
    SQRT_AD_MINUS_ONE = {0x2e1b, 0x497b, 0xf6a0, 0x7e97, 0x54bd, 0x1b78, 0x8e0c, 0xaf9d, 0xd1fd, 0x31f5, 0xfcc9, 0xf3c, 0x48ac, 0x2b83, 0x31bf, 0x3769};


/*-----------------------------------*/

void unpack25519(field_elem out, const u8 *in);
void carry25519(field_elem elem);
void fadd(field_elem out, const field_elem a, const field_elem b); /* out = a + b */
void fsub(field_elem out, const field_elem a, const field_elem b); /* out = a - b */
void fmul(field_elem out, const field_elem a, const field_elem b); /* out = a * b */
void swap25519(field_elem p, field_elem q, int bit);
void pack25519(u8 *out, const field_elem in);
int feq( const field_elem a,  const field_elem b);
int bytes_eq_32( const u8 a[BYTES_ELEM_SIZE],  const u8 b[BYTES_ELEM_SIZE]);
void fcopy(field_elem out, const field_elem in);
void b_copy(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]);
void fneg(field_elem out, const field_elem in);
int is_neg(const field_elem in);
int is_neg_bytes(const u8 in[BYTES_ELEM_SIZE]);
void fabsolute(field_elem out, const field_elem in);
void pow2(field_elem out, const field_elem a);
void pow3(field_elem out, const field_elem a);
void pow7(field_elem out, const field_elem a);
void pow2523(field_elem o, const field_elem i);

#endif //_TWEETNACLLIB_H