#include "helpers.h"
#include "tweetNaClLib.h"
#include "utils.h"


// tweetNaCl function, bcs. we are using Bernstein's representation of big numbers
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L335
void unpack25519(field_elem out, const u8 *in){
    int i;
    FOR(i,0,16) out[i] = in[2*i] + ((i64) in[2*i + 1] << 8);
    out[15] &= 0x7fff;
 }


// MEMORY ALLOCATION: 1x i64
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L273
// nie je to identikce, asi to keleppmann prepisal trosku 
void carry25519(field_elem elem){
    int i;
    i64 carry;
    FOR(i,0,16) {
        carry = elem[i] >> 16;
        elem[i] -= carry << 16;
        if (i < 15) elem[i + 1] += carry; else elem[0] += 38 * carry;
    }
}

// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L342
void fadd(field_elem out, const field_elem a, const field_elem b){ /* out = a + b */
    int i;
    FOR(i,0,16) out[i] = a[i] + b[i];
}


// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L348
void fsub(field_elem out, const field_elem a, const field_elem b){ /* out = a - b */
    int i;
    FOR(i,0,16) out[i] = a[i] - b[i];
}


// MEMORY ALLOCATION: 2x i64, 1x i64[31]
// TOTAL: 4x i64, 1x i64[31]
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L354
void fmul(field_elem out, const field_elem a, const field_elem b){ /* out = a * b */
    i64 i, j, product[31];
    FOR(i,0,31) product[i] = 0;
    FOR(i,0,16){
        FOR(j,0,16) product[i+j] += a[i] * b[j];
    }
    FOR(i,0,15) product[i] += 38 * product[i+16];
    FOR(i,0,16) out[i] = product[i];
    carry25519(out);
    carry25519(out);
}


// if bit = 1 -> swap, otherwise do nothing
// MEMORY ALLOCATION: 3x i64
// TOTAL: 3x i64
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L285
void swap25519(field_elem p, field_elem q, int bit){
    i64 t, i, c = ~(bit - 1);
    FOR(i,0,16) {
        t = c & (p[i] ^ q[i]);
        p[i] ^= t;
        q[i] ^= t;
    }
}


// MEMORY ALLOCATION: 2x field_elem
// TOTAL: 2x field_elem + 6x i64
// Same as unpack
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L295
void pack25519(u8 *out, const field_elem in){
    int i, j, carry;
    field_elem m, t;
    FOR(i,0,16) t[i] = in[i];
    carry25519(t); carry25519(t); carry25519(t);

    FOR(j,0,2) {
        m[0] = t[0] - 0xffed;
     
        FOR(i,1,15) {
            m[i] = t[i] - 0xffff - ((m[i-1] >> 16) & 1);
            m[i-1] &= 0xffff;
        }

        m[15] = t[15] - 0x7fff - ((m[14] >> 16) & 1);
        carry = (m[15] >> 16) & 1;
        m[14] &= 0xffff;
        swap25519(t, m, 1 - carry);
    }


    FOR(i,0,16) {
        out[2*i] = t[i] & 0xff;
        out[2*i + 1] = t[i] >> 8;
    }
}


/* ------------------------------------------------*/
/* -------------------my code----------------------*/
/* ------------------------------------------------*/

// return 1 if they're equal

// feq() => checks, if 2 field_elems are eq. in constant time
// MEMORY ALLOCATION: 2x u8[32]
// TOTAL: 2x u8[32] + 4x field_elem + 12x i64
int feq( const field_elem a,  const field_elem b){
    int result = 1;
    u8 a_32[BYTES_ELEM_SIZE],b_32[BYTES_ELEM_SIZE];
    pack25519(a_32,a);
    pack25519(b_32,b);

    // constant time
    result &= bytes_eq_32(a_32,b_32);

    // returns 1 if two are eq., 0 otherwise
    return result;
}

// return 1 if they're equal
// checking if two u8[32] are eq
int bytes_eq_32( const u8 a[BYTES_ELEM_SIZE],  const u8 b[BYTES_ELEM_SIZE]){
    int result = 1;

    for (int i = 0; i < BYTES_ELEM_SIZE; ++i){
        result &= a[i] == b[i];
    }

    return result;
}

// copy in to out
void fcopy(field_elem out, const field_elem in){
    int i;
    FOR(i,0,FIELED_ELEM_SIZE) out[i] = in[i];
}

void b_copy(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]){
    int i;
    FOR(i,0,BYTES_ELEM_SIZE) out[i]  = in[i];
}

// negation of element a -> p-a = -a 
// inspired by: https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_51.h#L94

// alebo zmenit na :
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_25_5.h#L116

// fneg() => creates negative input and "save" it to variable out
// MEMORY ALLOCATION: 1x i64
void fneg(field_elem out, const field_elem in){
    fsub(out, F_MODULUS, in);
    carry25519(out);
};


// Check if field_elem is negative, if LSB is set to 1, then it is negative, otherwise positive 
// inspired by: https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_25_5.h#L302
// MEMORY ALLOCATION: 1x u8[32]
// TOTAL: 1x u8[32] + 2x field_elem + 6x i64
int is_neg(const field_elem in){
    u8 temp[BYTES_ELEM_SIZE];
    pack25519(temp, in);
    return temp[0] & 1;
}

int is_neg_bytes(const u8 in[BYTES_ELEM_SIZE]){
    return in[0] & 1;
}


// fabsolute() => get absolute value of input in constant time 
// MEMORY ALLOCATION: 1x field_elem
// TOTAL: : 1x field_elem + 4x i64 + 1x u8[32]
void fabsolute(field_elem out, const field_elem in){
    field_elem temp;
    fcopy(temp,in); // temp=in, so I dont rewrite in
    fneg(out,temp); // out = ~in
    // CT_SWAP if it is neg.
    swap25519(out,temp,is_neg(in));
}


// square a^2
//MEMORY ALLOCATION: 4x i64, 1x i64[31]
void pow2(field_elem out, const field_elem a){
    fmul(out,a,a);
 }

 
// to the power of 3:  a^3

//MEMORY ALLOCATION: 8x i64, 2x i64[31]
void pow3(field_elem out, const field_elem a){
    pow2(out,a);
    // a^3
    fmul(out,out,a);
}

//MEMORY ALLOCATION: 12x i64, 3x i64[31]
void pow7(field_elem out, const field_elem a){
    pow3(out,a);
    pow2(out,out);
    fmul(out,out,a);
}

// In previous version of this library we used curve25519_pow_two5mtwo0_two250mtwo0 
// which was inspired by ristretto_dona: https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h
// Now we use pow2523 instead, which is implementation taken from tweetNaCl and slightly edited due to our field_elem representation.
// tweetNaCl impl.: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L382
void pow2523(field_elem o, const field_elem i){
    field_elem c;
    int a;
    FOR(a,0,16) c[a]=i[a];
    for(a=250;a>=0;a--) {
        pow2(c,c);
        if(a!=1) fmul(c,c,i);
    }
    FOR(a,0,16) o[a]=c[a];
}
