#include "ristretto255.h"
#include "helpers.h"
#include "utils.h"

/*-----------------------------------*/

#define FOR_T(i, start, end) for (i = (start); i < (end); i++)
#define FOR(i, start, end)         FOR_T(i, start, end)
#define COPY(i, dst, src, size)    FOR(i, 0, size) (dst)[i] = (src)[i]
#define ZERO(i, buf, size)         FOR(i, 0, size) (buf)[i] = 0

/*-----------------------------------*/


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

// tweetNaCl function, bcs. we are using Bernstein's representation of big numbers
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L335
void unpack25519(field_elem out, const u8 *in){
    int i;
    FOR(i,0,16) out[i] = in[2*i] + ((i64) in[2*i + 1] << 8);
    out[15] &= 0x7fff;
 }


// MEMORY ALLOCATION: 1x i64
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L273
static void carry25519(field_elem elem){
    int i;
    i64 carry;
    FOR(i,0,16) {
        carry = elem[i] >> 16;
        elem[i] -= carry << 16;
        if (i < 15) elem[i + 1] += carry; else elem[0] += 38 * carry;
    }
}

// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L342
static void fadd(field_elem out, const field_elem a, const field_elem b){ /* out = a + b */
    int i;
    FOR(i,0,16) out[i] = a[i] + b[i];
}


// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L348
static void fsub(field_elem out, const field_elem a, const field_elem b){ /* out = a - b */
    int i;
    FOR(i,0,16) out[i] = a[i] - b[i];
}


// MEMORY ALLOCATION: 2x i64, 1x i64[31]
// TOTAL: 4x i64, 1x i64[31]
// Taken from TweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L354
static void fmul(field_elem out, const field_elem a, const field_elem b){ /* out = a * b */
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
static void swap25519(field_elem p, field_elem q, int bit){
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
static void pow2(field_elem out, const field_elem a){
    fmul(out,a,a);
 }

 
// to the power of 3:  a^3

//MEMORY ALLOCATION: 8x i64, 2x i64[31]
static void pow3(field_elem out, const field_elem a){
    pow2(out,a);
    // a^3
    fmul(out,out,a);
}

//MEMORY ALLOCATION: 12x i64, 3x i64[31]
static void pow7(field_elem out, const field_elem a){
    pow3(out,a);
    pow2(out,out);
    fmul(out,out,a);
}

// In previous version of this library we used curve25519_pow_two5mtwo0_two250mtwo0 
// which was inspired by ristretto_dona: https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h
// Now we use pow2523 instead, which is implementation taken from tweetNaCl and slightly edited due to our field_elem representation.
// tweetNaCl impl.: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L382
static void pow2523(field_elem o, const field_elem i){
    field_elem c;
    int a;
    FOR(a,0,16) c[a]=i[a];
    for(a=250;a>=0;a--) {
        pow2(c,c);
        if(a!=1) fmul(c,c,i);
    }
    FOR(a,0,16) o[a]=c[a];
}

// function that calc ristretto255 inverse square root
// Inspired by: https://ristretto.group/formulas/invsqrt.html
// MEMORY ALLOCATION: 3x field_elem
// TOTAL: 59x i64 + 13x i64[31] + 6x u8[32] + 16x field_elem + undefined (curve25519_pow_two252m3)
static int inv_sqrt(field_elem out,const field_elem u, const field_elem v){
    field_elem tmp1, tmp2, tmp3;

    int r_is_negative;
    int correct_sign_sqrt;
    int flipped_sign_sqrt;
    int flipped_sign_sqrt_i;
    int was_nonzero_square;
    int should_rotate;


    #define v3 tmp1
    #define v7 tmp2
    pow3(v3,v);                             //v^3
    pow7(v7,v);                             //v^7

    #define st_bracket tmp3
    fmul(st_bracket,u,v3);                  // (u*v^3)
    #define nd_bracket tmp1
    fmul(nd_bracket,u,v7);                  // (u*v^7)

    
    #define r tmp2
    pow2523(r,nd_bracket);                  // r = (u*v^7) ^ {(p-5)/8}
    #define r2 tmp1
    fmul(r2,r,st_bracket);                  // r2 = (u*v^3) * (u*v^7) ^ {(p-5)/8}


    
    #define temporary tmp2
    pow2(temporary,r2);                         // tmp = r2 ^ 2 -> needed for check
    #define check tmp3
    fmul(check,v,temporary);                    // check = (r2 ^ 2) * v

    
    #define u_neg tmp2

    fneg(u_neg,u);

    correct_sign_sqrt = feq(check, u);
    flipped_sign_sqrt = feq(check, u_neg);

    fmul(u_neg,u_neg,SQRT_M1); // u_neg_i = u_neg * SQRT_M1
    flipped_sign_sqrt_i = feq(check, u_neg);


    
    #define r_prime tmp3
    fmul(r_prime, r2, SQRT_M1);
    should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i;
    swap25519(r2, r_prime, should_rotate);

    
    // Choose the non-negative square root

    #define r_negative tmp3
    r_is_negative = is_neg(r2);
    fneg(r_negative, r2);

    swap25519(r2, r_negative, r_is_negative);

    was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt;


    // output
    fcopy(out,r2);
    return was_nonzero_square;

}


// MAP function from draft
// MAP is used in hash_to_group() fuction to get ristretto point from hash
// Input parameter is field_elem "t"
// Output is point on Edward's (interpreted as ristretto255_point) curve with coords X,Y,Z,T
// MAP is also called ristretto255_elligator

// MEMORY ALLOCATION: 5x field_elem
// TOTAL: 14x field_elem + 70x i64, 15x i64[31] + 5x u8[32] 
static int MAP(ristretto255_point* ristretto_out, const field_elem t){ 
    field_elem tmp1, tmp2, tmp3, tmp4, tmp5;
    int was_square, wasnt_square;

    #define _r tmp1
    #define u tmp2
    #define c tmp3
    #define rpd tmp4
    #define v tmp5
    fmul(_r,t,t);
    fmul(_r,SQRT_M1,_r);
    fadd(u,_r,F_ONE);
    fmul(u,u,ONE_MINUS_D_SQ);
    fneg(c,F_ONE);

    fadd(rpd,_r,EDWARDS_D);
    fmul(v,_r,EDWARDS_D);
    fsub(v,c,v);
    fmul(v,v,rpd);

    #define s tmp4
    #define s_prime tmp2
    was_square = inv_sqrt(s,u,v);
    wasnt_square = 1 - was_square;


    fmul(s_prime,s,t);
    fabsolute(s_prime,s_prime);
    fneg(s_prime,s_prime);

    swap25519(s,s_prime,wasnt_square);
    #define rc tmp2
    fcopy(rc,_r);
    swap25519(c,rc,wasnt_square);

    #define n tmp2
    fsub(n,_r,F_ONE);
    fmul(n,n,c);
    fmul(n,n,D_MINUS_ONE_SQ);
    fsub(n,n,v);

    #define ss tmp1
    #define w0 tmp3
    fmul(ss,s,s);

    fadd(w0,s,s);
    fmul(w0,w0,v);
    #define w1 tmp4
    fmul(w1,n,SQRT_AD_MINUS_ONE);
    #define w2 tmp2
    #define w3 tmp5
    fsub(w2,F_ONE,ss);
    fadd(w3,F_ONE,ss);

    fmul(ristretto_out->x,w0,w3);
    fmul(ristretto_out->y,w2,w1);
    fmul(ristretto_out->z,w1,w3);
    fmul(ristretto_out->t,w0,w2);

    return 0;
}




// Ristretto255 point addition (Edward's coords | X,Y,Z,T ) 
// ge25519_addition is identical to tweetNaCl sv add(gf p[4],gf q[4]),
// Basically I changed name to smth more slef-explaining and chaged gf[4] to our ristretto255_point representation
// of ristretto points. NOTE: gf is the same as our field_elem (just naming is different)
// https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L590
void ge25519_addition(ristretto255_point* r,const ristretto255_point* p,const ristretto255_point* q){
    field_elem temp_1,temp_2,temp_3, temp_4, temp_5;

    #define a temp_1
    #define _t temp_2
    #define b temp_3
    
    fsub(a, p->y, p->x);
    fsub(_t, q->y, q->x);
    fmul(a, a, _t);
    fadd(b, p->x, p->y);
    fadd(_t, q->x, q->y);
    fmul(b, b, _t);

    #define e temp_5
    #define h temp_2
    fadd(h, b, a);
    fsub(e, b, a);
    
    
    #define d temp_1
    #define _c temp_3
    fmul(d, p->z, q->z);
    fadd(d, d, d);

    fmul(_c, p->t, q->t);
    fmul(_c, _c, EDWARDS_D2);
    #define f temp_4
    fsub(f, d, _c);
    #define g temp_3
    fadd(g, d, _c);

    
    

    fmul(r->x, e, f);
    fmul(r->y, h, g);
    fmul(r->z, g, f);
    fmul(r->t, e, h);
}

// cswap() => conditional swap, if b is set to 1, then swap 2 ristretto points,
// inspired by tweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L615
static void cswap(ristretto255_point* p, ristretto255_point* q,u8 b){
    swap25519(p->x,q->x,b);
    swap25519(p->y,q->y,b);
    swap25519(p->z,q->z,b);
    swap25519(p->t,q->t,b);
}


// Inspired by ristretto draft
// generate ristretto point from bytes[32]
// MEMORY ALLOCATION: 6x field_elem, 11x u8 checked_bytes[32]
// TOTAL: 31x field_elem + 1x u8 checked_bytes[32] + 109x i64 + 23x i64[31] + undefined part of inv_sqrt (depends on input vlaues)
int ristretto255_decode(ristretto255_point *ristretto_out, const u8 bytes_in[BYTES_ELEM_SIZE]){
  
  int was_square, is_canonical, is_negative;

  field_elem temp1,temp2,temp3,temp4,temp5,temp6;

  u8 checked_bytes[BYTES_ELEM_SIZE];

  // Step 1: Check that the encoding of the field element is canonical
  #define _s temp1
  unpack25519(_s, bytes_in);
  pack25519(checked_bytes,_s);
  

  // check if bytes_in == checked_bytes, else abort
  is_canonical = bytes_eq_32(checked_bytes,bytes_in);
  is_negative = is_neg_bytes(bytes_in);

  if (is_canonical == 0 || is_negative==1){
    #ifdef DDEBUG_FLAG
        printf("ristretto255_decode: Bad encoding or neg bytes passed to the ristretto255_decode function! is_canonical=%d, is_negative=%d\n", is_canonical, is_negative);
    #endif
    return 1;
  }

  // Step 2 calc ristretto255/ge25519 point
  // a = ± 1
  #define _ss temp2
  #define u1 temp3
  #define u2 temp4
  

  pow2(_ss,_s);
  fsub(u1,F_ONE,_ss); // 1 + as^2
  fadd(u2,F_ONE,_ss); // 1 - as^2

  
  #define uu1 temp5
  pow2(uu1,u1);
  


  // d -> from ristretto darft
  #define duu1_positive temp6
  fmul(duu1_positive,EDWARDS_D,uu1); // d*u1^2
  #define duu1 temp5
  fneg(duu1,duu1_positive); // -(D * u1^2) 

  #define uu2 temp6
  pow2(uu2,u2);
  #define _v temp2
  fsub(_v,duu1,uu2); // adu1^2 - u2^2
  #define vuu2 temp5
  fmul(vuu2,_v,uu2); // v*u2^2

  #define _I temp6
  was_square = inv_sqrt(_I,F_ONE,vuu2);

  #define Dx temp5
  fmul(Dx,_I,u2); // I*u2
  #define Dxv temp4
  fmul(Dxv, Dx, _v); // Dx * v


  // x: |2sDx|
  #define sDx temp2
  fmul(sDx,_s,Dx);
  fadd(ristretto_out->x,sDx,sDx); // 2sDx
  fabsolute(ristretto_out->x,ristretto_out->x);// |2sDx|

  #define Dy temp5
  fmul(Dy, _I, Dxv); // I*Dx*v

  // y: u1Dy
  fmul(ristretto_out->y,u1,Dy);

  // t: x,y
  fmul(ristretto_out->t,ristretto_out->x,ristretto_out->y);

  fcopy(ristretto_out->z, F_ONE);

  if (was_square == 0){
    #ifdef DEBUG_FLAG
      printf("\n\n\n ristretto255_decode: Bad encoding! was_square=%d \n\n\n",was_square);
    #endif
    return 1;
  }
  if (is_neg(ristretto_out->t)){
    #ifdef DEBUG_FLAG
      printf("\n\n\n ristretto255_decode: Bad encoding! is_neg(t)=%d \n\n\n",is_neg(t));
    #endif
    return 1;
  }
  if (feq(ristretto_out->t,F_ZERO)){
    #ifdef DEBUG_FLAG
      printf("\n\n\n ristretto255_decode: Bad encoding! feq(y,F_ZERO)=%d \n\n\n",feq(y,F_ZERO));
    #endif
    return 1;
  }
  return 0;
}


// Inspired by ristretto draft
// MEMORY ALLOCATION: 7x field_elem
// TOTAL: : 10x field_elem + 79x i64, 14x i64[31] + 3x u8[32]
int ristretto255_encode(unsigned char bytes_out[BYTES_ELEM_SIZE], const ristretto255_point* ristretto_in){

  
  field_elem _temp1,_temp2,_temp3,_temp4,_temp5,_temp6,_temp7;


  #define temp_zy1 _temp1
  #define temp_zy2 _temp2
  #define u1_ _temp3
  fadd(temp_zy1,ristretto_in->z,ristretto_in->y);
  fsub(temp_zy2,ristretto_in->z,ristretto_in->y);
  fmul(u1_,temp_zy1,temp_zy2); // (Z+Y)(Z-Y)

  #define u2_ _temp1
  fmul(u2_,ristretto_in->x,ristretto_in->y); // XY

  #define uu2_ _temp2
  pow2(uu2_,u2_); //u2^2 
  #define u1uu2 _temp4
  fmul(u1uu2,u1_,uu2_); // u1 * u2^2
  #define I_ _temp2
  inv_sqrt(I_,F_ONE,u1uu2);

  #define D1_ _temp4
  fmul(D1_,u1_,I_); // u1 * I
  #define D2_ _temp3
  fmul(D2_,u2_,I_); // u2 * I

  #define D1D2 _temp1
  fmul(D1D2,D1_,D2_); // D1 * D2
  #define Zinv _temp2
  fmul(Zinv,D1D2,ristretto_in->t); // Zinv = D1*D2*T
  

  #define enchanted_denominator _temp1
  fmul(enchanted_denominator,D1_,INVSQRT_A_MINUS_D);
  #define tZinv _temp4
  fmul(tZinv,ristretto_in->t,Zinv);
  
  int is_tZinv_neg = 1-is_neg(tZinv);
  #define _X _temp4
  #define _Y _temp5
  fcopy(_X,ristretto_in->x);
  fcopy(_Y,ristretto_in->y);

  #define iX _temp6
  #define iY _temp7
  fmul(iX, ristretto_in->x, SQRT_M1);
  fmul(iY, ristretto_in->y, SQRT_M1);
  swap25519(iY,_X,is_tZinv_neg);
  swap25519(iX,_Y,is_tZinv_neg);
  swap25519(enchanted_denominator,D2_,is_tZinv_neg);

  

  #define XZ_inv _temp3
  fmul(XZ_inv,iY,Zinv);
  #define n_Y _temp2
  fneg(n_Y,iX);
  fcopy(iX,iX);
  swap25519(iX,n_Y,is_neg(XZ_inv));
  

  #define _Z _temp3
  fcopy(_Z,ristretto_in->z);

  #define Z_Y _temp2
  fsub(Z_Y,_Z,iX);

  #define temp_s _temp3
  fmul(temp_s,enchanted_denominator,Z_Y);


  //fcopy(s,temp_s);
  fabsolute(temp_s,temp_s);
  
  pack25519(bytes_out,temp_s);


  return 0;
}


// Inspired by ristretto draft
// hash_to_group or element derivation function takes hash input and turn it into 
// valid ristretto point. In this implementation, input is hexa-string (see more below)
// hash_to_group function consists of 3 steps:
// 1) divide input into 2 halves and mask both hlaves
// 2) MAP both halves so u get 2 points represented with X,Y,Z,T coords (Extended edward's coords)
// 3) perform addition of 2 edward's point, note that we need to add 2 edwards points so fe25519 arithmetics won't fit there
// we need to use function that adds 2 edwards points


// MEMORY ALLOCATION: 3x field_elem
// TOTAL: 14x field_elem
int hash_to_group(u8 bytes_out[BYTES_ELEM_SIZE], const u8 bytes_in[HASH_BYTES_SIZE]){
  ristretto255_point a_point;
  ristretto255_point *a = &a_point;

  ristretto255_point b_point;
  ristretto255_point *b = &b_point;

  ristretto255_point r_point;
  ristretto255_point *r = &r_point;
  

  // make halves
  u8 t1[32], t2[32];
  b_copy(t1,bytes_in);
  b_copy(t2,bytes_in+32);

  
  // MASK LSB for each half, this is equivalent to modulo 2**255
  t1[31] &= 0x7F;
  t2[31] &= 0x7F;

  // encode t1,t2 to field_elem
  field_elem ft1,ft2;
  unpack25519(ft1,t1);
  unpack25519(ft2,t2);
  MAP(a,ft1); // map(ristretto_elligator) first half
  MAP(b,ft2); // map(ristretto_elligator) second hal

  ge25519_addition(r,a,b); // addition of 2 Edward's point

  ristretto255_encode(bytes_out, r);

  return 0;
}

// ristretto255_scalarmult() => scalar multiplication ristretto_point q * scalar s
// Note that scalar "s" is represented as u8[32]
// Inspired by tweetNaCl: https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L632
void ristretto255_scalarmult(ristretto255_point* p, ristretto255_point* q,const u8 *s){
  fcopy(p->x,F_ZERO);
  fcopy(p->y,F_ONE);
  fcopy(p->z,F_ONE);
  fcopy(p->t,F_ZERO);
  for (int i = 255;i >= 0;--i) {
    u8 b = (s[i/8]>>(i&7))&1;
    cswap(p,q,b);
    ge25519_addition(q,q,p);
    ge25519_addition(p,p,p);
    cswap(p,q,b);
  }
}


//////////////////
// MONOCYPHER - needed for efficient calculation of a mod L
/////////////////
/// Utilities ///
/////////////////
static void store32_le(u8 out[4], u32 in){
    out[0] =  in        & 0xff;
    out[1] = (in >>  8) & 0xff;
    out[2] = (in >> 16) & 0xff;
    out[3] = (in >> 24) & 0xff;
}

static void store32_le_buf(u8 *dst, const u32 *src, size_t size) {
    size_t i;
    FOR(i, 0, size) { store32_le(dst + i*4, src[i]); }
}

// get bit from scalar at position i
static int scalar_bit(const u8 s[BYTES_ELEM_SIZE], int i){
    if (i < 0) { return 0; } // handle -1 for sliding windows
    return (s[i>>3] >> (i&7)) & 1;
}


///////////////////////////
/// Arithmetic modulo L ///
///////////////////////////
static const u32 L[8] = {
    0x5cf5d3ed, 0x5812631a, 0xa2f79cd6, 0x14def9de,
    0x00000000, 0x00000000, 0x00000000, 0x10000000,
};

//  p = a*b + p
static void multiply(u32 p[16], const u32 a[8], const u32 b[8]){
   size_t i, j;
    FOR (i, 0, 8) {
        u64 carry = 0;
        FOR (j, 0, 8) {
            carry  += p[i+j] + (u64)a[i] * b[j];
            p[i+j]  = (u32)carry;
            carry >>= 32;
        }
        p[i+8] = (u32)carry;
    }
}

static int is_above_l(const u32 x[8]){
   size_t i;
    // We work with L directly, in a 2's complement encoding
    // (-L == ~L + 1)
    u64 carry = 1;
    FOR (i, 0, 8) {
        carry  += (u64)x[i] + (~L[i] & 0xffffffff);
        carry >>= 32;
    }
    return (int)carry; // carry is either 0 or 1
}

// Final reduction modulo L, by conditionally removing L.
// if x < l     , then r = x
// if l <= x 2*l, then r = x-l
// otherwise the result will be wrong
static void remove_l(u32 r[8], const u32 x[8]){
   size_t i;
    u64 carry = (u64)is_above_l(x);
    u32 mask  = ~(u32)carry + 1; // carry == 0 or 1
    FOR (i, 0, 8) {
        carry += (u64)x[i] + (~L[i] & mask);
        r[i]   = (u32)carry;
        carry >>= 32;
    }
}

// Full reduction modulo L (Barrett reduction)
static void mod_l(u8 reduced[32], const u32 x[16]){
   size_t i, j;
    static const u32 r[9] = {
        0x0a2c131b,0xed9ce5a3,0x086329a7,0x2106215d,
        0xffffffeb,0xffffffff,0xffffffff,0xffffffff,0xf,
    };
    // xr = x * r
    u32 xr[25] = {0};
    FOR (i, 0, 9) {
        u64 carry = 0;
        FOR (j, 0, 16) {
            carry  += xr[i+j] + (u64)r[i] * x[j];
            xr[i+j] = (u32)carry;
            carry >>= 32;
        }
        xr[i+16] = (u32)carry;
    }
    // xr = floor(xr / 2^512) * L
    // Since the result is guaranteed to be below 2*L,
    // it is enough to only compute the first 256 bits.
    // The division is performed by saying xr[i+16]. (16 * 32 = 512)
    ZERO(i, xr, 8);
    FOR (i, 0, 8) {
        u64 carry = 0;
        FOR (j, 0, 8-i) {
            carry   += xr[i+j] + (u64)xr[i+16] * L[j];
            xr[i+j] = (u32)carry;
            carry >>= 32;
        }
    }
    // xr = x - xr
    u64 carry = 1;
    FOR (i, 0, 8) {
        carry  += (u64)x[i] + (~xr[i] & 0xffffffff);
        xr[i]   = (u32)carry;
        carry >>= 32;
    }
    // Final reduction modulo L (conditional subtraction)
    remove_l(xr, xr);
    store32_le_buf(reduced, xr, 8);

//  WIPE_BUFFER(xr);
}

/********************* MD *******************************************************************/
// !!! fnkcne len na CPU s Litle Endian architekturou!!! Pre Big Endian treba doplnit kopirovanie dat!!!
static void multiply_mod_l(u32 r[8], const u32 a[8], const u32 b[8]){
    u32 c[16] = {0};
        multiply(c, a, b);
        mod_l((u8*) r, c);
}


// kod s vyuzitim TMP priestoru na zasobniku 
void inverse_mod_l(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]){
    static const  u8 Lm2[BYTES_ELEM_SIZE] = { // L - 2
        0xeb, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10,
    };
    u32 m_inv [8] = {1};

    // Compute the inverse
    for (int i = 252; i >= 0; i--) {
        multiply_mod_l(m_inv, m_inv, m_inv);

        if (scalar_bit(Lm2, i)) {
            multiply_mod_l(m_inv, m_inv, (u32*) in);

        }
    }
    int i;
    COPY(i ,out, (u8*) m_inv, BYTES_ELEM_SIZE);
}

































