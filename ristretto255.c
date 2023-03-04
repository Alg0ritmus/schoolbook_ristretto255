#include "ristretto255.h"

/*		UTILS		*/
void print(field_elem o){

//	printf("PRINT: ");
	for (int i=0;i<16;i++){
		printf("%02hhx ", o[i]);
		
	}
	printf("\n");
}

void print_32(const u8* o){

//	printf("PRINT: ");
	for (int i=0;i<32;i++){
		printf("%02hhx ", o[i]);
		
	}
	printf("\n");
}

/*-----------------------------------*/


// little-endian order --> a = a0*2^0 + a1*2^16 + a2*2^32 + ... + a15*2^240 (vzdy o ax*2^ o15 vyssie)
const field_elem F_ZERO = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
const field_elem F_ONE = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

const field_elem F_BIGGEST = {	32767,32767,32767,32767,
								32767,32767,32767,32767,
								32767,32767,32767,32767,
								32767,32767,32767,32767
};

const field_elem _121665 = {0xDB41, 1};



void unpack25519(field_elem out, const u8 *in)
{
int i;
for (i = 0; i < 16; ++i) out[i] = in[2*i] + ((i64) in[2*i + 1] << 8);
out[15] &= 0x7fff;
 }

 void carry25519(field_elem elem)
 {
 int i;
 i64 carry;
 for (i = 0; i < 16; ++i) {
 carry = elem[i] >> 16;
 elem[i] -= carry << 16;
 if (i < 15) elem[i + 1] += carry; else elem[0] += 38 * carry;
 }
 }

 void fadd(field_elem out, const field_elem a, const field_elem b) /* out = a + b */
 {
 int i;
 for (i = 0; i < 16; ++i) out[i] = a[i] + b[i];
 }

 void fsub(field_elem out, const field_elem a, const field_elem b) /* out = a - b */
 {
 int i;
 for (i = 0; i < 16; ++i) out[i] = a[i] - b[i];
 }

 void fmul(field_elem out, const field_elem a, const field_elem b) /* out = a * b */
 {
 i64 i, j, product[31];
 for (i = 0; i < 31; ++i) product[i] = 0;
 for (i = 0; i < 16; ++i) {
 for (j = 0; j < 16; ++j) product[i+j] += a[i] * b[j];
 }
 for (i = 0; i < 15; ++i) product[i] += 38 * product[i+16];
 for (i = 0; i < 16; ++i) out[i] = product[i];
 carry25519(out);
 carry25519(out);
 }

void finverse(field_elem out, const field_elem in)
{
field_elem c;
int i;
for (i = 0; i < 16; ++i) c[i] = in[i];
 for (i = 253; i >= 0; i--) {
 fmul(c,c,c);
 if (i != 2 && i != 4) fmul(c,c,in);
 }
 for (i = 0; i < 16; ++i) out[i] = c[i];
 }

 void swap25519(field_elem p, field_elem q, int bit)
 {
 i64 t, i, c = ~(bit - 1);
 for (i = 0; i < 16; ++i) {
 t = c & (p[i] ^ q[i]);
 p[i] ^= t;
 q[i] ^= t;
 }
 }

 void pack25519(u8 *out, const field_elem in)
 {
 int i, j, carry;
 field_elem m, t;
 for (i = 0; i < 16; ++i) t[i] = in[i];
 carry25519(t); carry25519(t); carry25519(t);
 for (j = 0; j < 2; ++j) {
 m[0] = t[0] - 0xffed;
 for(i = 1; i < 15; i++) {
 m[i] = t[i] - 0xffff - ((m[i-1] >> 16) & 1);
 m[i-1] &= 0xffff;
 }
 m[15] = t[15] - 0x7fff - ((m[14] >> 16) & 1);
 carry = (m[15] >> 16) & 1;
 m[14] &= 0xffff;
 swap25519(t, m, 1 - carry);
 }
 for (i = 0; i < 16; ++i) {
 out[2*i] = t[i] & 0xff;
 out[2*i + 1] = t[i] >> 8;
 }
 }



void scalarmult(u8 *out, const u8 *scalar, const u8 *point)
{
u8 clamped[32];
i64 bit, i;
field_elem a, b, c, d, e, f, x;
 for (i = 0; i < 31; ++i) clamped[i] = scalar[i];
 clamped[0] &= 0xf8;
 clamped[31] = (clamped[31] & 0x7f) | 0x40;
 unpack25519(x, point);
 for (i = 0; i < 16; ++i) {
 b[i] = x[i];
 d[i] = a[i] = c[i] = 0;
 }
 a[0] = d[0] = 1;
 for (i = 254; i >= 0; --i) {
 bit = (clamped[i >> 3] >> (i & 7)) & 1;
 swap25519(a, b, bit);
 swap25519(c, d, bit);
 fadd(e, a, c);
 fsub(a, a, c);
 fadd(c, b, d);
 fsub(b, b, d);
 fmul(d, e, e);
 fmul(f, a, a);
 fmul(a, c, a);
 fmul(c, b, e);
 fadd(e, a, c);
 fsub(a, a, c);
 fmul(b, a, a);
 fsub(c, d, f);
 fmul(a, c, _121665);
 fadd(a, a, d);
 fmul(c, c, a);
 fmul(a, d, f);
 fmul(d, b, x);
 fmul(b, e, e);
 swap25519(a, b, bit);
 swap25519(c, d, bit);
 }
 finverse(c, c);
 fmul(a, a, c);
 pack25519(out, a);
 }

/* ------------------------------------------------*/
/* -------------------my code----------------------*/
/* ------------------------------------------------*/

// return 1 if they're equal
int feq( const field_elem a,  const field_elem b){
	int result = 1;

	// constant time
	result &= a[0] == b[0];
	result &= a[1] == b[1];
	result &= a[2] == b[2];
	result &= a[3] == b[3];

	result &= a[4] == b[4];
	result &= a[5] == b[5];
	result &= a[6] == b[6];
	result &= a[7] == b[7];

	result &= a[8] == b[8];
	result &= a[9] == b[9];
	result &= a[10] == b[10];
	result &= a[11] == b[11];

	result &= a[12] == b[12];
	result &= a[13] == b[13];
	result &= a[14] == b[14];
	result &= a[15] == b[15];


	return result;
}

// copy in to out
void fcopy(field_elem out, const field_elem in){
	out[0]  = in[0];
	out[1]  = in[1];
	out[2]  = in[2];
	out[3]  = in[3];
	out[4]  = in[4];
	out[5]  = in[5];
	out[6]  = in[6];
	out[7]  = in[7];
	out[8]  = in[8];
	out[9]  = in[9];
	out[10] = in[10];
	out[11] = in[11];
	out[12] = in[12];
	out[13] = in[13];
	out[14] = in[14];
	out[15] = in[15];

}

// negation of element a -> 0-a = -a
void fneg(field_elem out, const field_elem in){
	fsub(out, F_ZERO, in);
};

 // square a˘2

 void pow2(field_elem out, const field_elem a){
 	fmul(out,a,a);
 }

 
// square a˘3
void pow3(field_elem out, const field_elem a){
 	field_elem a2;
 	// temporary: a^2
 	pow2(a2,a);

 	// a^3
 	fmul(out,a2,a);
}


void pow7(field_elem out, const field_elem a){
 	field_elem a3, a6;
 	pow3(a3,a);
 	pow2(a6,a3);
 	fmul(out,a6,a);
}

// not efficient but can calc a^ (2 * n)
void pow_xtimes(field_elem out, const field_elem a, int n){
	fcopy(out,a);
	for (int i = 0; i < n; ++i)
	{
		pow2(out,out);
	}
	

}

// * In:  b =   2^5 - 2^0
// * Out: b = 2^250 - 2^0
// this function is needed to calc z^(2^252 - 3)
// more specific, output of this function is z^(2^250 - 2^0)
// note that input of function is NOT "z", but rather z^(2^5 - 2^0)
// in other words, input z^31 == z^(2^5 - 2^0) 
// inspired by:
// https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h
void curve25519_pow_two5mtwo0_two250mtwo0(field_elem b){
	field_elem t0,c;

	/* 2^5  - 2^0 */ /* b */
	/* 2^10 - 2^5 */ pow_xtimes(t0, b, 5);
	/* 2^10 - 2^0 */ fmul(b, t0, b);
	/* 2^20 - 2^10 */ pow_xtimes(t0, b, 10);
	/* 2^20 - 2^0 */ fmul(c, t0, b);
	/* 2^40 - 2^20 */ pow_xtimes(t0, c, 20);
	/* 2^40 - 2^0 */ fmul(t0, t0, c);
	/* 2^50 - 2^10 */ pow_xtimes(t0, t0, 10);
	/* 2^50 - 2^0 */ fmul(b, t0, b);
	/* 2^100 - 2^50 */ pow_xtimes(t0, b, 50);
	/* 2^100 - 2^0 */ fmul(c, t0, b);
	/* 2^200 - 2^100 */ pow_xtimes(t0, c, 100);
	/* 2^200 - 2^0 */ fmul(t0, t0, c);
	/* 2^250 - 2^50 */ pow_xtimes(t0, t0, 50);
	/* 2^250 - 2^0 */ fmul(b, t0, b);

}

// this function calc: 
// z^((p-5)/8) = z^(2^252 - 3)
// needed for ristretto255 inv_sqrt
void curve25519_pow_two252m3(field_elem two252m3, const field_elem z){
	field_elem b,c,t0;

	/* 2 */ pow_xtimes(c, z, 1); /* c = 2 */
	/* 8 */ pow_xtimes(t0, c, 2); /* t0 = 8 */
	/* 9 */ fmul(b, t0, z); /* b = 9 */
	/* 11 */ fmul(c, b, c); /* c = 11 */
	/* 22 */ pow_xtimes(t0, c, 1);
	/* 2^5 - 2^0 = 31 */ fmul(b, t0, b);
	/* 2^250 - 2^0 */ curve25519_pow_two5mtwo0_two250mtwo0(b);
	/* 2^252 - 2^2 */ pow_xtimes(b, b, 2);
	/* 2^252 - 3 */ fmul(two252m3, b, z);
}



// function that calc ristretto255 inverse square root
// https://ristretto.group/formulas/invsqrt.html


// output: (u*v^3) * (u*v^7) ^ {(p-5)/8}
// (p-5)/8 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD => 7237005577332262213973186563042994240829374041602535252466099000494570602493
// but also:
// z^((p-5)/8) = z^(2^252 - 3)
// code inspired by:
// https://github.com/isislovecruft/ristretto-donna/blob/master/src/ristretto-donna.c#L97
void inv_sqrt(field_elem out,const field_elem u, const field_elem v){
	field_elem temp,temp2, v3, v7, p58, st_bracket, nd_bracket, r,r2, check, u_neg,u_neg_i, r_prime;
	int r_is_negative;
	int correct_sign_sqrt;
	int flipped_sign_sqrt;
	int flipped_sign_sqrt_i;
	int was_nonzero_square;
	int should_rotate;



	pow3(v3,v); 							//v^3
	pow7(v7,v); 							//v^7
	fmul(st_bracket,u,v3);					// (u*v^3)
	fmul(nd_bracket,u,v7);					// (u*v^7)
	curve25519_pow_two252m3(r,nd_bracket); 	// r = (u*v^7) ^ {(p-5)/8}

	fmul(r2,r,st_bracket);					// r2 = (u*v^3) * (u*v^7) ^ {(p-5)/8}


	

	pow2(temp2,r2);							// tmp = r2 ^ 2 -> needed for check
	fmul(check,v,temp2);					// check = (r2 ^ 2) * v

	// Check if everything is correct: check == u ?
	printf("\nPrint check:");
	print(check);

	printf("\nPrint expected_check:");
	print(u);


	fneg(u_neg,u);
	//fmul(u_neg_i,u_neg,SQRT_M1);

	//correct_sign_sqrt = feq(check, u);
  	//flipped_sign_sqrt = feq(check, u_neg);
  	//flipped_sign_sqrt_i = feq(check, u_neg_i);

  	//fmul(r_prime, r2, SQRT_M1);
	//should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i;
	//curve25519_swap_conditional(r2, r_prime, should_rotate);




	// dummy output
	fmul(out,v,v);

	/*

	// additional vars
	field_elem u_neg, u_neg_i;

	// assaign vars
	correct_sign_sqrt = bignum25519_ct_eq(check, u);
  	flipped_sign_sqrt = bignum25519_ct_eq(check, u_neg);
  	flipped_sign_sqrt_i = bignum25519_ct_eq(check, u_neg_i);

	
	// negation
	// fneg(u_neg, u);
	// fmul(u_neg_i, u_neg, SQRT_M1);

	// fmul(r_prime, r, SQRT_M1);
	// conditional swap
	// swap25519(r, r_prime, should_rotate)

	*/
}








































