#define DEBUGGING  1
#include "ristretto255.h"

/*		UTILS		*/
void print(field_elem o){

//	printf("PRINT: ");
	for (int i=0;i<16;i++){
		printf("%llx ", o[i]);
		
	}
	printf("\n");
}


void pack_and_print_32(field_elem o){
	u8 temp[32];
	pack25519(temp,o);
	print_32(temp);
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
const field_elem F_TWO = {2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

const field_elem F_BIGGEST = {	32767,32767,32767,32767,
								32767,32767,32767,32767,
								32767,32767,32767,32767,
								32767,32767,32767,32767
};

const field_elem _121665 = {0xDB41, 1};

// 2^255 - 19 -> from script convert_num_to_bernstein.py
const field_elem F_MODULUS = {
	0xffed, 0xffff, 0xffff, 0xffff,
	0xffff, 0xffff, 0xffff, 0xffff,
	0xffff, 0xffff, 0xffff, 0xffff,
	0xffff, 0xffff, 0xffff, 0x7fff
};


const field_elem SQRT_M1 = {
	0xa0b0, 0x4a0e, 0x1b27, 0xc4ee,
	0xe478, 0xad2f, 0x1806, 0x2f43,
	0xd7a7, 0x3dfb, 0x99, 	0x2b4d,
	0xdf0b, 0x4fc1, 0x2480, 0x2b83
};

const field_elem EDWARDS_D = {
	0x78a3, 0x1359, 0x4dca, 0x75eb,
	0xd8ab, 0x4141, 0xa4d, 0x70,
	0xe898, 0x7779, 0x4079, 0x8cc7,
	0xfe73, 0x2b6f, 0x6cee, 0x5203

};

const field_elem INVSQRT_A_MINUS_D  = {
	0x40ea, 0x805d, 0xfdaa, 0x99c8,
	0x72be, 0x5a41, 0x1617, 0x9d2f,
	0xd840, 0xfe01, 0x7b91, 0x16c2,
	0xfca2, 0xcfaf, 0x8905, 0x786c
};


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

// if bit = 1 -> swap, otherwise do nothing
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


// return 1 if they're equal
int bytes_eq_32( const u8 a[32],  const u8 b[32]){
	int result = 1;

	for (int i = 0; i < 32; ++i){
		result &= a[i] == b[i];
		//printf("a=%02hhx| b=%02hhx \n",a[i],b[i] );
	}

	return result;
}

// copy in to out
void fcopy(field_elem out, const field_elem in){
	for (int i = 0; i < 16; ++i)
	{
		printf("%llx ", in[i]);
	}
	printf(" ---\n");
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

// negation of element a -> p-a = -a 
// inspired by: https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_51.h#L94

// alebo zmenit na :
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_25_5.h#L116
void fneg(field_elem out, const field_elem in){
	fsub(out, F_MODULUS, in);
	carry25519(out);
};

// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/include/sodium/private/ed25519_ref10_fe_25_5.h#L302
int is_neg(const field_elem in){
	u8 temp[32];
	pack25519(temp, in);

	return temp[0] & 1;
}

int is_neg_bytes(const u8 in[32]){

	return in[0] & 1;
}



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


// copied from:
// https://github.com/sbp/tweetnacl-tools/blob/master/tweetnacl.c#L382
#define FOR(i,n) for (i = 0;i < n;++i)
#define sv static void

sv car25519(field_elem o)
{
  int i;
  i64 c;
  FOR(i,16) {
    o[i]+=(1LL<<16);
    c=o[i]>>16;
    o[(i+1)*(i<15)]+=c-1+37*(c-1)*(i==15);
    o[i]-=c<<16;
  }
}

sv M(field_elem o,const field_elem a,const field_elem b)
{
  i64 i,j,t[31];
  FOR(i,31) t[i]=0;
  FOR(i,16) FOR(j,16) t[i+j]+=a[i]*b[j];
  FOR(i,15) t[i]+=38*t[i+16];
  FOR(i,16) o[i]=t[i];
  car25519(o);
  car25519(o);
}

sv S(field_elem o,const field_elem a)
{
  M(o,a,a);
}
sv pow2523(field_elem o,const field_elem i)
{
  field_elem c;
  int a;
  FOR(a,16) c[a]=i[a];
  for(a=250;a>=0;a--) {
    S(c,c);
    if(a!=1) M(c,c,i);
  }
  FOR(a,16) o[a]=c[a];
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
	field_elem temp,temp2, v3, v7, p58, st_bracket, nd_bracket, r,r2, check, u_neg,u_neg_i, r_prime,r_negative;
	u8 r_bytes;
	int r_is_negative;
	int correct_sign_sqrt;
	int flipped_sign_sqrt;
	int flipped_sign_sqrt_i;
	int was_nonzero_square;
	int should_rotate;
	u8 temp_bytes[32];


	pow3(v3,v); 							//v^3
	pow7(v7,v); 							//v^7


	fmul(st_bracket,u,v3);					// (u*v^3)
	fmul(nd_bracket,u,v7);					// (u*v^7)

	
	curve25519_pow_two252m3(r,nd_bracket); 	// r = (u*v^7) ^ {(p-5)/8}
	//pow2523(r,nd_bracket);
	
	fmul(r2,r,st_bracket);					// r2 = (u*v^3) * (u*v^7) ^ {(p-5)/8}


	

	pow2(temp2,r2);							// tmp = r2 ^ 2 -> needed for check
	fmul(check,v,temp2);					// check = (r2 ^ 2) * v

	// Check if everything is correct: check == u ?
	
	printf("\nPrint check IN:");
	print(v);
	
	printf("\nPrint check:");
	print(check);

	//printf("\nPrint expected_check:");
	//print(u);
	

	fneg(u_neg,u);
	fmul(u_neg_i,u_neg,SQRT_M1);

	correct_sign_sqrt = feq(check, u);
  	flipped_sign_sqrt = feq(check, u_neg);
  	flipped_sign_sqrt_i = feq(check, u_neg_i);

  	/*
  	printf("correct_sign_sqrt = %d\n", correct_sign_sqrt);
  	printf("flipped_sign_sqrt = %d\n", flipped_sign_sqrt);
  	printf("flipped_sign_sqrt_i = %d\n", flipped_sign_sqrt_i);
	*/

  	fmul(r_prime, r2, SQRT_M1);
	should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i;
	swap25519(r2, r_prime, should_rotate);

	
	// Choose the non-negative square root


	r_is_negative = is_neg(r2);
	fneg(r_negative, r2);
	/*
	printf("\nPrint check r:\n");
	print(r2);
	print(r_negative);
	printf("r:%d r_negative ma byt pozitivne:%d\n",r_is_negative,is_neg(r_negative));
	*/
	swap25519(r2, r_negative, r_is_negative);


	//printf("\n r2 after swap %d:\n",r_is_negative);
	//print(r2);

	was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt;
	//printf("\nwas_nonzero_square = %d\n", was_nonzero_square);

	// output
	fcopy(out,r2);

}


/*

typedef struct ge_point25519{
	field_elem x,y,z,t;
}ristretto255_point;

*/

// generate ristretto point from bytes[32]
int ristretto255_decode(ristretto255_point *ristretto_out, const u8 bytes_in[32]){
	
	// temp vars
	field_elem s, ss, u1, u2, uu1, uu2,duu1_positive, duu1, v, vuu2, I, Dx, Dxv,Dy, sDx;

	// cords
	field_elem x, abs_x ,y,t;

	int is_canonical, is_negative;

	u8 checked_bytes[32], temporary[32];

	// Step 1: Check that the encoding of the field element is canonical
	unpack25519(s, bytes_in);
	pack25519(checked_bytes,s);
	

	// check if bytes_in == checked_bytes, else abort
	is_canonical = bytes_eq_32(checked_bytes,bytes_in);
	is_negative = is_neg_bytes(bytes_in);
	printf("je negative???: %d\n", bytes_in[0]&1);

	if (is_canonical == 0 || is_negative==1){
		printf("Bad encoding or neg bytes passed to the ristretto255_decode function! is_canonical=%d, is_negative=%d\n",is_canonical,is_negative);
		return 0;
	}

	// Step 2 calc ristretto255/ge25519 point
	// a = ± 1
	pow2(ss,s);
	fadd(u1,F_ONE,ss); // 1 + as^2
	fsub(u2,F_ONE,ss); // 1 - as^2

	

	pow2(uu1,u1);
	pow2(uu2,u2);


	// d -> from ristretto darft

	fmul(duu1_positive,EDWARDS_D,uu1); // d*u1^2
	fneg(duu1,duu1_positive); // -(D * u1^2) 
	fsub(v,duu1,uu2); // adu1^2 - u2^2

	fmul(vuu2,v,uu2); // v*u2^2

	
	inv_sqrt(I,F_ONE,vuu2);

	fmul(Dx,I,u2); // I*u2
	fmul(Dxv, Dx, v); // Dx * v
	fmul(Dy, I, Dxv); // I*Dx*v

	

	// cords
	// x: |2sDx|
	fmul(sDx,s,Dx);
	fadd(x,sDx,sDx); // 2sDx

	/*
	pack25519(temporary,x);
	printf("\nx:\n");
	print_32(temporary);
	*/

	if (is_neg(x)){
		fneg(abs_x,x); // |2sDx|

	}
	else{
		fcopy(abs_x,x);
	}

	// y: u1Dy
	fmul(y,u1,Dy);

	// t: x,y
	fmul(t,x,y);

	fcopy(ristretto_out->x, abs_x);
	fcopy(ristretto_out->y, y);
	fcopy(ristretto_out->z, F_ONE);
	fcopy(ristretto_out->t, t); 
	
	return 1;
}



int ristretto255_encode(unsigned char bytes_out[32], const ristretto255_point* ristretto_in){

	field_elem temp_zy1,temp_zy2,u1,u2,uu2,u1uu2,I,D1,D2,D, D1D2,Zinv, TZinv, ix,iy,xZiv,_y, z_y,z_yD,s;

	printf("in coding&&&&&&&c\n");
	fadd(temp_zy1,ristretto_in->z,ristretto_in->y);
	fsub(temp_zy2,ristretto_in->z,ristretto_in->y);
	fmul(u1,temp_zy1,temp_zy2); // (Z+Y)(Z-Y)

	fmul(u2,ristretto_in->x,ristretto_in->y); // XY

	pow2(uu2,u2); //u2^2 
	fmul(u1uu2,u1,uu2); // u1 * u2^2
	inv_sqrt(I,F_ONE,u1uu2);

	fmul(D1,u1,I); // u1 * I
	fmul(D2,u2,I); // u2 * I


	fmul(D1D2,D1,D2); // D1 * D2
	fmul(Zinv,D1D2,ristretto_in->t); // Zinv = D1*D2*T

	fmul(TZinv,ristretto_in->t,Zinv);

	// mozno chybne??
	if (is_neg(TZinv)){
		fmul(ix,ristretto_in->x,SQRT_M1); // ix = x*sqrt(-1)
		fmul(iy,ristretto_in->y,SQRT_M1); // iy = y*sqrt(-1)
		fmul(D,D1,INVSQRT_A_MINUS_D); // D = D1/√(a-d)  -> 1/√(a-d) ->  INVSQRT_A_MINUS_D

	}
	
	else{
		fcopy(ristretto_in->x,ristretto_in->x);
		fcopy(ristretto_in->y,ristretto_in->y);
		fcopy(ristretto_in->z,ristretto_in->z);
		fcopy(D,D2);


		// y = y
		// D = D2
	}
	

	fmul(xZiv,ristretto_in->x,Zinv);
	if (is_neg(xZiv)){
		fneg(_y,ristretto_in->y); // -y
		fcopy(ristretto_in->y,_y); // y = -y
	}

	fsub(z_y ,ristretto_in->z,ristretto_in->y);
	fmul(z_yD,z_y,D);

	if (is_neg(z_yD)){
		fneg(s,z_yD);
	}

	pack25519(bytes_out,s);
	printf("in coding-------\n");
	print_32(bytes_out);
	printf("in coding-------\n");

	// if X*Zinv in_neg
	// then y = - y



	return 1;
}





































