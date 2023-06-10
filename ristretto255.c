#define DEBUGGING  1
#include "ristretto255.h"
#include "constants.h"

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
		printf("%02hx ", o[i]);
		
	}
	printf("\n");
}

/*-----------------------------------*/

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
	u8 a_32[32],b_32[32];
	pack25519(a_32,a);
	pack25519(b_32,b);

	// constant time
	result &= bytes_eq_32(a_32,b_32);

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

void b_copy(u8 out[32], const u8 in[32]){
	for (int i = 0; i < 32; ++i)
	{
		out[i]  = in[i];
	}
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

void fabsolute(field_elem out, const field_elem in){
	if (is_neg(in)){
		fneg(out,in);
	}
	else{
		fcopy(out,in);
	}
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

// function that calc ristretto255 inverse square root
// https://ristretto.group/formulas/invsqrt.html


// output: (u*v^3) * (u*v^7) ^ {(p-5)/8}
// (p-5)/8 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD => 7237005577332262213973186563042994240829374041602535252466099000494570602493
// but also:
// z^((p-5)/8) = z^(2^252 - 3)
// code inspired by:
// https://github.com/isislovecruft/ristretto-donna/blob/master/src/ristretto-donna.c#L97
int inv_sqrt(field_elem out,const field_elem u, const field_elem v){
	field_elem temp2, v3, v7, st_bracket, nd_bracket, r,r2, check, u_neg,u_neg_i, r_prime,r_negative;
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

	

	fneg(u_neg,u);
	fmul(u_neg_i,u_neg,SQRT_M1);

	correct_sign_sqrt = feq(check, u);
  	flipped_sign_sqrt = feq(check, u_neg);
  	flipped_sign_sqrt_i = feq(check, u_neg_i);


  	

  	fmul(r_prime, r2, SQRT_M1);
	should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i;
	swap25519(r2, r_prime, should_rotate);

	
	// Choose the non-negative square root


	r_is_negative = is_neg(r2);
	fneg(r_negative, r2);

	swap25519(r2, r_negative, r_is_negative);

	was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt;


	// output
	fcopy(out,r2);
	return was_nonzero_square;

}



// generate ristretto point from bytes[32]
int ristretto255_decode(ristretto255_point *ristretto_out, const u8 bytes_in[32]){
	
	int was_square;
	// temp vars
	field_elem s, ss, u1, u2, uu1, uu2,duu1_positive, duu1, v, vuu2, I, Dx, Dxv,Dy, sDx;

	// cords
	field_elem x, abs_x ,y,t;

	int is_canonical, is_negative;

	u8 checked_bytes[32];

	// Step 1: Check that the encoding of the field element is canonical
	unpack25519(s, bytes_in);
	pack25519(checked_bytes,s);
	

	// check if bytes_in == checked_bytes, else abort
	is_canonical = bytes_eq_32(checked_bytes,bytes_in);
	is_negative = is_neg_bytes(bytes_in);

	if (is_canonical == 0 || is_negative==1){
		printf("ristretto255_decode: Bad encoding or neg bytes passed to the ristretto255_decode function! is_canonical=%d, is_negative=%d\n",is_canonical,is_negative);
		return 0;
	}

	// Step 2 calc ristretto255/ge25519 point
	// a = ± 1
	pow2(ss,s);
	fsub(u1,F_ONE,ss); // 1 + as^2
	fadd(u2,F_ONE,ss); // 1 - as^2

	

	pow2(uu1,u1);
	pow2(uu2,u2);


	// d -> from ristretto darft

	fmul(duu1_positive,EDWARDS_D,uu1); // d*u1^2
	fneg(duu1,duu1_positive); // -(D * u1^2) 
	fsub(v,duu1,uu2); // adu1^2 - u2^2

	fmul(vuu2,v,uu2); // v*u2^2

	
	was_square = inv_sqrt(I,F_ONE,vuu2);

	fmul(Dx,I,u2); // I*u2
	fmul(Dxv, Dx, v); // Dx * v
	fmul(Dy, I, Dxv); // I*Dx*v

	

	// cords
	// x: |2sDx|
	fmul(sDx,s,Dx);
	fadd(x,sDx,sDx); // 2sDx

	if (is_neg(x)){
		fneg(abs_x,x); // |2sDx|

	}
	else{
		fcopy(abs_x,x);
	}

	// y: u1Dy
	fmul(y,u1,Dy);

	// t: x,y
	fmul(t,abs_x,y);

	fcopy(ristretto_out->x, abs_x);
	fcopy(ristretto_out->y, y);
	fcopy(ristretto_out->z, F_ONE);
	fcopy(ristretto_out->t, t); 

	if (was_square == 0){
		printf("\n\n\n ristretto255_decode: Bad encoding! was_square=%d \n\n\n",was_square);
		return 0;
	}
	if (is_neg(t)){
		printf("\n\n\n ristretto255_decode: Bad encoding! is_neg(t)=%d \n\n\n",is_neg(t));
		return 0;
	}
	if (feq(y,F_ZERO)){
		printf("\n\n\n ristretto255_decode: Bad encoding! feq(y,F_ZERO)=%d \n\n\n",feq(y,F_ZERO));
		return 0;
	}
	return 1;
}



int ristretto255_encode(unsigned char bytes_out[32], const ristretto255_point* ristretto_in){

	field_elem temp_zy1,temp_zy2,u1,u2,uu2,u1uu2,I,D1,D2, D1D2,Zinv, iX, iY, enchanted_denominator, tZinv;
	field_elem _X, _Y,_Z, D_inv,XZ_inv,temp_s,s,Z_Y;

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
	

	fmul(iX, ristretto_in->x, SQRT_M1);
	fmul(iY, ristretto_in->y, SQRT_M1);

	fmul(enchanted_denominator,D1,INVSQRT_A_MINUS_D);
	fmul(tZinv,ristretto_in->t,Zinv);


	if (is_neg(tZinv)){
		fcopy(_X,iY);
		fcopy(_Y,iX);
		fcopy(D_inv,enchanted_denominator);
		
	}
	else{
		fcopy(_X,ristretto_in->x);
		fcopy(_Y,ristretto_in->y);
		fcopy(D_inv,D2);
	}
	fcopy(_Z,ristretto_in->z);

	fmul(XZ_inv,_X,Zinv);

	if (is_neg(XZ_inv)){
		fneg(_Y,_Y);
	}
	else{
		fcopy(_Y,_Y);
	}


	fsub(Z_Y,_Z,_Y);

	fmul(temp_s,D_inv,Z_Y);

	if(is_neg(temp_s)){
		fneg(s,temp_s);
	}
	else{
		fcopy(s,temp_s);
	}

	pack25519(bytes_out,s);


	return 1;
}





int MAP(ristretto255_point* ristretto_out, const field_elem t){
	field_elem r,u,c,rpd,v,s,s_prime,n, w0,w1,w2,w3,ss,x1,y1,z1,t1;
	int was_square, wasnt_square;

	fmul(r,t,t);
	fmul(r,SQRT_M1,r);
	fadd(u,r,F_ONE);
	fmul(u,u,ONE_MINUS_D_SQ);
	fneg(c,F_ONE);

	fadd(rpd,r,EDWARDS_D);
	fmul(v,r,EDWARDS_D);
	fsub(v,c,v);
	fmul(v,v,rpd);

	was_square = inv_sqrt(s,u,v);
	wasnt_square = 1 - was_square;

	fmul(s_prime,s,t);
	fabsolute(s_prime,s_prime);
	fneg(s_prime,s_prime);

	if (wasnt_square){
		fcopy(s,s_prime);
		fcopy(c,r);
	}

	fsub(n,r,F_ONE);
	fmul(n,n,c);
	fmul(n,n,D_MINUS_ONE_SQ);
	fsub(n,n,v);

	fadd(w0,s,s);
	fmul(w0,w0,v);
	fmul(w1,n,SQRT_AD_MINUS_ONE);
	fmul(ss,s,s);
	fsub(w2,F_ONE,ss);
	fadd(w3,F_ONE,ss);

	

	fmul(x1,w0,w3);
	fmul(y1,w2,w1);
	fmul(z1,w1,w3);
	fmul(t1,w0,w2);

	fcopy(ristretto_out->x, x1);
	fcopy(ristretto_out->y, y1);
	fcopy(ristretto_out->z, z1);
	fcopy(ristretto_out->t, t1); 


	return 1;
}


// LIBSODIUM ADD
void ge25519_p3_add(ristretto255_point* r,const ristretto255_point* p,const ristretto255_point* q){
	ristretto255_point q_cached_point;
	ristretto255_point *q_cached = &q_cached_point;

	ristretto255_point p1p1_point;
	ristretto255_point *p1p1 = &p1p1_point;

	ge25519_p3_to_cached(q_cached,q);
	ge25519_add_cached(p1p1,p,q_cached);
	ge25519_p1p1_to_p3(r,p1p1);
}

void ge25519_p3_to_cached(ristretto255_point* p_out,const ristretto255_point* p){


	fadd(p_out->x, p->y,p->x);
	fsub(p_out->y, p->y,p->x);
	fcopy(p_out->z, p->z);
	fmul(p_out->t, p->t,EDWARDS_D2); 
}


void ge25519_add_cached(ristretto255_point* r,const ristretto255_point* p,const ristretto255_point* q){
	field_elem t;

	fadd(r->x,p->y,p->x);
	fsub(r->y,p->y,p->x);

	fmul(r->z,r->x,q->x);
	fmul(r->y,r->y,q->y);
	fmul(r->t,q->t,p->t);
	fmul(r->x,p->z,q->z);

	fadd(t,r->x,r->x);
	fsub(r->x,r->z,r->y);
	fadd(r->y,r->z,r->y);
	fadd(r->z,t,r->t);
	fsub(r->t,t,r->t);
}
	
void ge25519_p1p1_to_p3(ristretto255_point* r,const ristretto255_point* p){
	fmul(r->x,p->x,p->t);
	fmul(r->y,p->y,p->z);
	fmul(r->z,p->z,p->t);
	fmul(r->t,p->x,p->y);

}
///


int hash_to_group(u8 bytes_out[32], const u8 bytes_in[64]){
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

	
	// MASK LSB for each half
	t1[31] &= 0x7F;
	t2[31] &= 0x7F;

	// encode t1,t2 to field_elem
	field_elem ft1,ft2;
	unpack25519(ft1,t1);
	unpack25519(ft2,t2);
	MAP(a,ft1); // map(ristretto_elligator) first half
	MAP(b,ft2); // map(ristretto_elligator) second hal

	ge25519_p3_add(r,a,b);

	ristretto255_encode(bytes_out, r);

	return 1;
}






























