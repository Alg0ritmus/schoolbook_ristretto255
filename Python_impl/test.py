from convertLib import * # all conversions needed are implemented here
from constants import * # constants from draft
import vectors

from libsodium_add import ge25519_p3_add

# This is an school_book/prototype implementation of ristretto255 writen purely in Python.
# Warning: this implementation is not suitable for production by any means.
# It is not secure nor fast, but can it serves as learning tool for better understanding of ristretto255.

# In this implementation we use python int for internal calculations. 
# We also implemented conversion routines to convert inteegers into bytes and bernstein's field_elem from tweetNaCl -> (typedef long long i64; typedef i64 field_elem[16];)
# Theese converion routines (hexToNum/numToHex) can be used to convert between int/bytes/field_elem


# Helpful links:
# TweetNaCl in datails by Martin Kleppmann : https://martin.kleppmann.com/papers/curve25519.pdf
# Ristretto255 Draft: https://datatracker.ietf.org/doc/draft-irtf-cfrg-ristretto255-decaf448/
# Libsodium edward's point addition: https://github.com/jedisct1/libsodium/blob/b7aebe5a1ef46bbb1345e8570fd2e8cea64e587f/src/libsodium/crypto_core/ed25519/ref10/ed25519_ref10.c#L2965
# Ristretto255-Dona | (u*v^7) ^ {(p-5)/8}: https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h


####################################################################
####################################################################
####################################################################



# curve25519 == fe25519 ARITHMETIC, ristretto255 uses this internally
# that means we are using arithmetic modulo P = 2^255 -19
# fe25519 arithmetic is inspired by Bernstein's TweetNaCl, hence hexToNum/numToHex 
# conversion were needed to implement.

# negate a, assuming that a is from <0, P-1>
def fneg(a):
	return (P-a)

# return absolute value of x
def fabs(x):
	return fneg(x) if is_neg(x) else x
	
# true if a is negative, flase otherwise	
def is_neg(a):
	arr_a = numToHex(a,NUMBER_INTERPRETATION_CHOICES["32x8"],False)
	return arr_a[0] & 1

# a*b % (2^255-19)
def fmul(a,b):
	out = (a * b) % P
	return out

# a^n % (2^255-19)
def pow_xtimes(a,n):
	out = a
	for i in range(n):
		out=fmul(out,out)
	return out



# curve25519 == fe25519 SPECIAL FUNCTIONS


# * In:  b =   2^5 - 2^0
# * Out: b = 2^250 - 2^0
# this function is needed to calc z^(2^252 - 3) which is equivalent to z^((p-5)/8)
# more specific, output of this function is z^(2^250 - 2^0)
# note that input of function is NOT "z", but rather z^(2^5 - 2^0)
# in other words, input z^31 == z^(2^5 - 2^0) 
# inspired by:
# https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h
def curve25519_pow_two5mtwo0_two250mtwo0(b):
	t0 = 0
	c = 0

	t0 = pow_xtimes(b,5)
	b = fmul(t0,b)
	t0 = pow_xtimes(b,10)
	c = fmul(t0,b)
	t0 = pow_xtimes(c,20)
	t0 = fmul(t0,c)
	t0 = pow_xtimes(t0,10)
	b = fmul(t0,b)
	t0 = pow_xtimes(b,50)
	c = fmul(t0,b)
	t0 = pow_xtimes(c,100)
	t0 = fmul(t0,c)
	t0 = pow_xtimes(t0,50)
	b = fmul(t0,b)
	return b


# this function calc: 
# z^((p-5)/8) = z^(2^252 - 3)
# needed for ristretto255 inv_sqrt
def curve25519_pow_two252m3(z):
	b=0
	c=0
	t0=0
	two252m3 = 0

	c = pow_xtimes(z,1)
	t0 = pow_xtimes(c,2)
	b = fmul(t0,z)
	c = fmul(b,c)
	t0 = pow_xtimes(c,1)
	b = fmul(t0,b)
	b=curve25519_pow_two5mtwo0_two250mtwo0(b)
	b = pow_xtimes(b,2)
	two252m3 = fmul(b,z)
	return two252m3



# our representation of SQRT_RATIO_M1 from draft
def inv_sqrt(u,v):
	# v^3
	v3 = (v**3) % P
	# v^7
	v7 = (v**7) % P

	# (u*v^3)
	st_bracket = (u*v3) % P 

	# (u*v^7)

	nd_bracket = (u*v7) % P 

	# r = (u*v^7) ^ {(p-5)/8}
	r = curve25519_pow_two252m3(nd_bracket)

	# r2 = [(u*v^7) ^ {(p-5)/8}] * (u*v^3)
	r2 = (r * st_bracket) % P

	temp2 = (r2**2) % P
	check = (v*temp2) % P

	u_neg = fneg(u)
	u_neg_i = (u_neg*SQRT_M1) % P


	correct_sign_sqrt = (check == u);
	flipped_sign_sqrt = (check == u_neg);
	flipped_sign_sqrt_i = (check == u_neg_i);


	r_prime = (r2*SQRT_M1) % P

	should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i

	if should_rotate:
		r2 = r_prime

	r_is_negative = is_neg(r2)
	r_negative = fneg(r2)

	if r_is_negative:
		r2 = r_negative

	was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt
	return (was_nonzero_square,r2)


# Input: inteeger
# Output: ristretto point represented as point with X,Y,Z,T (extended edward's coords)
def ristretto255_decode(s):
	ss = (s*s) % P
	u1 = (1 - ss) % P
	u2 = (1 + ss) % P

	uu1 = (u1*u1) % P
	uu2 = (u2*u2) % P

	duu1 = (EDWARDS_D*uu1) % P
	duu1_neg = fneg(duu1)
	v = (duu1_neg - uu2) % P 
	vuu2 = (v*uu2) % P

	was_square,I = inv_sqrt(1,vuu2)

	Dx = (I*u2) % P
	Dxv = (Dx*v) % P
	Dy = (I*Dxv) % P

	
	sDx = (s*Dx) % P
	x = (sDx+sDx) % P


	if is_neg(x):
		abs_x = fneg(x)
	else:
		abs_x = x

	y = (u1*Dy) % P
	t = (abs_x*y) % P

	if was_square == False:
		print("Decoding fails")
		MSG(f'was_square = {was_square}')
		#raise ValueError

	if is_neg(t):
		print("Decoding fails")
		MSG(f't = {is_neg(t)}')
		#raise ValueError

	if y==0:
		print("Decoding fails")
		MSG(f'y = {y}')
		#raise ValueError


	return (abs_x,y,1,t)



# Input: ristretto point represented as point with X,Y,Z,T (extended edward's coords)
# Output: inteeger, which can be then converted into bytes (check numToHex() in convertLib.py)
def ristretto255_encode(X,Y,Z,T):

	u1 = ( ((Z+Y)% P ) * ((Z-Y)%P) ) % P
	u2 = (X*Y) % P
	u2_2 = (u2*u2) % P
	u1u2_2 = (u1*u2_2) % P
	_,I = inv_sqrt(1,u1u2_2) 

	D1 = (u1*I) % P
	D2 = (u2*I) % P
	Zinv = (((D1*D2) %P) * T) % P

	iX = (X * SQRT_M1) % P
	iY = (Y * SQRT_M1) % P

	enchanted_denominator = (D1 * INVSQRT_A_MINUS_D) % P

	tZinv = (T*Zinv) % P


	if is_neg(tZinv):
		X = iY
		Y = iX
		D_inv = enchanted_denominator
	else:
		X = X 
		Y = Y
		D_inv = D2

	Z = Z

	XZ_inv = (X * Zinv) % P
	if is_neg(XZ_inv):
		Y = fneg(Y)
	else:
		Y = Y

	temp_s = (D_inv * ((Z-Y)%P) ) % P

	if is_neg(temp_s):
		s = fneg(temp_s)
	else:
		s = temp_s
	

	return s 


# MAP function from draft
# MAP is used in hash_to_group() fuction to get ristretto point from hash
# Input parameter is inteeger "t"
# Output is point on Edward's curve with coords X,Y,Z,T, which in our internal representation is
# touple(int,int,int,int)
def MAP(t): # also known as ristretto255_elligator

	r = (t**2) % P
	r =  (SQRT_M1 * r) % P
	u =  (r +1) % P
	u = (u * ONE_MINUS_D_SQ) % P
	c = fneg(1)


	rpd = (r+EDWARDS_D) % P
	v = (r*EDWARDS_D) % P 
	v = (c - v) % P
	v = (v* rpd) % P
	

	was_square, s = inv_sqrt(u,v)

	wasnt_square = 1 - was_square
	
	s_prime = (s * t) % P
	s_prime = fneg(fabs(s_prime))

	if wasnt_square:
		s = s_prime
		c = r

	
	n = (r - 1) % P
	n = (n * c) % P
	n = (n * D_MINUS_ONE_SQ) % P
	n = (n - v) % P
	



	w0 = (s + s) % P
	w0 = (w0 * v) % P
	w1 =  (n * SQRT_AD_MINUS_ONE ) % P
	ss = (s * s) % P
	w2 = (1 - ss ) % P
	w3 = (1 + ss) % P

	return (
		((w0*w3) % P),
		((w2*w1) % P),
		((w1*w3) % P),
		((w0*w2) % P) 
		)



# hash_to_group or element derivation function takes hash input and turn it into 
# valid ristretto point. In this implementation, input is hexa-string (see more below)
# hash_to_group function consists of 3 steps:
# 1) divide input into 2 halves and mask both hlaves
# 2) MAP both halves so u get 2 points represented with X,Y,Z,T coords (Extended edward's coords)
# 3) perform addition of 2 edward's point, note that we need to add 2 edwards points so fe25519 arithmetics won't fit there
# we need to use function that adds 2 edwards points
def hash_to_group(input):
	# divide input into 2 halves, note that "input" paarameter needs to be in hexa-string
	# format, e.g. "5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6"
	# for more examples please check "INPUT_VECTORS_HASH_TO_GROUP_STATIC_HEXSTRING" in vectors.py
	
	t1,t2 = hash_to_num(input) # hash_to_num defined in convertLib.py -> returns 2 inteegers (first and second half)
	t1,t2 = t1% (2**255), t2%(2**255) # masking LSB: according to draft (4.3.4. Element derivation) we need to mask LSB which is equivalent to modulo 2**255

	
	a = MAP(t1) # map(ristretto_elligator) first half
	b = MAP(t2) # map(ristretto_elligator) second half 

	# a + b, this is addition of 2 points on edward's curve 
	# inspired by C cryptographic library "libsodium"
	# note that output "r" has coords X,Y,Z,T
	r = ge25519_p3_add(a,b) # add from libsodium -> libsodium_add.py

	R = ristretto255_encode(*r) # encode "r"" so you get bytes from point
	
	# convert and print result, note that "R" is inteeger in our internal representation
	# so convertion is needed to get bytes [32x8] or [16x16]
	# for more details about convertion please check convertionLib.py 
	numToHex(R,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	numToHex(R,NUMBER_INTERPRETATION_CHOICES["16x16"],True)

	return r


HASH_VECTOR = vectors.INPUT_VECTORS_HASH_TO_GROUP_STATIC_HEXSTRING[0] # vector from vectors.py	
hash_to_group(HASH_VECTOR)

