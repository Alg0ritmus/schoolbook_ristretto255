from convertLib import * # all conversions needed are implemented here
from constants import * # constants from draft
from libsodium_add import ge25519_p3_add
import vectors
# This is an school_book/prototype implementation of ristretto255 writen purely in Python.
# Warning: this implementation is not suitable for production by any means.
# It is not secure nor fast, but can it serves as learning tool for better understanding of ristretto255.

# In this implementation we use python int for internal calculations. 
# We also implemented conversion routines to convert inteegers into bytes and bernsteins field_elem from tweetNaCl -> (typedef long long i64; typedef i64 field_elem[16];)
# Theese converion routines (hexToNum/numToHex) can be used to convert between int/bytes/field_elem


# Helpful links:
# TweetNaCl in datails by Martin Kleppmann : https://martin.kleppmann.com/papers/curve25519.pdf
# Ristretto255 Draft: https://datatracker.ietf.org/doc/draft-irtf-cfrg-ristretto255-decaf448/
# Libsodium edwards point addition: https://github.com/jedisct1/libsodium/blob/b7aebe5a1ef46bbb1345e8570fd2e8cea64e587f/src/libsodium/crypto_core/ed25519/ref10/ed25519_ref10.c#L2965
# Ristretto255-Dona | (u*v^7) ^ {(p-5)/8}: https://github.com/floodyberry/ed25519-donna/blob/master/curve25519-donna-helpers.h


####################################################################
####################################################################
####################################################################



# curve25519 == fe25519 ARITHMETIC, ristretto255 uses this internally
# that means we are using arithmetic modulo P = 2^255 -19
# fe25519 arithmetic is inspired by Bernsteins TweetNaCl, hence hexToNum/numToHex 
# conversion were needed to implement.


def emil_red(a):
	return a % REDUCE


# negate a, assuming that a is from <0, P-1>
def fneg(a):
	if use25519:
		return P-a
	a = a % REDUCE
	return (REDUCE-a)


# return absolute value of x
def fabs(x):
	x = emil_red(x)
	return fneg(x) if is_neg(x) else x
	
# true if a is negative, flase otherwise	
def is_neg(a):
	#if not use25519:
		#a = a % REDUCE
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

	check = emil_red(check)
	u = emil_red(u)
	u_neg_i = emil_red(u_neg_i)

	correct_sign_sqrt = (check == u)
	flipped_sign_sqrt = (check == u_neg)
	flipped_sign_sqrt_i = (check == u_neg_i)


	r_prime = (r2*SQRT_M1) % P

	should_rotate = flipped_sign_sqrt | flipped_sign_sqrt_i

	if should_rotate:
		r2 = r_prime

	r2 = emil_red(r2)
	r_is_negative = is_neg(r2)
	r_negative = fneg(r2)

	if r_is_negative:
		r2 = r_negative

	was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt
	return (was_nonzero_square,r2)


# Input: inteeger
# Output: ristretto point represented as point with X,Y,Z,T (extended edwards coords)
def ristretto255_decode(s):
	#check if cannonical
	a =  numToHex(s)
	aa = hexToNum(a)
	if aa!=s or is_neg(s):
		print("Decoding fails, is can:",aa==s,"neg:",is_neg(s))
		MSG(f'non-canonical input')
		raise ValueError


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

	x = emil_red(x)
	if is_neg(x):
		abs_x = fneg(x)
	else:
		abs_x = x

	y = (u1*Dy) % P
	t = (abs_x*y) % P

	if was_square == False:
		print("Decoding fails")
		MSG(f'was_square = {was_square}')
		raise ValueError

	t = emil_red(t)
	if is_neg(t):
		print("Decoding fails")
		MSG(f't = {is_neg(t)}')
		raise ValueError

	if y==0:
		print("Decoding fails")
		MSG(f'y = {y}')
		raise ValueError

	if not use25519:
		return (abs_x,y%REDUCE,1,t%REDUCE)

	return (abs_x,y,1,t)



# Input: ristretto point represented as point with X,Y,Z,T (extended edwards coords)
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

	tZinv = emil_red(tZinv)
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
	XZ_inv = emil_red(XZ_inv)
	if is_neg(XZ_inv):
		Y = fneg(Y)
	else:
		Y = Y

	temp_s = (D_inv * ((Z-Y)%P) ) % P

	temp_s = emil_red(temp_s)
	if is_neg(temp_s):
		s = fneg(temp_s)
	else:
		s = temp_s
	
	if not use25519:
		return s % REDUCE
	return s


# MAP function from draft
# MAP is used in hash_to_group() fuction to get ristretto point from hash
# Input parameter is inteeger "t"
# Output is point on Edwards curve with coords X,Y,Z,T, which in our internal representation is
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

	if not use25519:
		return (
			((w0*w3) % REDUCE),
			((w2*w1) % REDUCE),
			((w1*w3) % REDUCE),
			((w0*w2) % REDUCE) 
			)
	
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
# 2) MAP both halves so u get 2 points represented with X,Y,Z,T coords (Extended edwards coords)
# 3) perform addition of 2 edwards point, note that we need to add 2 edwards points so fe25519 arithmetics wont fit there
# we need to use function that adds 2 edwards points
def hash_to_group(input):
	# divide input into 2 halves, note that "input" paarameter needs to be in hexa-string
	# format, e.g. "5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6"
	# for more examples please check "INPUT_VECTORS_HASH_TO_GROUP_STATIC_HEXSTRING" in vectors.py
	
	t1,t2 = hash_to_num(input) # hash_to_num defined in convertLib.py -> returns 2 inteegers (first and second half)
	t1,t2 = t1% (2**255), t2%(2**255) # masking LSB: according to draft (4.3.4. Element derivation) we need to mask LSB which is equivalent to modulo 2**255
	
	a = MAP(t1) # map(ristretto_elligator) first half
	b = MAP(t2) # map(ristretto_elligator) second half 

	# a + b, this is addition of 2 points on edwards curve 
	# inspired by C cryptographic library "libsodium"
	# note that output "r" has coords X,Y,Z,T
	r = ge25519_p3_add(a,b) # add from libsodium -> libsodium_add.py

	R = ristretto255_encode(*r) # encode "r"" so you get bytes from point
	
	# convert and print result, note that "R" is inteeger in our internal representation
	# so conversion is needed to get bytes [32x8] or [16x16]
	# for more details about conversion please check convertionLib.py 
	r_bytes = numToHex(R,NUMBER_INTERPRETATION_CHOICES["32x8"])
	#numToHex(R,NUMBER_INTERPRETATION_CHOICES["16x16"],True)

	return r,r_bytes



def swap25519(a,b,c):
	a,b = (b,a) if c else (a,b)
	return a,b

# NOTE: that p,q are ristretto255 points (X,Y,Z,T)
def cswaps(p,q,b):
	X1,Y1,Z1,T1 = p
	X2,Y2,Z2,T2 = q
	X1,X2 = swap25519(X1,X2,b)
	Y1,Y2 = swap25519(Y1,Y2,b)
	Z1,Z2 = swap25519(Z1,Z2,b)
	T1,T2 = swap25519(T1,T2,b)
	return (X1,Y1,Z1,T1),(X2,Y2,Z2,T2)

# NOTE: that s is int, and q is ristretto point (X,Y,Z,T)
def ristretto255_scalarmult(q,s):
	s = numToHex(s,NUMBER_INTERPRETATION_CHOICES["32x8"],False)
	p = 0,1,1,0
	b=[]
	for i in range(255-1,-1,-1):
		b = (s[i//8]>>(i&7))&1
		p,q = cswaps(p,q,b)
		q=ge25519_p3_add(q,p)
		p=ge25519_p3_add(p,p)
		p,q = cswaps(p,q,b)
	return p

L =  [0xeb, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10]

def decorator(f):
	print("-------------")
	print(f)
	print("-------------")

def scalar_bit(in_a, i):
	in_a_str = bin(in_a)[2::][::-1]
	if i<0:
		return 0
	return int(in_a_str[i])

def Xoruj(a,b):
	# a xor b
	if len(a)!=len(b):
		raise ValueError(f"a:{len(a)} in not eq to b{len(b)}")
	if not isinstance(a[0],int):
		a = [int(i,16) for i in a]
	if not isinstance(b[0],int):
		b = [int(i,16) for i in b]


	skuska = [a[i]^b[i] for i in range(len(a))]
	
	decorator(skuska)
	return skuska

Lnum = hexToNum(L,NUMBER_INTERPRETATION_CHOICES["32x8"])

def mul_l(a,b):
		return (a*b) % (Lnum + 2)

def inverse_mod_l(a_in):
	m_inv = 1

	for i in range(252,-1,-1):
		m_inv = mul_l(m_inv,m_inv)

		if scalar_bit(Lnum,i):
			m_inv = mul_l(m_inv,a_in)

	return m_inv

"""
IN_hex = [255 for i in range(32)]
k = [0x57 for i in range(32)]
k_num = hexToNum(k,NUMBER_INTERPRETATION_CHOICES["32x8"],False)
IN = hexToNum(IN_hex,NUMBER_INTERPRETATION_CHOICES["32x8"],False)

print("TEST")
for i in range(1000):
	r = inverse_mod_l(IN)

	r_32 = numToHex(r,NUMBER_INTERPRETATION_CHOICES["32x8"],False)
	#print(r_32)
	
	IN = hexToNum(Xoruj(r_32,k),NUMBER_INTERPRETATION_CHOICES["32x8"],False)

numToHex(IN,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
print("TEST OUTPUT")


#a = mul_l(IN,IN)
#r = inverse_mod_l(IN)
#numToHex(r,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
"""

## TESTS

# TESTING A.1 -> small multiples of generator
TEST_RESULT = 1
GEN = [0xe2, 0xf2, 0xae, 0xa, 0x6a, 0xbc, 0x4e, 0x71, 0xa8, 0x84, 0xa9, 0x61,0xc5, 0x0, 0x51, 0x5f,0x58, 0xe3, 0xb, 0x6a, 0xa5, 0x82, 0xdd, 0x8d, 0xb6, 0xa6, 0x59, 0x45, 0xe0, 0x8d, 0x2d, 0x76]
GEN_int = hexToNum(GEN,NUMBER_INTERPRETATION_CHOICES["32x8"])
for i in range(16):
	GEN_ristretto255_point = ristretto255_decode(GEN_int)
	q = ristretto255_scalarmult(GEN_ristretto255_point,i)
	q_encoded = ristretto255_encode(*q)
	TEST_RESULT  &= numToHex(q_encoded,NUMBER_INTERPRETATION_CHOICES["32x8"]) == vectors.SMALL_MULTIPLES_OF_GENERATOR_VECTORS[i]
	print(TEST_RESULT)
print("vysledok po small mult",TEST_RESULT)


# testing hash_to_group
a = 0
HASH_VECTOR = vectors.INPUT_VECTORS_HASH_TO_GROUP_STATIC_HEXSTRING # vector from vectors.py	
for i,vec in enumerate(HASH_VECTOR):
	#print([hex(j) for j in bytes.fromhex(i)],",")
	_, r_bytes=hash_to_group(vec)
	TEST_RESULT  &= r_bytes == vectors.OUTPUT_VECTORS_HASH_TO_GROUP_BYTES[i]
	print(TEST_RESULT)
print("vysledok po hash to group",TEST_RESULT)


# testing if vector is negative 
for i in range(8):
	vec = hexToNum(vectors.test_negative_vectors[i])
	TEST_RESULT &= is_neg(vec)
	print(TEST_RESULT)
print("vysledok po testing if vector is negativ",TEST_RESULT)

#testing if vector is P-negative is positive 
for i in range(8):
	vec = hexToNum(vectors.test_negative_vectors_compl[i])
	TEST_RESULT &= 1-is_neg(vec)
	print(TEST_RESULT)
print("vysledok po testing if vector is P-negative is positive",TEST_RESULT)

#testing non-canonical vectors -> should result in error during decoding
for i in range(4):
	vec = hexToNum(vectors.non_canonical_vectors[i])
	try:
		ristretto255_decode(vec)
		TEST_RESULT &= 0
	except Exception as e:
		if e == ValueError:
			TEST_RESULT &= 1
	print(TEST_RESULT)
print("vysledok po testing non-canonical vectors",TEST_RESULT)
	
#testing Non-square x^2 vectors -> should result in error during decoding
for i in range(8):
	vec = hexToNum(vectors.non_square_x2[i])
	try:
		ristretto255_decode(vec)
		TEST_RESULT &= 0
	except Exception as e:
		if e == ValueError:
			TEST_RESULT &= 1
	print(TEST_RESULT)
print("vysledok po testing Non-square x^2 vectors",TEST_RESULT)

#testing Negative xy value vectors -> should result in error during decoding
for i in range(8):
	vec = hexToNum(vectors.negative_xy[i])
	try:
		ristretto255_decode(vec)
		TEST_RESULT &= 0
	except Exception as e:
		if e == ValueError:
			TEST_RESULT &= 1
	print(TEST_RESULT)
print("vysledok po testing Negative xy value vectors",TEST_RESULT)

MSG(f'Final Test result: {"SUCCESS" if TEST_RESULT else "FAIL"}')