from convertLib import *
from constants import *

from bernstein import Fadd
from libsodium_add import ge25519_p3_add
from edwardsc_paper_add import edwards_addition



# curve25519 == fe25519 ARITHMETIC, ristretto255 uses this internally
def fneg(a):
	return (P-a)

def is_neg(a):
	arr_a = numToHex(a,NUMBER_INTERPRETATION_CHOICES["32x8"],False)
	return arr_a[0] & 1

def fmul(a,b):
	out = (a * b) % P
	return out

def pow_xtimes(a,n):
	out = a
	for i in range(n):
		out=fmul(out,out)
	return out

def fabs(x):
	return fneg(x) if is_neg(x) else x

# curve25519 == fe25519 SPECIAL FUNCTIONS

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



# Ristretto255 functions

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


# MAP == MAP2 len inak napisana
def MAP(t): # also known as ristretto255_elligator
	
	tt = (t*t) % P
	r =  (SQRT_M1 * tt) % P
	r1 =  (r +1) % P
	u = (r1 * ONE_MINUS_D_SQ) % P
	rD = (r*EDWARDS_D) % P
	r_plus_D = (r+EDWARDS_D) % P
	_1 = fneg(1) 
	_1_rD = (_1 - rD) % P
	v = (_1_rD * r_plus_D) % P

	was_square, s = inv_sqrt(u,v)

	st = (s*t) % P
	

	if is_neg(st):
		s_prime = fneg(st)
	else:
		s_prime = st

	s_prime = fneg(s_prime)
	

	if was_square:
		s = s
		c = _1
	else:
		s = s_prime
		c = r

	
	r_1 = (r - 1) % P
	cr_1 = (c*r_1) % P
	temp = (cr_1 * D_MINUS_ONE_SQ) % P
	N = (temp-v) % P


	ss = (s*s) % P
	s_ = (s+s) % P

	w0 =  (s_ * v) % P
	w1 =  (N * SQRT_AD_MINUS_ONE ) % P
	w2 = (1 - ss ) % P
	w3 = (1 + ss) % P

	
	return (
		((w0*w3) % P),
		((w2*w1) % P),
		((w1*w3) % P),
		((w0*w2) % P) 
		)



def MAP2(t): # also known as ristretto255_elligator 
	
	r = (SQRT_M1 * t * t) %P 
	u = ((r + 1) * ONE_MINUS_D_SQ) % P
	v = ((-1 - r*EDWARDS_D) * (r + EDWARDS_D)) % P

	was_square, s = inv_sqrt(u,v)
	st = (s*t) % P
	s_prime = fneg(fabs(st))
	

	if  was_square:
		s = s
		c = -1
	else:
		s = s_prime
		c = r 

	N = (c* (r -1 )) % P
	N = (N * D_MINUS_ONE_SQ) % P
	N = ( N - v) % P
	
	w0 = (2 * s * v) % P
	w1 = (N * SQRT_AD_MINUS_ONE) % P
	w2 = (1 - s*s) % P
	w3 = (1 + s*s) % P


	return (
		((w0*w3) % P),
		((w2*w1) % P),
		((w1*w3) % P),
		((w0*w2) % P) 
		)



# source: unknown
def special_addition(W,V):

	A = ((W[1] - W[0]) * (V[1] - V[0])) % P
	B = ((W[1] + W[0]) * (V[1] + V[0])) % P
	C = (2 * EDWARDS_D * W[2] * V[2]) % P
	D = (2 * W[3] * V[3]) % P
	E = (B - A) % P
	F = (D - C) % P

	# Compute the output point
	x3 = ((E*3 - F*F) * (W[2] * V[2])) % P
	y3 = (E * (F * (W[0] * V[1] + W[1] * V[0]) + A * B)) % P
	z3 = (F * E * W[2] * V[2]) % P
	t3 = (E * (F * (W[1] * V[1] + W[0] * V[0]) - A * B)) % P

	R = (x3, y3, t3, z3)
	return R



def hash_to_group(input):
	# rozdel hash_input na 2 polky a kazdu polku premen na cislo (python podporuje aj aritmetuku s velkymi cislami)
	t1,t2 = hash_to_num(input)
	t1,t2 = t1%P, t2%P
	
	a = MAP2(t1) # MAP == elligator
	b = MAP2(t2)

	# scitanie
	r = Fadd(a,b) # ADD z TweetNaCl -> bernstein.py
	r2 = edwards_addition(a,b) # scitanie z paper -> edwardsc_paper_add.py
	r3 = ge25519_p3_add(a,b) # add z libsodia -> libsodium_add.py

	R = ristretto255_encode(*r) 
	
	# iba vypis 
	numToHex(R,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	numToHex(R,NUMBER_INTERPRETATION_CHOICES["16x16"],True)
	
HASH_VECTOR = "5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6"
hash_to_group(HASH_VECTOR)




