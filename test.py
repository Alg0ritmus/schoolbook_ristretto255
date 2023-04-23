from convertLib import *
MINUS_ONE = 115565932813024562229384322928592814283244066726840484812818018414147674303743
P = (2**255)-19 
P_HEX = [
	0xec, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0x7f,
]

# Number test vectors

# Hex test vectors

HEX_1 = [
0x01, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00,
]

HEX_2 = [
0x01, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0xff, 
0xff, 0xff, 0xff, 0x7f,
]

HEX_3 = [
0xed, 0x57, 0xff, 0xd8, 
0xc9, 0x14, 0xfb, 0x20, 
0x14, 0x71, 0xd1, 0xc3, 
0xd2, 0x45, 0xce, 0x3c, 
0x74, 0x6f, 0xcb, 0xe6, 
0x3a, 0x36, 0x79, 0xd5, 
0x1b, 0x6a, 0x51, 0x6e, 
0xbe, 0xbe, 0x0e, 0x20, 
]

HEX_4 = [
0xc3, 0x4c, 0x4e, 0x18, 
0x26, 0xe5, 0xd4, 0x03, 
0xb7, 0x8e, 0x24, 0x6e, 
0x88, 0xaa, 0x05, 0x1c, 
0x36, 0xcc, 0xf0, 0xaa, 
0xfe, 0xbf, 0xfe, 0x13, 
0x7d, 0x14, 0x8a, 0x2b, 
0xf9, 0x10, 0x45, 0x62, 
]

HEX_5 = [
0xc9, 0x40, 0xe5, 0xa4, 
0x40, 0x41, 0x57, 0xcf, 
0xb1, 0x62, 0x8b, 0x10, 
0x8d, 0xb0, 0x51, 0xa8, 
0xd4, 0x39, 0xe1, 0xa4, 
0x21, 0x39, 0x4e, 0xc4, 
0xeb, 0xcc, 0xb9, 0xec, 
0x92, 0xa8, 0xac, 0x78, 
]

HEX_6 = [
0x47, 0xcf, 0xc5, 0x49, 
0x7c, 0x53, 0xdc, 0x8e, 
0x61, 0xc9, 0x1d, 0x17, 
0xfd, 0x62, 0x6f, 0xfb, 
0x1c, 0x49, 0xe2, 0xbc, 
0xa9, 0x4e, 0xed, 0x05, 
0x22, 0x81, 0xb5, 0x10, 
0xb1, 0x11, 0x7a, 0x24, 
]

HEX_7 = [
0xf1, 0xc6, 0x16, 0x5d, 
0x33, 0x36, 0x73, 0x51, 
0xb0, 0xda, 0x8f, 0x6e, 
0x45, 0x11, 0x01, 0x0c, 
0x68, 0x17, 0x4a, 0x03, 
0xb6, 0x58, 0x12, 0x12, 
0xc7, 0x1c, 0x0e, 0x1d, 
0x02, 0x6c, 0x3c, 0x72, 
]

HEX_8 = [
0x87, 0x26, 0x0f, 0x7a, 
0x2f, 0x12, 0x49, 0x51, 
0x18, 0x36, 0x0f, 0x02, 
0xc2, 0x6a, 0x47, 0x0f, 
0x45, 0x0d, 0xad, 0xf3, 
0x4a, 0x41, 0x3d, 0x21, 
0x04, 0x2b, 0x43, 0xb9, 
0xd9, 0x3e, 0x13, 0x09,
]


HEX_ALL = [
	HEX_1,HEX_2,HEX_3,HEX_4,
	HEX_5,HEX_6,HEX_7,HEX_8,
]

DEC_ALL = [
57896044618658097711785492504343953926634992332820282019728792003956564819948,
236,
43395981139273074876507941764684574595537234604803325968249716566753441523712,
13447355864474782990064147865452522366103789347284224292039261867954817577770,
3313441639055640068739630253660246146533292736357745798877292175242976476964,
41397104624399426060398833257223205409790519106274891578395515334690214654118,
6225623610693817694070834656612723216338662007248481655499821746921891510524,
53791225109085577843242291480306666131143723373880160896372464448377688480102,
]



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

	# r2 
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



"""
print("\nend:\n")
u = 1
v = DEC_ALL[7]
r2 = inv_sqrt(u,v)
numToHex(r2,NUMBER_INTERPRETATION_CHOICES["32x8"])
"""

def ristretto255_decode(s):
	ss = (s*s) % P
	u1 = (1 + ss) % P
	u2 = (1 - ss) % P

	uu1 = (u1*u1) % P
	uu2 = (u2*u2) % P

	duu1 = (EDWARDS_D*uu1) % P
	duu1_neg = fneg(duu1)
	v = (duu1_neg - uu2) % P # add minus here due to draft
	vuu2 = (v*uu2) % P

	was_square,I = inv_sqrt(1,vuu2)

	Dx = (I*u2) % P
	Dxv = (Dx*v) % P
	Dy = (I*Dxv) % P

	
	sDx = (s*Dx) % P
	x = (sDx+sDx) % P

	"""
	print("\nx:\n")
	numToHex(x,NUMBER_INTERPRETATION_CHOICES["32x8"])
	print("\nx-------:\n")
	"""

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



x,y,z,t = ristretto255_decode(DEC_ALL[7])

print("\nX:\n")
numToHex(x,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
print("\nY:\n")
numToHex(y,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
print("\nZ:\n")
numToHex(z,NUMBER_INTERPRETATION_CHOICES["32x8"],True)
print("\nT:\n")
numToHex(t,NUMBER_INTERPRETATION_CHOICES["32x8"],True)


vys = ristretto255_encode(x,y,z,t)
print("\n RESULT after encode:\n")
numToHex(vys,NUMBER_INTERPRETATION_CHOICES["32x8"],True)

