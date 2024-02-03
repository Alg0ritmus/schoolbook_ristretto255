from constants import *
from convertLib import *

#ristretto255_from_hash - libsodium:
# https://github.com/jedisct1/libsodium/blob/b7aebe5a1ef46bbb1345e8570fd2e8cea64e587f/src/libsodium/crypto_core/ed25519/ref10/ed25519_ref10.c#L2965
# X,Y,Z,T
# 0,1,2,3
# edward's ARITHMETIC
def ge25519_p3_to_cached(p):
	YplusX = (p[1] + p[0]) % P
	YminusX = (p[1] - p[0]) % P
	z = p[2] 
	T2d = (p[3]* D2) % P
	return (YplusX,YminusX,z,T2d)

def ge25519_add_cached(p,q):
	r = [0,0,0,0]
		#X,Y,Z,T
	r[0] = (p[1] + p[0]) % P
	r[1] = (p[1] - p[0]) % P
	r[2] = (r[0] * q[0]) % P
	r[1] = (r[1] * q[1]) % P
	r[3] = (q[3] * p[3]) % P
	r[0] = (p[2] * q[2]) % P
	t0   = (r[0] + r[0]) % P
	r[0] = (r[2] - r[1]) % P
	r[1] = (r[2] + r[1]) % P
	r[2] = (t0   + r[3]) % P
	r[3] = (t0   - r[3]) % P
	
	return tuple(r)

def ge25519_p1p1_to_p3(p):
	r = [0,0,0,0]
	r[0] = (p[0] * p[3]) % REDUCE
	r[1] = (p[1] * p[2]) % REDUCE
	r[2] = (p[2] * p[3]) % REDUCE
	r[3] = (p[0] * p[1]) % REDUCE
	return tuple(r)

# input body (x,y,z,t)
def ge25519_p3_add(p,q):
	q_cached = ge25519_p3_to_cached(q)
	p1p1 = ge25519_add_cached(p,q_cached)
	r = ge25519_p1p1_to_p3(p1p1)
	
	#print("----TESTq_cached---\n")
	#numToHex(r[0],NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	#numToHex(r[1],NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	#numToHex(r[2],NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	#numToHex(r[3],NUMBER_INTERPRETATION_CHOICES["32x8"],True)
	return r #-> (x,y,z,t)


