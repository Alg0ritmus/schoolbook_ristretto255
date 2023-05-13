from constants import *
# check 3.1 section
# https://eprint.iacr.org/2008/522.pdf
def edwards_addition(W,V):
	X1,Y1,Z1,T1, =  W
	X2,Y2,Z2,T2, =  V
	a = 1

	X3 = ((X1*Y2 + Y1*X2)*(Z1*Z2 - EDWARDS_D*T1*T2)) % P
	Y3 = ((Y1*Y2 - a*X1*X2)*(Z1*Z2 + EDWARDS_D*T1*T2)) % P
	T3 = ((Y1*Y2 - a*X1*X2)*(X1*Y2 + Y1*X2)) % P
	Z3 = ((Z1*Z2 - EDWARDS_D*T1*T2)*(Z1*Z2 + EDWARDS_D*T1*T2)) % P

	return (X3,Y3,Z3,T3)
