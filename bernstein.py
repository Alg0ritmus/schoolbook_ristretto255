from constants import *

# ADD z tweetNaCl
# https://github.com/dominictarr/tweetnacl/blob/master/tweetnacl.c#L590
def Fadd(p,q):
 
  a = Z(p[1], p[0])
  t = Z(q[1], q[0])
  a = M(a, t)
  b = A(p[0], p[1])
  t = A(q[0], q[1])
  b = M(b, t)
  c = M(p[3], q[3])
  c = M(c, D2)
  d =M(p[2], q[2])
  d =A(d, d)
  e = Z(b, a)
  f = Z(d, c)
  g = A(d, c)
  h = A(b, a)


  p0 = M(e, f);
  p1 = M(h, g);
  p2 = M(g, f);
  p3 = M(e, h);
  return (p0,p1,p2,p3)


def A(a,b):
  return ((a+b) % P)


def Z(a,b):
  return ((a-b) % P)

def M(a,b):
  return ((a*b) % P)

