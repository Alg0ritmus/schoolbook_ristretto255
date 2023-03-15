# ristretto pre-computed constants: 
# https://datatracker.ietf.org/doc/id/draft-irtf-cfrg-ristretto255-00.html#name-internal-utility-functions-
SQRT_M1 = 57896044618658097711785492504343953926634992332820282019728792003956564819949
D = 37095705934669439343138083508754565189542113879843219016388785533085940283555
SQRT_AD_MINUS_ONE = 25063068953384623474111414158702152701244531502492656460079210482610430750235
INVSQRT_A_MINUS_D = 54469307008909316920995813868745141605393597292927456921205312896311721017578
ONE_MINUS_D_SQ = 1159843021668779879193775521855586647937357759715417654439879720876111806838
D_MINUS_ONE_SQ = 40440834346308536858101042469323190826248399146238708352240133220865137265952

MINUS_ONE_NUMBER = 115565932813024562229384322928592814283244066726840484812818018414147674303743

MINUS_ONE = [
  0xecff, 0xffff, 0xffff, 0xffff,
  0xffff, 0xffff, 0xffff, 0xffff,
  0xffff, 0xffff, 0xffff, 0xffff,
  0xffff, 0xffff, 0xffff, 0xff7f
  ]

One = 1

# variable "input" is hex number representation of 32xu8 or 16*u16
# if input is 32xu8 you need to select interpret. [32,8]
# if input is 16xu16 you need to select interpret. [16,16]
input = [
  0xec, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0x7f
]


# choose between 16 or 32 element representation (just uncomment preferable option) 
interpret= [32,8]
#interpret= [16,16]
sum = 0

for i in range(interpret[0]):
  m = 2**(interpret[1]*i)

  sum += input[i]*m


print(sum)


  
