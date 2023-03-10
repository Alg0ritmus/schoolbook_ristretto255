# ristretto pre-computed constants: 
# https://datatracker.ietf.org/doc/id/draft-irtf-cfrg-ristretto255-00.html#name-internal-utility-functions-
SQRT_M1 = 19681161376707505956807079304988542015446066515923890162744021073123829784752
D = 37095705934669439343138083508754565189542113879843219016388785533085940283555
SQRT_AD_MINUS_ONE = 25063068953384623474111414158702152701244531502492656460079210482610430750235
INVSQRT_A_MINUS_D = 54469307008909316920995813868745141605393597292927456921205312896311721017578
ONE_MINUS_D_SQ = 1159843021668779879193775521855586647937357759715417654439879720876111806838
D_MINUS_ONE_SQ = 40440834346308536858101042469323190826248399146238708352240133220865137265952

P = (2**255)-19 

# select number you want to convert to hex repr.
BIG_NUMBER = INVSQRT_A_MINUS_D

# choose between 16 or 32 element representation (just uncomment preferable option) 
#interpret= [32,8]
interpret= [16,16]
sum = 0

result=[]

for i in range(interpret[0]-1,-1,-1):
  m = 2**(interpret[1]*i)

  a = BIG_NUMBER // m
  BIG_NUMBER = BIG_NUMBER % m

  result.append(a)

result = result[::-1]
#########
# PRINT #
#########
for i in range(len(result)):
  if i%4==0:
    print()
  print(result[i], end=", ")

print("\n\n\n")

for i in range(len(result)):
  if i%4==0:
    print()
  print(hex(result[i]), end=", ")