SQRT_M1 = 57896044618658097711785492504343953926634992332820282019728792003956564819949


# choose between 16 or 32 element representation (just uncomment preferable option) 
interpret= [32,8]
#interpret= [16,16]
sum = 0

result=[]

for i in range(interpret[0]-1,-1,-1):
  m = 2**(interpret[1]*i)

  a = SQRT_M1 // m
  SQRT_M1 = SQRT_M1 % m

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