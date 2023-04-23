NUMBER_INTERPRETATION_CHOICES={
	"32x8" : [32,8],
	"16x16" : [16,16]
}

P = (2**255)-19 

SQRT_M1 = 19681161376707505956807079304988542015446066515923890162744021073123829784752	
EDWARDS_D = 37095705934669439343138083508754565189542113879843219016388785533085940283555
INVSQRT_A_MINUS_D = 54469307008909316920995813868745141605393597292927456921205312896311721017578
def hexToNum(input,interpret,b):
	sum = 0

	for i in range(interpret[0]):
		m = 2**(interpret[1]*i)

		sum += input[i]*m


	#print(sum,"\ninverse number:", P-sum)
	return sum


def numToHex(BIG_NUMBER,interpret,b):
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
	"""
	for i in range(len(result)):
	  if i%4==0:
	    print()
	  print(result[i], end=", ")

	print("\n\n\n")

	"""
	if b == True:
		for i in range(len(result)):
		  if i%4==0:
		    print()
		  print(hex(result[i]), end=", ")
		print("\n")

	return result






def MSG(message):
	print("\n\n")
	print("*"*(len(message)+3))
	print("|  "+message+"|")
	print("*"*(len(message)+3))
	print("\n\n")
