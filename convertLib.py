NUMBER_INTERPRETATION_CHOICES={
	"32x8" : [32,8],
	"16x16" : [16,16]
}

P = (2**255)-19 

	

def hexToNum(input,interpret):
	sum = 0

	for i in range(interpret[0]):
		m = 2**(interpret[1]*i)

		sum += input[i]*m


	print(sum,"\ninverse number:", P-sum)


def numToHex(BIG_NUMBER,interpret):
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





