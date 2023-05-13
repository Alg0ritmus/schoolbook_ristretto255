NUMBER_INTERPRETATION_CHOICES={
	"32x8" : [32,8],
	"16x16" : [16,16]
}


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



# helper function -> make 2 halves from hash_string and convert it into 2 arrays of bytes
def hash_to_bytes(str_hash):
	P1,P2 = str_hash[:64],str_hash[64:]
	P1_bytes, P2_bytes = [],[]
	for i in range(2,66,2):
		P1_bytes.append((int(P1[i-2:i],16)))
		P2_bytes.append((int(P2[i-2:i],16)))

	
	return (P1_bytes,P2_bytes)


# helper function -> make 2 halves from hash_string and convert it into 2 numbers
def hash_to_num(str_hash):
	x1,x2 = hash_to_bytes(str_hash)
	x1,x2 = hexToNum(x1,NUMBER_INTERPRETATION_CHOICES["32x8"],False),hexToNum(x2,NUMBER_INTERPRETATION_CHOICES["32x8"],False)

	return (x1,x2)



# poznamky pre mna:
# t1= [93, 27, 224, 158, 61, 12, 130, 252, 83, 129, 18, 73, 14, 53, 112, 25, 121, 217, 158, 6, 202, 62, 43, 91, 84, 191, 254, 139, 77, 199, 114, 193]



# TEST VECTOR hash to group:
#5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6

# IN = bytes.fromhex("5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6")
# import ge25519 as g
# g.ge25519_p3.from_hash_ristretto255(IN).hex()