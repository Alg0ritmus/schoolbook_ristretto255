# Here are some functions used for conversion. We wanted to be able to convert between
# python's inteeger and field_elem / bytes. Note that field_elem is type from C impl. of Bernstein's TweetNacl.

NUMBER_INTERPRETATION_CHOICES = {
    "32x8": [32, 8],  # set this option to perform int/bytes conversion
    "16x16": [16, 16]  # set this option to perform int/field_elem conversion
}


# USAGE
# for bytes[32] to int -> use NUMBER_INTERPRETATION_CHOICES["32x8"]
# hexToNum(input,NUMBER_INTERPRETATION_CHOICES["32x8"],False) -> set b=False if you dont want to print result
# hexToNum(input,NUMBER_INTERPRETATION_CHOICES["32x8"],True) -> set b=True if you want to print result

# for field_elem bytes[16] to int -> use NUMBER_INTERPRETATION_CHOICES["16x16"]
# hexToNum(input,NUMBER_INTERPRETATION_CHOICES["16x16"],False) -> set b=False if you dont want to print result
# hexToNum(input,NUMBER_INTERPRETATION_CHOICES["16x16"],True) -> set b=True if you want to print result
# if useed properly we can use it for unpack25519/pack25519 convertion compatible with TweetNaCl
def hexToNum(input, interpret = NUMBER_INTERPRETATION_CHOICES["32x8"]):
    sum = 0

    for i in range(interpret[0]):
        m = 2**(interpret[1]*i)

        sum += input[i]*m

    return sum

# USAGE
# for int to bytes[32] -> use NUMBER_INTERPRETATION_CHOICES["32x8"]
# numToHex(input,NUMBER_INTERPRETATION_CHOICES["32x8"],False) -> set b=False if you dont want to print result
# numToHex(input,NUMBER_INTERPRETATION_CHOICES["32x8"],True) -> set b=True if you want to print result

# for int to field_elem bytes[16] -> use NUMBER_INTERPRETATION_CHOICES["16x16"]
# numToHex(input,NUMBER_INTERPRETATION_CHOICES["16x16"],False) -> set b=False if you dont want to print result
# numToHex(input,NUMBER_INTERPRETATION_CHOICES["16x16"],True) -> set b=True if you want to print result


def numToHex(BIG_NUMBER, interpret, b):
    result = []
    for i in range(interpret[0]-1, -1, -1):
        m = 2**(interpret[1]*i)

        a = BIG_NUMBER // m
        BIG_NUMBER = BIG_NUMBER % m

        result.append(a)

    result = result[::-1]

    if b == True:
        for i in range(len(result)):
            if i % 100 == 0:
                print()
            print(hex(result[i]), end=", ")
        print("\n")

    return result


def numToBytes(BIG_NUMBER, interpret):
    result = []
    for i in range(interpret[0]-1, -1, -1):
        m = 2**(interpret[1]*i)

        a = BIG_NUMBER // m
        BIG_NUMBER = BIG_NUMBER % m

        result.append(a)

    result = result[::-1]

    return result


def MSG(message):
    print("\n\n")
    print("*"*(len(message)+3))
    print("|  "+message+"|")
    print("*"*(len(message)+3))
    print("\n\n")


# helper function -> make 2 halves from hash_string and convert it into 2 arrays of bytes
def hash_to_bytes(str_hash):
    P1, P2 = str_hash[:64], str_hash[64:]
    P1_bytes, P2_bytes = [], []
    for i in range(2, 66, 2):
        P1_bytes.append((int(P1[i-2:i], 16)))
        P2_bytes.append((int(P2[i-2:i], 16)))

    return (P1_bytes, P2_bytes)


# helper function -> make 2 halves from hash_string and convert it into 2 numbers
def hash_to_num(str_hash):
    a = str_hash[:64]
    b = str_hash[64:]

    a = bytes.fromhex(a)
    b = bytes.fromhex(b)

    a = [i for i in a]
    b = [i for i in b]

    a = hexToNum(a, NUMBER_INTERPRETATION_CHOICES["32x8"])
    b = hexToNum(b, NUMBER_INTERPRETATION_CHOICES["32x8"])

    return (a, b)
