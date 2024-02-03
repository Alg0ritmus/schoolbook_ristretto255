# Here are some functions used for conversion. We wanted to be able to convert between
# python's inteeger and field_elem / bytes. Note that field_elem is type from C impl. of Bernstein's TweetNacl.

NUMBER_INTERPRETATION_CHOICES = {
    "32x8": [32, 8],  # set this option to perform int/bytes conversion
    "16x16": [16, 16]  # set this option to perform int/field_elem conversion
}

def hexToNum(IN, interpret = NUMBER_INTERPRETATION_CHOICES["32x8"]):
    sum = int.from_bytes(b''.join(number.to_bytes(1, 'little') for number in IN),'little')

    return sum

def numToHex(BIG_NUMBER, interpret=NUMBER_INTERPRETATION_CHOICES["32x8"], b=False):
    return [i for i in BIG_NUMBER.to_bytes(32,'little')]

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


def eq(a):
    return hexToNum([int(f'0x{i}',16) for i in a.split()])