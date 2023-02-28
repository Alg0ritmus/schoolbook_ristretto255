typedef unsigned char u8;
typedef long long i64;
typedef i64 field_elem[16];
static const field_elem _121665 = {0xDB41, 1},

static void unpack25519(field_elem out, const u8 *in);
static void carry25519(field_elem elem);
static void fadd(field_elem out, const field_elem a, const field_elem b);
static void fsub(field_elem out, const field_elem a, const field_elem b);
static void fmul(field_elem out, const field_elem a, const field_elem b);/* out = a * b */
static void finverse(field_elem out, const field_elem in);
static void swap25519(field_elem p, field_elem q, int bit);
static void pack25519(u8 *out, const field_elem in);
void scalarmult(u8 *out, const u8 *scalar, const u8 *point);