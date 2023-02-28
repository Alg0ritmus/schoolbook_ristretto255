#include "ristretto255.h"

static void unpack25519(field_elem out, const u8 *in)
{
int i;
for (i = 0; i < 16; ++i) out[i] = in[2*i] + ((i64) in[2*i + 1] << 8);
out[15] &= 0x7fff;
 }

 static void carry25519(field_elem elem)
 {
 int i;
 i64 carry;
 for (i = 0; i < 16; ++i) {
 carry = elem[i] >> 16;
 elem[i] -= carry << 16;
 if (i < 15) elem[i + 1] += carry; else elem[0] += 38 * carry;
 }
 }

 static void fadd(field_elem out, const field_elem a, const field_elem b) /* out = a + b */
 {
 int i;
 for (i = 0; i < 16; ++i) out[i] = a[i] + b[i];
 }

 static void fsub(field_elem out, const field_elem a, const field_elem b) /* out = a - b */
 {
 int i;
 for (i = 0; i < 16; ++i) out[i] = a[i] - b[i];
 }

 static void fmul(field_elem out, const field_elem a, const field_elem b) /* out = a * b */
 {
 i64 i, j, product[31];
 for (i = 0; i < 31; ++i) product[i] = 0;
 for (i = 0; i < 16; ++i) {
 for (j = 0; j < 16; ++j) product[i+j] += a[i] * b[j];
 }
 for (i = 0; i < 15; ++i) product[i] += 38 * product[i+16];
 for (i = 0; i < 16; ++i) out[i] = product[i];
 carry25519(out);
 carry25519(out);
 }

static void finverse(field_elem out, const field_elem in)
{
field_elem c;
int i;
for (i = 0; i < 16; ++i) c[i] = in[i];
 for (i = 253; i >= 0; i--) {
 fmul(c,c,c);
 if (i != 2 && i != 4) fmul(c,c,in);
 }
 for (i = 0; i < 16; ++i) out[i] = c[i];
 }

 static void swap25519(field_elem p, field_elem q, int bit)
 {
 i64 t, i, c = ~(bit - 1);
 for (i = 0; i < 16; ++i) {
 t = c & (p[i] ^ q[i]);
 p[i] ^= t;
 q[i] ^= t;
 }
 }

 static void pack25519(u8 *out, const field_elem in)
 {
 int i, j, carry;
 field_elem m, t;
 for (i = 0; i < 16; ++i) t[i] = in[i];
 carry25519(t); carry25519(t); carry25519(t);
 for (j = 0; j < 2; ++j) {
 m[0] = t[0] - 0xffed;
 for(i = 1; i < 15; i++) {
 m[i] = t[i] - 0xffff - ((m[i-1] >> 16) & 1);
 m[i-1] &= 0xffff;
 }
 m[15] = t[15] - 0x7fff - ((m[14] >> 16) & 1);
 carry = (m[15] >> 16) & 1;
 m[14] &= 0xffff;
 swap25519(t, m, 1 - carry);
 }
 for (i = 0; i < 16; ++i) {
 out[2*i] = t[i] & 0xff;
 out[2*i + 1] = t[i] >> 8;
 }
 }


void scalarmult(u8 *out, const u8 *scalar, const u8 *point)
{
u8 clamped[32];
i64 bit, i;
field_elem a, b, c, d, e, f, x;
 for (i = 0; i < 31; ++i) clamped[i] = scalar[i];
 clamped[0] &= 0xf8;
 clamped[31] = (clamped[31] & 0x7f) | 0x40;
 unpack25519(x, point);
 for (i = 0; i < 16; ++i) {
 b[i] = x[i];
 d[i] = a[i] = c[i] = 0;
 }
 a[0] = d[0] = 1;
 for (i = 254; i >= 0; --i) {
 bit = (clamped[i >> 3] >> (i & 7)) & 1;
 swap25519(a, b, bit);
 swap25519(c, d, bit);
 fadd(e, a, c);
 fsub(a, a, c);
 fadd(c, b, d);
 fsub(b, b, d);
 fmul(d, e, e);
 fmul(f, a, a);
 fmul(a, c, a);
 fmul(c, b, e);
 fadd(e, a, c);
 fsub(a, a, c);
 fmul(b, a, a);
 fsub(c, d, f);
 fmul(a, c, _121665);
 fadd(a, a, d);
 fmul(c, c, a);
 fmul(a, d, f);
 fmul(d, b, x);
 fmul(b, e, e);
 swap25519(a, b, bit);
 swap25519(c, d, bit);
 }
 finverse(c, c);
 fmul(a, a, c);
 pack25519(out, a);
 }
