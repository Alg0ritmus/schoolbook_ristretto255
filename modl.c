#include "helpers.h" // included definition of types like u8, u32... + FOR macros, etc.
#include "modl.h" 

//////////////////
// MONOCYPHER - needed for efficient calculation of a mod L
/////////////////
/// Utilities ///
/////////////////
static void store32_le(u8 out[4], u32 in){
    out[0] =  in        & 0xff;
    out[1] = (in >>  8) & 0xff;
    out[2] = (in >> 16) & 0xff;
    out[3] = (in >> 24) & 0xff;
}

static void store32_le_buf(u8 *dst, const u32 *src, size_t size) {
    size_t i;
    FOR(i, 0, size) { store32_le(dst + i*4, src[i]); }
}

// get bit from scalar at position i
static int scalar_bit(const u8 s[BYTES_ELEM_SIZE], int i){
    if (i < 0) { return 0; } // handle -1 for sliding windows
    return (s[i>>3] >> (i&7)) & 1;
}


///////////////////////////
/// Arithmetic modulo L ///
///////////////////////////
static const u32 L[8] = {
    0x5cf5d3ed, 0x5812631a, 0xa2f79cd6, 0x14def9de,
    0x00000000, 0x00000000, 0x00000000, 0x10000000,
};

//  p = a*b + p
static void multiply(u32 p[16], const u32 a[8], const u32 b[8]){
   size_t i, j;
    FOR (i, 0, 8) {
        u64 carry = 0;
        FOR (j, 0, 8) {
            carry  += p[i+j] + (u64)a[i] * b[j];
            p[i+j]  = (u32)carry;
            carry >>= 32;
        }
        p[i+8] = (u32)carry;
    }
}

static int is_above_l(const u32 x[8]){
   size_t i;
    // We work with L directly, in a 2's complement encoding
    // (-L == ~L + 1)
    u64 carry = 1;
    FOR (i, 0, 8) {
        carry  += (u64)x[i] + (~L[i] & 0xffffffff);
        carry >>= 32;
    }
    return (int)carry; // carry is either 0 or 1
}

// Final reduction modulo L, by conditionally removing L.
// if x < l     , then r = x
// if l <= x 2*l, then r = x-l
// otherwise the result will be wrong
static void remove_l(u32 r[8], const u32 x[8]){
   size_t i;
    u64 carry = (u64)is_above_l(x);
    u32 mask  = ~(u32)carry + 1; // carry == 0 or 1
    FOR (i, 0, 8) {
        carry += (u64)x[i] + (~L[i] & mask);
        r[i]   = (u32)carry;
        carry >>= 32;
    }
}

// Full reduction modulo L (Barrett reduction)
void mod_l(u8 reduced[32], const u32 x[16]){
   size_t i, j;
    static const u32 r[9] = {
        0x0a2c131b,0xed9ce5a3,0x086329a7,0x2106215d,
        0xffffffeb,0xffffffff,0xffffffff,0xffffffff,0xf,
    };
    // xr = x * r
    u32 xr[25] = {0};
    FOR (i, 0, 9) {
        u64 carry = 0;
        FOR (j, 0, 16) {
            carry  += xr[i+j] + (u64)r[i] * x[j];
            xr[i+j] = (u32)carry;
            carry >>= 32;
        }
        xr[i+16] = (u32)carry;
    }
    // xr = floor(xr / 2^512) * L
    // Since the result is guaranteed to be below 2*L,
    // it is enough to only compute the first 256 bits.
    // The division is performed by saying xr[i+16]. (16 * 32 = 512)
    ZERO(i, xr, 8);
    FOR (i, 0, 8) {
        u64 carry = 0;
        FOR (j, 0, 8-i) {
            carry   += xr[i+j] + (u64)xr[i+16] * L[j];
            xr[i+j] = (u32)carry;
            carry >>= 32;
        }
    }
    // xr = x - xr
    u64 carry = 1;
    FOR (i, 0, 8) {
        carry  += (u64)x[i] + (~xr[i] & 0xffffffff);
        xr[i]   = (u32)carry;
        carry >>= 32;
    }
    // Final reduction modulo L (conditional subtraction)
    remove_l(xr, xr);
    store32_le_buf(reduced, xr, 8);

//  WIPE_BUFFER(xr);
}

/********************* MD *******************************************************************/
// !!! fnkcne len na CPU s Litle Endian architekturou!!! Pre Big Endian treba doplnit kopirovanie dat!!!
static void multiply_mod_l(u32 r[8], const u32 a[8], const u32 b[8]){
    u32 c[16] = {0};
        multiply(c, a, b);
        mod_l((u8*) r, c);
}


// kod s vyuzitim TMP priestoru na zasobniku 
void inverse_mod_l(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]){
    static const  u8 Lm2[BYTES_ELEM_SIZE] = { // L - 2
        0xeb, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10,
    };
    u32 m_inv [8] = {1};

    // Compute the inverse
    for (int i = 252; i >= 0; i--) {
        multiply_mod_l(m_inv, m_inv, m_inv);

        if (scalar_bit(Lm2, i)) {
            multiply_mod_l(m_inv, m_inv, (u32*) in);

        }
    }
    int i;
    COPY(i ,out, (u8*) m_inv, BYTES_ELEM_SIZE);
}