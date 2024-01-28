#ifndef _MODL_H
#define _MODL_H

void mod_l(u8 reduced[BYTES_ELEM_SIZE], const u32 x[16]);
void inverse_mod_l(u8 out[BYTES_ELEM_SIZE], const u8 in[BYTES_ELEM_SIZE]);

#endif //_MODL_H