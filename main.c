#include "ristretto255.h"



int main(){

	const u8 in[32] = {2,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16};
	field_elem out;
	field_elem inv_sq_field_elem;

	unpack25519(out,in);
	
	printf("inv_sqrt( out, u, const, v)\n");
	inv_sqrt(inv_sq_field_elem,out,out);

	// TEST
	field_elem t2, t4, t8, tx2, tx8;

	//a^2*1  -> (a^2*1)^2*2 -> (a^2)^4 = a^8
	pow_xtimes(tx2,out,1);
	pow_xtimes(tx8,tx2,2);


	// ((a^2)^2)^2 = a^8
	fmul(t2,out,out);
	fmul(t4,t2,t2);
	fmul(t8,t4,t4);

	printf("POROVNAVAM!!!\n");
	print(tx8);
	print(t8);

	int skuska = feq(tx8,t8);
	printf("skuska:%d\n",skuska );

	// skuska modulus
	
	field_elem A;
	u8 A_pack[32];
	u8 A_pack_out[32];
	fcopy(A,_121665);
	pack25519(A_pack,A);
	print_32(A_pack);

	ristretto255_point *temp_rp;
	ristretto255_decode(&temp_rp, in);
	printf("point cords \n");
	print(&temp_rp->x);
	print(&temp_rp->y);
	print(&temp_rp->z);
	print(&temp_rp->t);


	printf("\nA_pack_out before encoding:\n");
	print_32(A_pack_out);
	ristretto255_encode(A_pack_out, temp_rp);
	printf("\nA_pack_out after encoding:\n");
	print_32(A_pack_out);

	

	return 0;
}


