#include "ristretto255.h"

int main(){

	const u8 in[32] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,255,255};
	u8 in_inv[32];
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

	return 0;
}


/*

echo "# schoolbook_ristretto255" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/Alg0ritmus/schoolbook_ristretto255.git
git push -u origin main

*/