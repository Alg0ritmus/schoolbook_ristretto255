#include "ristretto255.h"

int main(){

	const u8 in[32] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,255,255};
	u8 in_inv[32];
	field_elem out;
	field_elem out_inv;
	field_elem result;
	field_elem inv_sq_field_elem;
	unpack25519(out,in);
	// out inv
	finverse(out_inv,out);
	// out * out_inv = 1?
	fmul(result,out,out_inv);

	printf("\nPrint out:");
	print(out);

	printf("\nPrint out_inv:");
	print(out_inv);

	printf("\nPrint result:");
	print(result);

	printf("\nPrint IN:");
	print_32(in);


	pack25519(in_inv,result);
	printf("\nPrint IN_inv:");
	print_32(in_inv);

	// inv_sqrt( out, u, const, v)
	printf("inv_sqrt( out, u, const, v)\n");
	inv_sqrt(inv_sq_field_elem,out_inv,out);

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