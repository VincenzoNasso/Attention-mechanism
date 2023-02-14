#include <stdio.h>
#include <stdlib.h>
#include <mm_malloc.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	*n = rows;
	*k = cols;
	
	return data;
}

int main() {

   	int Np, nnp, Nn, nnn;
	MATRIX dsProf = load_data("test_2048_48_32.os",&Np,&nnp);
	MATRIX dsNostro = load_data("out32_2048_8_64_48.ds2",&Nn,&nnn);
	for(int i = 0; i < Np*nnp; i++) {
		printf("scarto[%d] = %f\n", i, dsProf[i]-dsNostro[i]);
	}

	dealloc_matrix(dsProf);
	dealloc_matrix(dsNostro);
}
