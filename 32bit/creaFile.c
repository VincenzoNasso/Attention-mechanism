#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
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

void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

void creaDSInput(){
	int r, c, fattore = 4000;
    MATRIX ds = alloc_matrix(r,c);
	ds = load_data("test_2048_48_32.ds",&r,&c);
	MATRIX tmp = alloc_matrix(r*fattore,c);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_2048_48_32.ds2", tmp, r*fattore, c);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaDSOutput(){
	int r, c, fattore = 4000;
    MATRIX ds = alloc_matrix(r,c);
	ds = load_data("test_2048_48_32.os",&r,&c);
	MATRIX tmp = alloc_matrix(r*fattore,c);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_2048_48_32.os2", tmp, r*fattore, c);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaPesiQ(){
	int r = 48, c = 32, fattore = 4;
    MATRIX ds = alloc_matrix(r,c);
	ds = load_data("test_48_32_32.wq",&r,&c);
	MATRIX tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_48_32_32.wq2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaPesiK(){
	int r = 48, c = 32, fattore = 4;
    MATRIX ds = alloc_matrix(r,c);
	ds = load_data("test_48_32_32.wk",&r,&c);
	MATRIX tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_48_32_32.wk2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaPesiV(){
	int r = 48, c = 32, fattore = 4;
    MATRIX ds = alloc_matrix(r,c);
	ds = load_data("test_48_32_32.wv",&r,&c);
	MATRIX tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_48_32_32.wv2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaBiasQ(){
	int r = 1, c = 32, fattore = 4;
    VECTOR ds = alloc_matrix(r,c);
	ds = load_data("test_32_32.bq",&r,&c);
	VECTOR tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_32_32.bq2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaBiasK(){
	int r = 1, c = 32, fattore = 4;
    VECTOR ds = alloc_matrix(r,c);
	ds = load_data("test_32_32.bk",&r,&c);
	VECTOR tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_32_32.bk2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaBiasV(){
	int r = 1, c = 32, fattore = 4;
    VECTOR ds = alloc_matrix(r,c);
	ds = load_data("test_32_32.bv",&r,&c);
	VECTOR tmp = alloc_matrix(r,c*fattore);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < r*c; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_32_32.bv2", tmp, r, c*fattore);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}


int main(){
	
	creaDSInput();
	/*
	creaPesiQ();
	creaPesiK();
	creaPesiV();
	creaBiasQ();
	creaBiasK();
	creaBiasV();
	creaDSOutput();
	*/
    return 0;
}