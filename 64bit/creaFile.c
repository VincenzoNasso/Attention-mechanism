#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX ds; 	// dataset
	MATRIX wq; 	// pesi WQ
	MATRIX wk; 	// pesi WK
	MATRIX wv; 	// pesi WV
	MATRIX out;	// matrice contenente risultato (N x nn)
	VECTOR bq; 	// pesi bq
	VECTOR bk; 	// pesi bk
	VECTOR bv; 	// pesi bv
	int N;		// numero di righe del dataset
	int s; 		// prima dimensione del tensore S
	int n; 		// seconda dimensione del tensore S
	int d; 		// terza dimensione del tensore S
	int ns; 	// numero di tensori nel dataset
	int nn;		// numero di neuroni
	int display;
	int silent;
} params;


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
	int N = 2048, d = 48, fattore = 2000;
    MATRIX ds = alloc_matrix(N,d);
	ds = load_data("test_2048_48_64.ds",&N,&d);
	MATRIX tmp = alloc_matrix(N*fattore,d);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < N*d; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_2048_48_64.ds2", tmp, N*fattore, d);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

void creaDSOutput(){
	int N = 2048, d = 48, fattore = 2100;
    MATRIX ds = alloc_matrix(N,d);
	ds = load_data("test_2048_48_64.os",&N,&d);
	MATRIX tmp = alloc_matrix(N*fattore,d);
	int pos = 0;
	for(int i = 0; i < fattore; i++){
		for(int j = 0; j < N*d; j++)
			tmp[pos++] = ds[j];
	}
	save_data("test_2048_48_64.os2", tmp, N*fattore, d);
	dealloc_matrix(ds);
	dealloc_matrix(tmp);
}

int main(){
	
	creaDSInput();
	//creaDSOutput();
    return 0;
}