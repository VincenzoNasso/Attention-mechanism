#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <omp.h>

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

// PROCEDURE ASSEMBLY
void extern	primaFase_U1(MATRIX fetta, MATRIX pesi, VECTOR bias, int n,  int d, int nn, MATRIX ret); //no unrilling
void extern	primaFase_U2(MATRIX fetta, MATRIX pesi, VECTOR bias, int n,  int d, int nn, MATRIX ret); //2 unrolling
void extern	primaFase_U3(MATRIX fetta, MATRIX pesi, VECTOR bias, int n,  int d, int nn, MATRIX ret); //2 unrolling
void extern	primaFase_U4(MATRIX fetta, MATRIX pesi, VECTOR bias, int n,  int d, int nn, MATRIX ret); //4 unrolling

void extern secondaFase_U1(MATRIX matriceQ, MATRIX matriceK, MATRIX matriceSecondaTerzaFase, int n, int nn, type fattore);//no unrolling
void extern secondaFase_U2(MATRIX matriceQ, MATRIX matriceK, MATRIX matriceSecondaTerzaFase, int n, int nn, type fattore);//2 unrolling
void extern secondaFase_U3(MATRIX matriceQ, MATRIX matriceK, MATRIX matriceSecondaTerzaFase, int n, int nn, type fattore);//4 unrolling
void extern secondaFase_U4(MATRIX matriceQ, MATRIX matriceK, MATRIX matriceSecondaTerzaFase, int n, int nn, type fattore);//4 unrolling

void extern terzaFase(int n, MATRIX matriceSecondaTerzaFase);//2 unrolling
void extern terzaFase_U1(int n, MATRIX matriceSecondaTerzaFase);//no unrolling

void extern quartaFase_U1(MATRIX matriceSecondaTerzaFase, MATRIX matriceV, MATRIX out, int n, int nn, int inizio); //no unrolling
void extern quartaFase_U2(MATRIX matriceSecondaTerzaFase, MATRIX matriceV, MATRIX out, int n, int nn, int inizio);//2 unrolling
void extern quartaFase_U3(MATRIX matriceSecondaTerzaFase, MATRIX matriceV, MATRIX out, int n, int nn, int inizio);//2 unrolling
void extern quartaFase_U4(MATRIX matriceSecondaTerzaFase, MATRIX matriceV, MATRIX out, int n, int nn, int inizio);//4 unrolling

void estraiFetta(int indiceTensore, int j, params *input, MATRIX ret) {
	int c1 = 0;

	int inizio = (indiceTensore * input->s * input->n * input->d) + (j * input->n * input->d);
	int fine = inizio + input->n * input->d;

	for (int i = inizio; i < fine; i++){
		ret[c1++] = input->ds[i];
	}
}

void trasposta(MATRIX k, MATRIX ret, int n, int nn){
	for(int l = 0; l < nn; l++){
		for(int m = 0; m < n; m++){
			ret[l*n+m] = k[m*nn+l];	
		}
	}
}


void calcolaTensori(params *input, int indiceTensore) {
	
	MATRIX fetta = alloc_matrix(input->n, input->d);
	MATRIX matriceQ = alloc_matrix(input->n, input->nn);
	MATRIX matriceK = alloc_matrix(input->n, input->nn);
	MATRIX matriceV = alloc_matrix(input->n, input->nn);
	MATRIX matriceSecondaTerzaFase = alloc_matrix(input->n, input->n);
	MATRIX k_t = alloc_matrix(input->nn, input->n);
	type fattoreScalare = (1/sqrt(input->d));
	
	for(int j = 0; j < input->s; j++) {
		estraiFetta(indiceTensore, j, input, fetta);
		int inizioFettaFinale = 8*(indiceTensore * input->s * input->n * input->nn) + 8*(j * input->n * input->nn);

		if(input->n % 16 == 0 && input->d % 16 == 0 && input->nn % 16 == 0){
			
			primaFase_U4(fetta, input->wq, input->bq, input->n, input->d, input->nn, matriceQ); //4 unrolling
			primaFase_U4(fetta, input->wk, input->bk, input->n, input->d, input->nn, matriceK); //4 unrolling
			primaFase_U4(fetta, input->wv, input->bv, input->n, input->d, input->nn, matriceV); //4 unrolling

			trasposta(matriceK, k_t, input->n, input->nn);

			secondaFase_U4(matriceQ, k_t, matriceSecondaTerzaFase, input->n, input->nn, fattoreScalare);//4 unrolling
		
			terzaFase(input->n, matriceSecondaTerzaFase);

			quartaFase_U4(matriceSecondaTerzaFase,matriceV,input->out,input->n,input->nn, inizioFettaFinale);
		}
		else if(input->n % 12 == 0 && input->d % 12 == 0 && input->nn % 12 == 0){
			primaFase_U3(fetta, input->wq, input->bq, input->n, input->d, input->nn, matriceQ); //4 unrolling
			primaFase_U3(fetta, input->wk, input->bk, input->n, input->d, input->nn, matriceK); //4 unrolling
			primaFase_U3(fetta, input->wv, input->bv, input->n, input->d, input->nn, matriceV); //4 unrolling

			trasposta(matriceK, k_t, input->n, input->nn);

			secondaFase_U3(matriceQ, k_t, matriceSecondaTerzaFase, input->n, input->nn, fattoreScalare);//4 unrolling
		
			terzaFase_U1(input->n, matriceSecondaTerzaFase);

			quartaFase_U3(matriceSecondaTerzaFase,matriceV,input->out,input->n,input->nn, inizioFettaFinale);
		}
		else if(input->n % 8 == 0 && input->d % 8 == 0 && input->nn % 8 == 0){
			primaFase_U2(fetta, input->wq, input->bq, input->n, input->d, input->nn, matriceQ); //4 unrolling
			primaFase_U2(fetta, input->wk, input->bk, input->n, input->d, input->nn, matriceK); //4 unrolling
			primaFase_U2(fetta, input->wv, input->bv, input->n, input->d, input->nn, matriceV); //4 unrolling

			trasposta(matriceK, k_t, input->n, input->nn);

			secondaFase_U2(matriceQ, k_t, matriceSecondaTerzaFase, input->n, input->nn, fattoreScalare);//4 unrolling
		
			terzaFase(input->n, matriceSecondaTerzaFase);

			quartaFase_U2(matriceSecondaTerzaFase,matriceV,input->out,input->n,input->nn, inizioFettaFinale);
		}
		else{
			primaFase_U1(fetta, input->wq, input->bq, input->n, input->d, input->nn, matriceQ); //4 unrolling
			primaFase_U1(fetta, input->wk, input->bk, input->n, input->d, input->nn, matriceK); //4 unrolling
			primaFase_U1(fetta, input->wv, input->bv, input->n, input->d, input->nn, matriceV); //4 unrolling

			trasposta(matriceK, k_t, input->n, input->nn);

			secondaFase_U1(matriceQ, k_t, matriceSecondaTerzaFase, input->n, input->nn, fattoreScalare);//4 unrolling
		
			terzaFase_U1(input->n, matriceSecondaTerzaFase);

			quartaFase_U1(matriceSecondaTerzaFase,matriceV,input->out,input->n,input->nn, inizioFettaFinale);
		}
	}

	dealloc_matrix(fetta);
	dealloc_matrix(matriceQ);
	dealloc_matrix(matriceK);
	dealloc_matrix(matriceV);
	dealloc_matrix(matriceSecondaTerzaFase);
	dealloc_matrix(k_t);
}

void att(params *input){
	// -------------------------------------------------
	// Codificare qui l'algoritmo Attention mechanism
	// -------------------------------------------------

	#pragma omp parallel for
	for (int i = 0; i < input->ns; i++){
		calcolaTensori(input, i);	
	}
}

int main(int argc, char** argv) {

	char fname[256];
	char* dsfilename = NULL;
	char* wqfilename = NULL;
	char* wkfilename = NULL;
	char* wvfilename = NULL;
	char* bqfilename = NULL;
	char* bkfilename = NULL;
	char* bvfilename = NULL;
	clock_t t;
	float time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->ds = NULL;
	input->wq = NULL;
	input->wk = NULL;
	input->wv = NULL;
	input->bq = NULL;
	input->bk = NULL;
	input->bv = NULL;
	input->s = -1;
	input->n = -1;
	input->d = -1;
	input->ns = -1;
	
	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -ds <DS> -wq <WQ> -wk <WK> -wv <WV> -bq <bQ> -bk <bK> -bv <bV> -si <si> -n <n> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tDS: il nome del file ds2 contenente il dataset\n");
		printf("\tWQ: il nome del file ds2 contenente i pesi WQ\n");
		printf("\tWK: il nome del file ds2 contenente i pesi WK\n");
		printf("\tWV: il nome del file ds2 contenente i pesi WV\n");
		printf("\tbQ: il nome del file ds2 contenente i pesi bQ\n");
		printf("\tbK: il nome del file ds2 contenente i pesi bK\n");
		printf("\tbV: il nome del file ds2 contenente i pesi bV\n");
		printf("\tN: numero di righe del dataset\n");
		printf("\tsi: prima dimensione del tensore\n");
		printf("\tn: seconda dimensione del tensore\n");
		printf("\tnn: numero di neuroni\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-ds") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			dsfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-wq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wq file name!\n");
				exit(1);
			}
			wqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-wk") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wk file name!\n");
				exit(1);
			}
			wkfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-wv") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wv file name!\n");
				exit(1);
			}
			wvfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-bq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing bq file name!\n");
				exit(1);
			}
			bqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-bk") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing bk file name!\n");
				exit(1);
			}
			bkfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-bv") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing bv file name!\n");
				exit(1);
			}
			bvfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-si") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing si value!\n");
				exit(1);
			}
			input->s = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-n") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing n value!\n");
				exit(1);
			}
			input->n = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(dsfilename == NULL || strlen(dsfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	if(wqfilename == NULL || strlen(wqfilename) == 0){
		printf("Missing wq file name!\n");
		exit(1);
	}

	if(wkfilename == NULL || strlen(wkfilename) == 0){
		printf("Missing wk file name!\n");
		exit(1);
	}

	if(wvfilename == NULL || strlen(wvfilename) == 0){
		printf("Missing wv file name!\n");
		exit(1);
	}

	if(bqfilename == NULL || strlen(bqfilename) == 0){
		printf("Missing bq file name!\n");
		exit(1);
	}

	if(bkfilename == NULL || strlen(bkfilename) == 0){
		printf("Missing bk file name!\n");
		exit(1);
	}

	if(bvfilename == NULL || strlen(bvfilename) == 0){
		printf("Missing bv file name!\n");
		exit(1);
	}

	input->ds = load_data(dsfilename, &input->N, &input->d);

	if(input->s <= 0){
		printf("Invalid value of si parameter!\n");
		exit(1);
	}

	if(input->n <= 0 || input->N % input->n != 0){
		printf("Invalid value of n parameter!\n");
		exit(1);
	}

	input->ns = (int) ceil((double) input->N / (input->s * input->n));

	// Caricamento matrici livelli
	int n, nn;
	input->wq = load_data(wqfilename, &n, &input->nn);

	if(input->d != n){
		printf("Invalid wq size!\n");
		exit(1);
	}

	input->wk = load_data(wkfilename, &n, &nn);

	if(input->d != n || input->nn != nn){
		printf("Invalid wk size!\n");
		exit(1);
	}

	input->wv = load_data(wvfilename, &n, &nn);

	if(input->d != n || input->nn != nn){
		printf("Invalid wv size!\n");
		exit(1);
	}

	// Caricamento bias
	input->bq = load_data(bqfilename, &n, &nn);

	if(n != 1 || input->nn != nn){
		printf("Invalid bq size!\n");
		exit(1);
	}

	input->bk = load_data(bkfilename, &n, &nn);

	if(n != 1 || input->nn != nn){
		printf("Invalid bk size!\n");
		exit(1);
	}

	input->bv = load_data(bvfilename, &n, &nn);

	if(n != 1 || input->nn != nn){
		printf("Invalid bv size!\n");
		exit(1);
	}

	input->out = alloc_matrix(input->N, input->nn);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", dsfilename);
		printf("WQ file name: '%s'\n", wqfilename);
		printf("WK file name: '%s'\n", wkfilename);
		printf("WV file name: '%s'\n", wvfilename);
		printf("bQ file name: '%s'\n", bqfilename);
		printf("bK file name: '%s'\n", bkfilename);
		printf("bV file name: '%s'\n", bvfilename);
		printf("Dataset row number: %d\n", input->N);
		printf("Tensor first dimention: %d\n", input->s);
		printf("Tensor second dimention: %d\n", input->n);
		printf("Tensor third dimention: %d\n", input->d);
		printf("Dataset block number: %d\n", input->ns);
		printf("Layer neuron number: %d\n", input->nn);
	}

	//
	// Attention Mechanism
	//

	double start = omp_get_wtime(), end;
	att(input);
	end = omp_get_wtime();
	time = end - start;


	if(!input->silent)
		printf("ATT time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname, "out64omp_%d_%d_%d_%d.ds2", input->N, input->s, input->n, input->d);
	save_data(fname, input->out, input->N, input->nn);
	if(input->display){
		if(input->out == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("out: [");
			for(i=0; i<input->N; i++){
				for(j=0; j<input->nn-1; j++)
					printf("%f,", input->out[input->d*i+j]);
				printf("%f\n", input->out[input->d*i+j]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	return 0;
}
