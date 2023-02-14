%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	p		equ     4
	dim		equ		4
section .bss			; Sezione contenente dati non inizializzati
section .text			; Sezione contenente il codice macchina

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global primaFaseQ

_4fetta		equ		8  ;fetta di dimensione nxd
_4pesi		equ		12 ;dimensione dxnn
_4bias		equ		16
_4n1 			equ		20
_4dimd		equ		24
_4nn1			equ		28
_4ris1		equ		32	
				
primaFaseQ:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
		push	ebp			; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push	ebx			; salva i registri da preservare
		push	esi
		push	edi
		
		; ------------------------------------------------------------
		; PRODOTTO TRA MATRICI
		; ------------------------------------------------------------		
		
		xor		eax, eax				; eax = i = 0
fori4:	xor		ebx, ebx				; ebx = j = 0
forj4:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		xorps	xmm7, xmm7	
		xor		ecx, ecx				; ecx = k = 0

fork4:	; prendo 1 elemento dalla fetta[i*d+k] e lo metto in xmm1
		; in edi mi calcolo i*d+k
		mov		edi, eax
		imul	edi, [ebp+_4dimd]
		add		edi, ecx
		add		edi, [ebp+_4fetta]
		movss	xmm1,[edi]
       	shufps  xmm1, xmm1, 0	
		
		; prendo 4 elementi da pesi[k*nn+j] e li metto in xmm2
		; in esi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_4nn1]	
		add		edi, ebx
		add		edi, [ebp+_4pesi]
		movaps	xmm2,[edi]
	
		; prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2		

        ;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	
       	; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2
		; UNROLLING=4
		movaps	xmm2, [edi+48]	
		mulps	xmm2, xmm1
		addps	xmm7, xmm2

		add		ecx, dim					; k++
		imul	edi, [ebp+_4dimd], dim
		cmp		ecx, edi					; (k < d) ?
		jb		fork4
		
		; sommo il vettore dei bias[j]
		mov		edi, [ebp+_4bias]
		addps	xmm4, [edi+ebx]
        addps	xmm5, [edi+ebx+16]
        addps	xmm6, [edi+ebx+32]
        addps	xmm7, [edi+ebx+48]

		; copio t in ris[i][j] = ris[i*nn+j]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_4nn1]
		add		edi, ebx
		add		edi, [ebp+_4ris1]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6
		movaps	[edi+48], xmm7

		add		ebx, dim*p*4		; j+=p
		imul	edi, [ebp+_4nn1], dim
		cmp		ebx, edi		; (j < nn) ?
		jb		forj4

		add		eax, dim		;i++
		imul	edi, [ebp+_4n1], dim	
		cmp		eax, edi		; i < n ?
		jb		fori4
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante


global primaFaseT

_3fetta		equ		8  ;fetta di dimensione nxd
_3pesi		equ		12 ;dimensione dxnn
_3bias		equ		16
_3n1 			equ		20
_3dimd		equ		24
_3nn1			equ		28
_3ris1		equ		32
				
primaFaseT:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
		push	ebp			; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push	ebx			; salva i registri da preservare
		push	esi
		push	edi
		
		; ------------------------------------------------------------
		; PRODOTTO TRA MATRICI
		; ------------------------------------------------------------		
		
		xor		eax, eax				; eax = i = 0
fori3:	xor		ebx, ebx				; ebx = j = 0
forj3:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		xor		ecx, ecx				; ecx = k = 0

fork3:	; prendo 1 elemento dalla fetta[i*d+k] e lo metto in xmm1
		; in edi mi calcolo i*d+k
		mov		edi, eax
		imul	edi, [ebp+_3dimd]
		add		edi, ecx
		add		edi, [ebp+_3fetta]
		movss	xmm1,[edi]
       	shufps  xmm1, xmm1, 0	
		
		; prendo 4 elementi da pesi[k*nn+j] e li metto in xmm2
		; in esi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_3nn1]	
		add		edi, ebx
		add		edi, [ebp+_3pesi]
		movaps	xmm2,[edi]
	
		; prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2		

        ;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	
       	; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2

		add		ecx, dim					; k++
		imul	edi, [ebp+_3dimd], dim
		cmp		ecx, edi					; (k < d) ?
		jb		fork3
		
		; sommo il vettore dei bias[j]
		mov		edi, [ebp+_3bias]
		addps	xmm4, [edi+ebx]
        addps	xmm5, [edi+ebx+16]
        addps	xmm6, [edi+ebx+32]

		; copio t in ris[i][j] = ris[i*nn+j]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_3nn1]
		add		edi, ebx
		add		edi, [ebp+_3ris1]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6

		add		ebx, dim*p*3		; j+=p
		imul	edi, [ebp+_3nn1], dim
		cmp		ebx, edi		; (j < nn) ?
		jb		forj3

		add		eax, dim		;i++
		imul	edi, [ebp+_3n1], dim	
		cmp		eax, edi		; i < n ?
		jb		fori3
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante


global primaFaseD

_2fetta		equ		8  ;fetta di dimensione nxd
_2pesi		equ		12 ;dimensione dxnn
_2bias		equ		16
_2n1 			equ		20
_2dimd		equ		24
_2nn1			equ		28
_2ris1		equ		32

primaFaseD:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
		push	ebp			; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push	ebx			; salva i registri da preservare
		push	esi
		push	edi
		
		; ------------------------------------------------------------
		; PRODOTTO TRA MATRICI
		; ------------------------------------------------------------		
		
		xor		eax, eax				; eax = i = 0
fori2:	xor		ebx, ebx				; ebx = j = 0
forj2:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xor		ecx, ecx				; ecx = k = 0

fork2:	; prendo 1 elemento dalla fetta[i*d+k] e lo metto in xmm1
		; in edi mi calcolo i*d+k
		mov		edi, eax
		imul	edi, [ebp+_2dimd]
		add		edi, ecx
		add		edi, [ebp+_2fetta]
		movss	xmm1,[edi]
       	shufps  xmm1, xmm1, 0	
		
		; prendo 4 elementi da pesi[k*nn+j] e li metto in xmm2
		; in esi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_2nn1]	
		add		edi, ebx
		add		edi, [ebp+_2pesi]
		movaps	xmm2,[edi]
	
		; prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2		

        ;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	

		add		ecx, dim					; k++
		imul	edi, [ebp+_2dimd], dim
		cmp		ecx, edi					; (k < d) ?
		jb		fork2
		
		; sommo il vettore dei bias[j]
		mov		edi, [ebp+_2bias]
		addps	xmm4, [edi+ebx]
        addps	xmm5, [edi+ebx+16]

		; copio t in ris[i][j] = ris[i*nn+j]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_2nn1]
		add		edi, ebx
		add		edi, [ebp+_2ris1]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5

		add		ebx, dim*p*2		; j+=p
		imul	edi, [ebp+_2nn1], dim
		cmp		ebx, edi		; (j < nn) ?
		jb		forj2

		add		eax, dim		;i++
		imul	edi, [ebp+_2n1], dim	
		cmp		eax, edi		; i < n ?
		jb		fori2
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante


global primaFaseU

_1fetta		equ		8  ;fetta di dimensione nxd
_1pesi		equ		12 ;dimensione dxnn
_1bias		equ		16
_1n1 			equ		20
_1dimd		equ		24
_1nn1			equ		28
_1ris1		equ		32
				
primaFaseU:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
		push	ebp			; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push	ebx			; salva i registri da preservare
		push	esi
		push	edi
		
		; ------------------------------------------------------------
		; PRODOTTO TRA MATRICI
		; ------------------------------------------------------------		
		
		xor		eax, eax				; eax = i = 0
fori1:	xor		ebx, ebx				; ebx = j = 0
forj1:	; t = 0
        xorps	xmm4, xmm4	
		xor		ecx, ecx				; ecx = k = 0

fork1:	; prendo 1 elemento dalla fetta[i*d+k] e lo metto in xmm1
		; in edi mi calcolo i*d+k
		mov		edi, eax
		imul	edi, [ebp+_1dimd]
		add		edi, ecx
		add		edi, [ebp+_1fetta]
		movss	xmm1,[edi]
       	shufps  xmm1, xmm1, 0	
		
		; prendo 4 elementi da pesi[k*nn+j] e li metto in xmm2
		; in esi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_1nn1]	
		add		edi, ebx
		add		edi, [ebp+_1pesi]
		movaps	xmm2,[edi]
	
		; prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2		

		add		ecx, dim					; k++
		imul	edi, [ebp+_1dimd], dim
		cmp		ecx, edi					; (k < d) ?
		jb		fork1
		
		; sommo il vettore dei bias[j]
		mov		edi, [ebp+_1bias]
		addps	xmm4, [edi+ebx]

		; copio t in ris[i][j] = ris[i*nn+j]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_1nn1]
		add		edi, ebx
		add		edi, [ebp+_1ris1]
		movaps	[edi], xmm4

		add		ebx, dim*p		; j+=p
		imul	edi, [ebp+_1nn1], dim
		cmp		ebx, edi		; (j < nn) ?
		jb		forj1

		add		eax, dim		;i++
		imul	edi, [ebp+_1n1], dim	
		cmp		eax, edi		; i < n ?
		jb		fori1
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante


global secondaFaseQ

_q4		equ		8		
_4k		equ		12
_4ris		equ		16
_4n		equ		20
_4nn		equ		24
_4fattore equ     28

secondaFaseQ:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		movss	xmm3, [ebp+_4fattore]
		shufps	xmm3, xmm3, 0

		mov		eax, 0			; i = 0
.fori4:	mov		ebx, 0			; j = 0

.forj4:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		xorps	xmm7, xmm7	
		mov		ecx, 0			; k = 0
		
.fork4:	; prendo 1 elementi da q e lo metto in xmm1
		;q[i][k] = q[i*nn+k] = [q+i*nn+k] = [edx+edi]
		;in edi mi calcolo i*nn+k
		mov		edi, eax
		imul	edi, [ebp+_4nn]
		add		edi, ecx
		add		edi, [ebp+_q4]
		movss	xmm1, [edi]
        shufps  xmm1, xmm1, 0
		
	; prendo 4 elementi da k_trasposta
		;k_t[k][j] = k_t[k*n+j] = [k_t+k*n+j] = [edx+esi]
		;in esi mi calcolo k*n+j
		mov		edi, ecx
		imul	edi, [ebp+_4n]
		add		edi, ebx
		add		edi, [ebp+_4k]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2	

		;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	
       	; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2
		; UNROLLING=4
		movaps	xmm2, [edi+48]	
		mulps	xmm2, xmm1
		addps	xmm7, xmm2
		
        add		ecx, dim	; k++
		imul	edx, [ebp+_4nn], dim
		cmp		ecx, edx	; (k < nn) ?
		jl		.fork4

	;moltiplico xmm0 per il fattore di scala 1/sqrt(d) che sta in xmm3
		mulps	xmm4, xmm3
        mulps	xmm5, xmm3
        mulps	xmm6, xmm3
        mulps	xmm7, xmm3
	
		; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in edi mi calcolo i*n+j
		mov		edi, eax
		imul	edi, [ebp+_4n]
		add		edi, ebx
		add		edi, [ebp+_4ris]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6
		movaps	[edi+48], xmm7

		add		ebx, dim*p*4		; j+=p
		imul 	edx, [ebp+_4n], dim
		cmp		ebx, edx	; (j < n) ?
		jb		.forj4

		add		eax, dim		; i++
		imul 	edx, [ebp+_4n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori4
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global secondaFaseT

_q3		equ		8		
_3k		equ		12
_3ris		equ		16
_3n		equ		20
_3nn		equ		24
_3fattore equ     28

secondaFaseT:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		movss	xmm3, [ebp+_3fattore]
		shufps	xmm3, xmm3, 0

		mov		eax, 0			; i = 0
.fori3:	mov		ebx, 0			; j = 0

.forj3:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		mov		ecx, 0			; k = 0
		
.fork3:	; prendo 1 elementi da q e lo metto in xmm1
		;q[i][k] = q[i*nn+k] = [q+i*nn+k] = [edx+edi]
		;in edi mi calcolo i*nn+k
		mov		edi, eax
		imul	edi, [ebp+_3nn]
		add		edi, ecx
		add		edi, [ebp+_q3]
		movss	xmm1, [edi]
        shufps  xmm1, xmm1, 0
		
	; prendo 4 elementi da k_trasposta
		;k_t[k][j] = k_t[k*n+j] = [k_t+k*n+j] = [edx+esi]
		;in esi mi calcolo k*n+j
		mov		edi, ecx
		imul	edi, [ebp+_3n]
		add		edi, ebx
		add		edi, [ebp+_3k]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2	

		;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	
       	; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2
		
        add		ecx, dim	; k++
		imul	edx, [ebp+_3nn], dim
		cmp		ecx, edx	; (k < nn) ?
		jl		.fork3

	;moltiplico xmm0 per il fattore di scala 1/sqrt(d) che sta in xmm3
		mulps	xmm4, xmm3
        mulps	xmm5, xmm3
        mulps	xmm6, xmm3
	
		; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in edi mi calcolo i*n+j
		mov		edi, eax
		imul	edi, [ebp+_3n]
		add		edi, ebx
		add		edi, [ebp+_3ris]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6

		add		ebx, dim*p*3		; j+=p
		imul 	edx, [ebp+_3n], dim
		cmp		ebx, edx	; (j < n) ?
		jb		.forj3

		add		eax, dim		; i++
		imul 	edx, [ebp+_3n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori3
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global secondaFaseD

_q2		equ		8		
_2k		equ		12
_2ris		equ		16
_2n		equ		20
_2nn		equ		24
_2fattore equ     28

secondaFaseD:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		movss	xmm3, [ebp+_2fattore]
		shufps	xmm3, xmm3, 0

		mov		eax, 0			; i = 0
.fori2:	mov		ebx, 0			; j = 0

.forj2:	; t = 0
        xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		mov		ecx, 0			; k = 0
		
.fork2:	; prendo 1 elementi da q e lo metto in xmm1
		;q[i][k] = q[i*nn+k] = [q+i*nn+k] = [edx+edi]
		;in edi mi calcolo i*nn+k
		mov		edi, eax
		imul	edi, [ebp+_2nn]
		add		edi, ecx
		add		edi, [ebp+_q2]
		movss	xmm1, [edi]
        shufps  xmm1, xmm1, 0
		
	; prendo 4 elementi da k_trasposta
		;k_t[k][j] = k_t[k*n+j] = [k_t+k*n+j] = [edx+esi]
		;in esi mi calcolo k*n+j
		mov		edi, ecx
		imul	edi, [ebp+_2n]
		add		edi, ebx
		add		edi, [ebp+_2k]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2	

		;UNROLLING=2
        movaps	xmm2, [edi+16]
		mulps	xmm2, xmm1
		addps	xmm5, xmm2	
		
        add		ecx, dim	; k++
		imul	edx, [ebp+_2nn], dim
		cmp		ecx, edx	; (k < nn) ?
		jl		.fork2

	;moltiplico xmm0 per il fattore di scala 1/sqrt(d) che sta in xmm3
		mulps	xmm4, xmm3
        mulps	xmm5, xmm3
	
		; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in edi mi calcolo i*n+j
		mov		edi, eax
		imul	edi, [ebp+_2n]
		add		edi, ebx
		add		edi, [ebp+_2ris]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5

		add		ebx, dim*p*2		; j+=p
		imul 	edx, [ebp+_2n], dim
		cmp		ebx, edx	; (j < n) ?
		jb		.forj2

		add		eax, dim		; i++
		imul 	edx, [ebp+_2n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori2
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global secondaFaseU

_q1		equ		8		
_1k		equ		12
_1ris		equ		16
_1n		equ		20
_1nn		equ		24
_1fattore equ     28

secondaFaseU:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		movss	xmm3, [ebp+_1fattore]
		shufps	xmm3, xmm3, 0

		mov		eax, 0			; i = 0
.fori1:	mov		ebx, 0			; j = 0

.forj1:	; t = 0
        xorps	xmm4, xmm4	
		mov		ecx, 0			; k = 0
		
.fork1:	; prendo 1 elementi da q e lo metto in xmm1
		;q[i][k] = q[i*nn+k] = [q+i*nn+k] = [edx+edi]
		;in edi mi calcolo i*nn+k
		mov		edi, eax
		imul	edi, [ebp+_1nn]
		add		edi, ecx
		add		edi, [ebp+_q1]
		movss	xmm1, [edi]
        shufps  xmm1, xmm1, 0
		
	; prendo 4 elementi da k_trasposta
		;k_t[k][j] = k_t[k*n+j] = [k_t+k*n+j] = [edx+esi]
		;in esi mi calcolo k*n+j
		mov		edi, ecx
		imul	edi, [ebp+_1n]
		add		edi, ebx
		add		edi, [ebp+_1k]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2	
		
        add		ecx, dim	; k++
		imul	edx, [ebp+_1nn], dim
		cmp		ecx, edx	; (k < nn) ?
		jl		.fork1

	;moltiplico xmm0 per il fattore di scala 1/sqrt(d) che sta in xmm3
		mulps	xmm4, xmm3
	
		; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in edi mi calcolo i*n+j
		mov		edi, eax
		imul	edi, [ebp+_1n]
		add		edi, ebx
		add		edi, [ebp+_1ris]
		movaps	[edi], xmm4

		add		ebx, dim*p		; j+=p
		imul 	edx, [ebp+_1n], dim
		cmp		ebx, edx	; (j < n) ?
		jb		.forj1

		add		eax, dim		; i++
		imul 	edx, [ebp+_1n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori1
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global terzaFase

n3 	equ 8 ; intero indica la dimensione della matrice
MAT equ 12 ; puntatore a float, 32 bit

align 16
due: dd 2.0,2.0,2.0,2.0

align 16
zeropuntocinque: dd 0.5, 0.5, 0.5, 0.5

align 16
menouno: dd -1.0, -1.0,-1.0,-1.0

align 16
zero: dd 0.0
align 16
maschera1 : dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
		
terzaFase: 

; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi

	mov eax, [ebp + n3] ;valore di n
	mov ecx, [ebp + n3]
	imul eax,ecx ;eax_edx = ecx*eax ; n*n (dimensione del vettore)
	imul eax,4   ; in eax ho n*n*4 ovvero quanti byte occupati dalla matrice
					
	movaps xmm1,[due] ; in xmm1 ho i valori da sommare	
	movaps xmm2,[zeropuntocinque]
		
	movaps xmm3, [menouno]
	movaps xmm5,[maschera1]
	movaps xmm6,[maschera1]
	movaps xmm7,[menouno]
			
	xor ecx,ecx            ; ecx = 0 azzero l'indice
	mov ebx,[ebp+MAT]
 
ciclo:	movaps xmm0 , [ebx+ecx] ;prelevo 4 valori alla volta dalla matrice, il contatore lo incremento di 16 devo prendere 4 float di 4 byte 
		movaps xmm4,[ebx+ecx+16]
	
		andps xmm5, xmm0
		divps xmm5, xmm0

		andps xmm6,xmm4
		divps xmm6,xmm4
		
		; si nega tutto e si aggiunge uno 
		addps xmm0,xmm1 ; elementoMat+2 ad ogni elemento di xmm0 gli sommo2 
		divps xmm3, xmm0; xmm0 = xmm3/xmm0 -> -1/(elem+2) 
		movaps xmm0,xmm3 
		addps xmm0,xmm2 ; -1/(elem+2) + 0.5
		movaps xmm3, [menouno]
		mulps xmm0, xmm5 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    movaps xmm5, [maschera1]
		addps xmm0, xmm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5
		
		;UNROLLING
		addps xmm4,xmm1	
		divps xmm7, xmm4; xmm0 = xmm3/xmm0 -> -1/(elem+2) 
		movaps xmm4,xmm7 
		addps xmm4,xmm2 ; -1/(elem+2) + 0.5
		movaps xmm7, [menouno]
		mulps xmm4, xmm6 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    movaps xmm6, [maschera1]
		addps xmm4, xmm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5	
		
		movaps [ebx+ecx],xmm0 ; sposto il contenuto di xmm0 nella locazione di memoria da dove è stato estratto inizialmente 
		movaps [ebx+ecx+16],xmm4 
		
		add ecx, 32 ; mi sposto di 16 posizioni in quanto leggo 4 valori di 4 byte alla volta
		cmp ecx, eax ; (ecx) < (eax) ? 
		jb ciclo
	
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante
	

global terzaFaseNoUnroll

n3_nu 	equ 8 ; intero indica la dimensione della matrice
MAT_nu equ 12 ; puntatore a float, 32 bit

align 16
_due: dd 2.0,2.0,2.0,2.0

align 16
_zeropuntocinque: dd 0.5, 0.5, 0.5, 0.5

align 16
_menouno: dd -1.0, -1.0,-1.0,-1.0

align 16
_zero: dd 0.0
align 16
_maschera1 : dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
			dd  011111111111111111111111111111111b
		
terzaFaseNoUnroll: 

; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi

	mov eax, [ebp + n3_nu] ;valore di n
	mov ecx, [ebp + n3_nu]
	imul eax,ecx ;eax_edx = ecx*eax ; n*n (dimensione del vettore)
	imul eax,4   ; in eax ho n*n*4 ovvero quanti byte occupati dalla matrice
					
	movaps xmm1,[_due] ; in xmm1 ho i valori da sommare	
	movaps xmm2,[_zeropuntocinque]
		
	movaps xmm3, [_menouno]
	movaps xmm5,[_maschera1]
	
			
	xor ecx,ecx            ; ecx = 0 azzero l'indice
	mov ebx,[ebp+MAT_nu]
 
cicloNU:movaps xmm0 , [ebx+ecx] ;prelevo 4 valori alla volta dalla matrice, il contatore lo incremento di 16 devo prendere 4 float di 4 byte 
		
		andps xmm5, xmm0
		divps xmm5, xmm0

		; si nega tutto e si aggiunge uno 
		addps xmm0,xmm1 ; elementoMat+2 ad ogni elemento di xmm0 gli sommo2 
		divps xmm3, xmm0; xmm0 = xmm3/xmm0 -> -1/(elem+2) 
		movaps xmm0,xmm3 
		addps xmm0,xmm2 ; -1/(elem+2) + 0.5
		movaps xmm3, [_menouno]
		mulps xmm0, xmm5 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    movaps xmm5, [_maschera1]
		addps xmm0, xmm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5
		
		movaps [ebx+ecx],xmm0 ; sposto il contenuto di xmm0 nella locazione di memoria da dove è stato estratto inizialmente 
		
		add ecx, 16 ; mi sposto di 16 posizioni in quanto leggo 4 valori di 4 byte alla volta
		cmp ecx, eax ; (ecx) < (eax) ? 
		jb cicloNU
	
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante
		
	
global quartaFaseQ

_4mat		equ		8		
_4v		equ		12
_4ris		equ		16
_4n		equ		20
_4nn		equ		24
_4inizio	equ		28


quartaFaseQ:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx		; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		eax, 0			; i = 0
.fori4:	mov		ebx, 0			; j = 0
.forj4:	; t = 0
		xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		xorps	xmm7, xmm7		
		mov		ecx, 0			; k = 0
.fork4:	
		; prendo 1 elemento da mat
		;mat[i][k] = mat[i*n+k] = [mat+i*n+k] = [edx+edi]
		;in edi mi calcolo i+k*n
		mov		edi, eax
		imul	edi, [ebp+_4n]
		add		edi, ecx
		add		edi,[ebp+_4mat]
		movss	xmm1, [edi]
		shufps	xmm1, xmm1, 0

		; prendo 4 elementi da v
		;v[k][j] = v[k*nn+j] = [v+k*nn+j] = [edx+edi]
		;in edi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_4nn]
		add		edi, ebx
		add		edi, [ebp+_4v]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2

		; UNROLLING=2
		movaps	xmm2, [edi+16]	
		mulps	xmm2, xmm1
		addps	xmm5, xmm2
		; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2
		; UNROLLING=4
		movaps	xmm2, [edi+48]	
		mulps	xmm2, xmm1
		addps	xmm7, xmm2

		add		ecx, dim		; k++
		imul	edx, [ebp+_4n], dim
		cmp		ecx, edx		; (k < n) ?
		jl		.fork4

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_4nn]
		add		edi, ebx  
		add		edi, [ebp+_4ris]
		add		edi, [ebp+_4inizio]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6
		movaps	[edi+48], xmm7

		add		ebx, dim*p*4		; j+=p
		imul 	edx, [ebp+_4nn], dim
		cmp		ebx, edx	; (j < nn) ?
		jb		.forj4

		add		eax, dim		; i ++
		imul 	edx, [ebp+_4n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori4	
; ------------------------------------------------------------	
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global quartaFaseT

_3mat		equ		8		
_3v		equ		12
_3ris		equ		16
_3n		equ		20
_3nn		equ		24
_3inizio	equ		28

quartaFaseT:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx		; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		eax, 0			; i = 0
.fori3:	mov		ebx, 0			; j = 0
.forj3:	; t = 0
		xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		xorps	xmm6, xmm6	
		mov		ecx, 0			; k = 0
.fork3:	
		; prendo 1 elemento da mat
		;mat[i][k] = mat[i*n+k] = [mat+i*n+k] = [edx+edi]
		;in edi mi calcolo i+k*n
		mov		edi, eax
		imul	edi, [ebp+_3n]
		add		edi, ecx
		add		edi,[ebp+_3mat]
		movss	xmm1, [edi]
		shufps	xmm1, xmm1, 0

		; prendo 4 elementi da v
		;v[k][j] = v[k*nn+j] = [v+k*nn+j] = [edx+edi]
		;in edi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_3nn]
		add		edi, ebx
		add		edi, [ebp+_3v]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2

		; UNROLLING=2
		movaps	xmm2, [edi+16]	
		mulps	xmm2, xmm1
		addps	xmm5, xmm2
		; UNROLLING=3
		movaps	xmm2, [edi+32]	
		mulps	xmm2, xmm1
		addps	xmm6, xmm2

		add		ecx, dim		; k++
		imul	edx, [ebp+_3n], dim
		cmp		ecx, edx		; (k < n) ?
		jl		.fork3

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_3nn]
		add		edi, ebx  
		add		edi, [ebp+_3ris]
		add		edi, [ebp+_3inizio]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5
		movaps	[edi+32], xmm6

		add		ebx, dim*p*3		; j+=p
		imul 	edx, [ebp+_3nn], dim
		cmp		ebx, edx	; (j < nn) ?
		jb		.forj3

		add		eax, dim		; i ++
		imul 	edx, [ebp+_3n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori3
			
; ------------------------------------------------------------	
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante
	

global quartaFaseD

_2mat		equ		8		
_2v		equ		12
_2ris		equ		16
_2n		equ		20
_2nn		equ		24
_2inizio	equ		28

quartaFaseD:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx		; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		eax, 0			; i = 0
.fori2:	mov		ebx, 0			; j = 0
.forj2:	; t = 0
		xorps	xmm4, xmm4	
		xorps	xmm5, xmm5	
		mov		ecx, 0			; k = 0
.fork2:	
		; prendo 1 elemento da mat
		;mat[i][k] = mat[i*n+k] = [mat+i*n+k] = [edx+edi]
		;in edi mi calcolo i+k*n
		mov		edi, eax
		imul	edi, [ebp+_2n]
		add		edi, ecx
		add		edi,[ebp+_2mat]
		movss	xmm1, [edi]
		shufps	xmm1, xmm1, 0

		; prendo 4 elementi da v
		;v[k][j] = v[k*nn+j] = [v+k*nn+j] = [edx+edi]
		;in edi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_2nn]
		add		edi, ebx
		add		edi, [ebp+_2v]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2

		; UNROLLING=2
		movaps	xmm2, [edi+16]	
		mulps	xmm2, xmm1
		addps	xmm5, xmm2

		add		ecx, dim		; k++
		imul	edx, [ebp+_2n], dim
		cmp		ecx, edx		; (k < n) ?
		jl		.fork2

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_2nn]
		add		edi, ebx  
		add		edi, [ebp+_2ris]
		add		edi, [ebp+_2inizio]
		movaps	[edi], xmm4
		movaps	[edi+16], xmm5

		add		ebx, dim*p*2		; j+=p
		imul 	edx, [ebp+_2nn], dim
		cmp		ebx, edx	; (j < nn) ?
		jb		.forj2

		add		eax, dim		; i ++
		imul 	edx, [ebp+_2n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori2
			
; ------------------------------------------------------------	
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante


global quartaFaseU

_1mat		equ		8		
_1v		equ		12
_1ris		equ		16
_1n		equ		20
_1nn		equ		24
_1inizio	equ		28

quartaFaseU:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
	push	ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx		; salva i registri da preservare
	push	esi
	push	edi
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		eax, 0			; i = 0
.fori1:	mov		ebx, 0			; j = 0
.forj1:	; t = 0
		xorps	xmm4, xmm4	
		mov		ecx, 0			; k = 0
.fork1:	
		; prendo 1 elemento da mat
		;mat[i][k] = mat[i*n+k] = [mat+i*n+k] = [edx+edi]
		;in edi mi calcolo i+k*n
		mov		edi, eax
		imul	edi, [ebp+_1n]
		add		edi, ecx
		add		edi,[ebp+_1mat]
		movss	xmm1, [edi]
		shufps	xmm1, xmm1, 0

		; prendo 4 elementi da v
		;v[k][j] = v[k*nn+j] = [v+k*nn+j] = [edx+edi]
		;in edi mi calcolo k*nn+j
		mov		edi, ecx
		imul	edi, [ebp+_1nn]
		add		edi, ebx
		add		edi, [ebp+_1v]
		movaps	xmm2, [edi]	
		
		;prodotto scalare
		mulps	xmm2, xmm1
		addps	xmm4, xmm2

		add		ecx, dim		; k++
		imul	edx, [ebp+_1n], dim
		cmp		ecx, edx		; (k < n) ?
		jl		.fork1

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in edi mi calcolo i*nn+j
		mov		edi, eax
		imul	edi, [ebp+_1nn]
		add		edi, ebx  
		add		edi, [ebp+_1ris]
		add		edi, [ebp+_1inizio]
		movaps	[edi], xmm4

		add		ebx, dim*p		; j+=p
		imul 	edx, [ebp+_1nn], dim
		cmp		ebx, edx	; (j < nn) ?
		jb		.forj1

		add		eax, dim		; i ++
		imul 	edx, [ebp+_1n], dim
		cmp		eax, edx		; (i < n) ?
		jb		.fori1
			
; ------------------------------------------------------------	
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante	