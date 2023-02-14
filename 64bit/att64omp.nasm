%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
	p		equ     4
	dim		equ		8
section .bss			; Sezione contenente dati non inizializzati
	fattore:	resq	1
section .text			; Sezione contenente il codice macchina

global primaFase_U1
;RDI = fetta
;RSI = pesi
;RDX = bias
;RCX = n
;R8 = d
;R9 = nn
risU1	equ	16
				
primaFase_U1:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		r13, [rbp+risU1]
		mov		rax, 0				; rax = i = 0
fori_1:	mov		rbx, 0				; rbx = j = 0

forj_1:	; t = 0
        vxorpd	ymm4, ymm4	
		mov		r10, 0				; r10 = k = 0
		
fork_1:	; in r11 mi calcolo i*d+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		; prendo 1 elemento dalla fetta[i*d+k] e lo metto in ymm1. Rircorda che [RDI] = fetta
		vmovsd xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà
		
		; in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r9	
		add		r11, rbx
		; prendo 4 elementi da pesi[k*nn+j] e li metto in ymm2. Ricorda che [RSI] = pesi
		vmovapd	ymm2,[rsi+r11] 
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2		

		add		r10, dim					; k++
		imul	r11, r8, dim
		cmp		r10, r11					; (k < d) ?
		jb		fork_1
		
		; sommo il vettore dei bias[j]. Ricorda che [rdx] = bias
		vaddpd	ymm4, [rdx+rbx]	

		; copio t in ris[i][j] = ris[i*nn+j] = r[13+r11]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r9
		add		r11, rbx
		vmovapd	[r13+r11], ymm4	

		add		rbx, dim*p	; j+=p
		imul	r11, r9, dim
		cmp		rbx, r11			; (j < nn) ?
		jb		forj_1

		add		rax, dim			;i++	
		imul	r11, rcx, dim	
		cmp		rax, r11			; i < n ?
		jb		fori_1
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global primaFase_U2
;RDI = fetta
;RSI = pesi
;RDX = bias
;RCX = n
;R8 = d
;R9 = nn
risU2	equ	16
				
primaFase_U2:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		r13, [rbp+risU2]
		mov		rax, 0				; rax = i = 0
fori_2:	mov		rbx, 0				; rbx = j = 0

forj_2:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
			
		mov		r10, 0				; r10 = k = 0
		
fork_2:	; in r11 mi calcolo i*d+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		; prendo 1 elemento dalla fetta[i*d+k] e lo metto in ymm1. Rircorda che [RDI] = fetta
		vmovsd xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà
		
		; in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r9	
		add		r11, rbx
		; prendo 4 elementi da pesi[k*nn+j] e li metto in ymm2. Ricorda che [RSI] = pesi
		vmovapd	ymm2,[rsi+r11] 
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2		

        ;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	

		add		r10, dim					; k++
		imul	r11, r8, dim
		cmp		r10, r11					; (k < d) ?
		jb		fork_2
		
		; sommo il vettore dei bias[j]. Ricorda che [rdx] = bias
		vaddpd	ymm4, [rdx+rbx]
		vaddpd	ymm5, [rdx+rbx+32]

		; copio t in ris[i][j] = ris[i*nn+j] = r[13+r11]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r9
		add		r11, rbx
		vmovapd	[r13+r11], ymm4
		vmovapd	[r13+r11+32], ymm5

		add		rbx, dim*p*2	; j+=p
		imul	r11, r9, dim
		cmp		rbx, r11			; (j < nn) ?
		jb		forj_2

		add		rax, dim			;i++	
		imul	r11, rcx, dim	
		cmp		rax, r11			; i < n ?
		jb		fori_2
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global primaFase_U3
;RDI = fetta
;RSI = pesi
;RDX = bias
;RCX = n
;R8 = d
;R9 = nn
risU3	equ	16
				
primaFase_U3:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		r13, [rbp+risU3]
		mov		rax, 0				; rax = i = 0
fori_3:	mov		rbx, 0				; rbx = j = 0

forj_3:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
			
		mov		r10, 0				; r10 = k = 0
		
fork_3:	; in r11 mi calcolo i*d+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		; prendo 1 elemento dalla fetta[i*d+k] e lo metto in ymm1. Rircorda che [RDI] = fetta
		vmovsd xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà
		
		; in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r9	
		add		r11, rbx
		; prendo 4 elementi da pesi[k*nn+j] e li metto in ymm2. Ricorda che [RSI] = pesi
		vmovapd	ymm2,[rsi+r11] 
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2		

        ;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	
       	;UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2

		add		r10, dim					; k++
		imul	r11, r8, dim
		cmp		r10, r11					; (k < d) ?
		jb		fork_3
		
		; sommo il vettore dei bias[j]. Ricorda che [rdx] = bias
		vaddpd	ymm4, [rdx+rbx]
		vaddpd	ymm5, [rdx+rbx+32]
        vaddpd	ymm6, [rdx+rbx+64]
        
		
		; copio t in ris[i][j] = ris[i*nn+j] = r[13+r11]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r9
		add		r11, rbx
		vmovapd	[r13+r11], ymm4
		vmovapd	[r13+r11+32], ymm5
		vmovapd	[r13+r11+64], ymm6
		
		
		add		rbx, dim*p*3; j+=p
		imul	r11, r9, dim
		cmp		rbx, r11			; (j < nn) ?
		jb		forj_3

		add		rax, dim			;i++	
		imul	r11, rcx, dim	
		cmp		rax, r11			; i < n ?
		jb		fori_3
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global primaFase_U4
;RDI = fetta
;RSI = pesi
;RDX = bias
;RCX = n
;R8 = d
;R9 = nn
risU4	equ	16
	
primaFase_U4:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		r13, [rbp+risU4]
		mov		rax, 0				; rax = i = 0
fori_4:	mov		rbx, 0				; rbx = j = 0

forj_4:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
		vxorpd	ymm7, ymm7
			
		mov		r10, 0				; r10 = k = 0
		
fork_4:	; in r11 mi calcolo i*d+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		; prendo 1 elemento dalla fetta[i*d+k] e lo metto in ymm1. Rircorda che [RDI] = fetta
		vmovsd xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà
		
		; in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r9	
		add		r11, rbx
		; prendo 4 elementi da pesi[k*nn+j] e li metto in ymm2. Ricorda che [RSI] = pesi
		vmovapd	ymm2,[rsi+r11] 
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2		

        ;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	
       	;UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2
		;UNROLLING=4
		vmovapd	ymm2, [rsi+r11+96]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm7, ymm2
		

		add		r10, dim					; k++
		imul	r11, r8, dim
		cmp		r10, r11					; (k < d) ?
		jb		fork_4
		
		; sommo il vettore dei bias[j]. Ricorda che [rdx] = bias
		vaddpd	ymm4, [rdx+rbx]
		vaddpd	ymm5, [rdx+rbx+32]
        vaddpd	ymm6, [rdx+rbx+64]
        vaddpd	ymm7, [rdx+rbx+96]

		; copio t in ris[i][j] = ris[i*nn+j] = r[13+r11]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r9
		add		r11, rbx
		vmovapd	[r13+r11], ymm4
		vmovapd	[r13+r11+32], ymm5
		vmovapd	[r13+r11+64], ymm6
		vmovapd	[r13+r11+96], ymm7

		add		rbx, dim*p*4	; j+=p
		imul	r11, r9, dim
		cmp		rbx, r11			; (j < nn) ?
		jb		forj_4

		add		rax, dim			;i++	
		imul	r11, rcx, dim	
		cmp		rax, r11			; i < n ?
		jb		fori_4
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global secondaFase_U1
;RDI = q			
;RSI = k		
;RDX = ris		
;RCX = n		
;R8 = nn	
;xmm0 = fattore

secondaFase_U1:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
			
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------
		movups	[fattore], xmm0
		vbroadcastsd ymm0, [fattore]

		mov		rax, 0			; i = 0
.fori:	mov		rbx, 0			; j = 0

.forj:	; t = 0
        vxorpd	ymm4, ymm4	
		mov		r10, 0			; k = 0
		
.fork:	;prendo 1 elemento da q
		;in r11 mi calcolo i*nn+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		vmovsd 	xmm2, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm2	                ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm1 nella seconda metà
		
		;prendo 4 elementi da k_trasposta
		;in r11 mi calcolo k*n+j
		mov		r11, r10
		imul	r11, rcx
		add		r11, rbx
		vmovapd	ymm2, [rsi+r11]	
		;prodotto scalare
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2	
	
        add		r10, dim	; k++
		imul	r11, r8, dim
		cmp		r10, r11	; (k < nn) ?
		jl		.fork

	;moltiplico per il fattore di scala 1/sqrt(d) che sta in ymm0
		vmulpd	ymm4, ymm0
	
	; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in r11 mi calcolo i*n+j
		mov		r11, rax
		imul	r11, rcx
		add		r11, rbx
		vmovapd	[rdx+r11], ymm4

		add		rbx, dim*p		; j+=p
		imul 	r11, rcx, dim
		cmp		rbx, r11	; (j < n) ?
		jb		.forj

		add		rax, dim		; i++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global secondaFase_U2
;RDI = q			
;RSI = k		
;RDX = ris		
;RCX = n		
;R8 = nn	
;xmm0 = fattore

secondaFase_U2:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
			
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------
		movups	[fattore], xmm0
		vbroadcastsd ymm0, [fattore]

		mov		rax, 0			; i = 0
.fori2:	mov		rbx, 0			; j = 0

.forj2:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		
		mov		r10, 0			; k = 0
		
.fork2:	;prendo 1 elemento da q
		;in r11 mi calcolo i*nn+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		vmovsd 	xmm2, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm2	                ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm1 nella seconda metà
		
		;prendo 4 elementi da k_trasposta
		;in r11 mi calcolo k*n+j
		mov		r11, r10
		imul	r11, rcx
		add		r11, rbx
		vmovapd	ymm2, [rsi+r11]	
		;prodotto scalare
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2	

		;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	
       	
	
        add		r10, dim	; k++
		imul	r11, r8, dim
		cmp		r10, r11	; (k < nn) ?
		jl		.fork2

	;moltiplico per il fattore di scala 1/sqrt(d) che sta in ymm0
		vmulpd	ymm4, ymm0
        vmulpd	ymm5, ymm0
        
	
	; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in r11 mi calcolo i*n+j
		mov		r11, rax
		imul	r11, rcx
		add		r11, rbx
		vmovapd	[rdx+r11], ymm4
		vmovapd	[rdx+r11+32], ymm5
		

		add		rbx, dim*p*2		; j+=p
		imul 	r11, rcx, dim
		cmp		rbx, r11	; (j < n) ?
		jb		.forj2

		add		rax, dim		; i++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori2
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global secondaFase_U3
;RDI = q			
;RSI = k		
;RDX = ris		
;RCX = n		
;R8 = nn	
;xmm0 = fattore

secondaFase_U3:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
			
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------
		movups	[fattore], xmm0
		vbroadcastsd ymm0, [fattore]

		mov		rax, 0			; i = 0
.fori3:	mov		rbx, 0			; j = 0

.forj3:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
	
		
		mov		r10, 0			; k = 0
		
.fork3:	;prendo 1 elemento da q
		;in r11 mi calcolo i*nn+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		vmovsd 	xmm2, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm2	                ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm1 nella seconda metà
		
		;prendo 4 elementi da k_trasposta
		;in r11 mi calcolo k*n+j
		mov		r11, r10
		imul	r11, rcx
		add		r11, rbx
		vmovapd	ymm2, [rsi+r11]	
		;prodotto scalare
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2	

		;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	
       	; UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2
		
	
        add		r10, dim	; k++
		imul	r11, r8, dim
		cmp		r10, r11	; (k < nn) ?
		jl		.fork3

	;moltiplico per il fattore di scala 1/sqrt(d) che sta in ymm0
		vmulpd	ymm4, ymm0
        vmulpd	ymm5, ymm0
        vmulpd	ymm6, ymm0
     
		
	
	; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in r11 mi calcolo i*n+j
		mov		r11, rax
		imul	r11, rcx
		add		r11, rbx
		vmovapd	[rdx+r11], ymm4
		vmovapd	[rdx+r11+32], ymm5
		vmovapd	[rdx+r11+64], ymm6
		

		add		rbx, dim*p*3		; j+=p
		imul 	r11, rcx, dim
		cmp		rbx, r11	; (j < n) ?
		jb		.forj3

		add		rax, dim		; i++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori3
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global secondaFase_U4
;RDI = q			
;RSI = k		
;RDX = ris		
;RCX = n		
;R8 = nn	
;xmm0 = fattore

secondaFase_U4:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
			
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------
		movups	[fattore], xmm0
		vbroadcastsd ymm0, [fattore]

		mov		rax, 0			; i = 0
.fori4:	mov		rbx, 0			; j = 0

.forj4:	; t = 0
        vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
		vxorpd	ymm7, ymm7	
		
		mov		r10, 0			; k = 0
		
.fork4:	;prendo 1 elemento da q
		;in r11 mi calcolo i*nn+k
		mov		r11, rax
		imul	r11, r8
		add		r11, r10
		vmovsd 	xmm2, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm2	                ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm1 nella seconda metà
		
		;prendo 4 elementi da k_trasposta
		;in r11 mi calcolo k*n+j
		mov		r11, r10
		imul	r11, rcx
		add		r11, rbx
		vmovapd	ymm2, [rsi+r11]	
		;prodotto scalare
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2	

		;UNROLLING=2
        vmovapd	ymm2, [rsi+r11+32]
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2	
       	; UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2
		; UNROLLING=4
		vmovapd	ymm2, [rsi+r11+96]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm7, ymm2
		
	
        add		r10, dim	; k++
		imul	r11, r8, dim
		cmp		r10, r11	; (k < nn) ?
		jl		.fork4

	;moltiplico per il fattore di scala 1/sqrt(d) che sta in ymm0
		vmulpd	ymm4, ymm0
        vmulpd	ymm5, ymm0
        vmulpd	ymm6, ymm0
        vmulpd	ymm7, ymm0
		
	
	; copio t in ris[i][j] = ris[i*n+j] = [ris+i*n+j] = [edx+edi]
		;in r11 mi calcolo i*n+j
		mov		r11, rax
		imul	r11, rcx
		add		r11, rbx
		vmovapd	[rdx+r11], ymm4
		vmovapd	[rdx+r11+32], ymm5
		vmovapd	[rdx+r11+64], ymm6
		vmovapd	[rdx+r11+96], ymm7
		

		add		rbx, dim*p*4		; j+=p
		imul 	r11, rcx, dim
		cmp		rbx, r11	; (j < n) ?
		jb		.forj4

		add		rax, dim		; i++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori4
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global terzaFase
;RDI = n
;RSI = MAT
align 32
due: dq 2.0,2.0,2.0,2.0

align 32
zeropuntocinque: dq 0.5, 0.5, 0.5, 0.5

align 32
menouno: dq -1.0, -1.0,-1.0,-1.0

align 32
zero: dq 0.0

align 32
maschera1 : dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			
terzaFase: 
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	mov rax, rdi ;valore di rdi=n
	mov rcx, rdi ;in rsi abbiamo il puntatore alla matrice
	imul rax,rcx ;n*n (dimensione della matrice)
	imul rax, 8  ;in rax ho n*n*8 (double) ovvero quanti byte occupati dalla matrice

	vmovapd ymm1,[due] ; in xmm1 ho i valori da sommare	
	vmovapd ymm2,[zeropuntocinque]
	
	vmovapd ymm3, [menouno]
	vmovapd ymm5,[maschera1]
	vmovapd ymm6,[maschera1]
	vmovapd ymm7,[menouno]
			
	xor rcx,rcx            ; ecx = 0 azzero l'indice
 
ciclo:	vmovapd ymm0 , [rsi+rcx] ;prelevo 4 valori alla volta dalla matrice, il contatore lo incremento di 16 devo prendere 4 float di 4 byte 
	    vmovapd ymm4,[rsi+rcx+32]

		vandpd ymm5, ymm0
		vdivpd ymm5, ymm0

		vandpd ymm6,ymm4
		vdivpd ymm6,ymm4
		
		; si nega tutto e si aggiunge uno 
		vaddpd ymm0, ymm1 ;ad ogni elemento di ymm0 gli sommo 2 
		vdivpd ymm3, ymm0 ;ymm0 = ymm3/ymm0 -> -1/(elem+2) 
		vmovapd ymm0,ymm3 
		vaddpd ymm0, ymm2 ; -1/(elem+2) + 0.5
		vmovapd ymm3, [menouno]
		vmulpd ymm0, ymm5 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    vmovapd ymm5, [maschera1]
		vaddpd ymm0, ymm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5

		;UNROLLING
		vaddpd ymm4,ymm1 ; elementoMat+2 ad ogni elemento di xmm0 gli sommo2 
		vdivpd ymm7, ymm4; xmm0 = xmm3/xmm0 -> -1/(elem+2) 
		vmovapd ymm4,ymm7 
		vaddpd ymm4,ymm2 ; -1/(elem+2) + 0.5
		vmovapd ymm7, [menouno]
		vmulpd ymm4, ymm6 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    vmovapd ymm6, [maschera1]
		vaddpd ymm4, ymm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5
		
		vmovapd [rsi+rcx],ymm0 ; sposto il contenuto di xmm0 nella locazione di memoria da dove è stato estratto inizialmente 
		vmovapd [rsi+rcx+32],ymm4 
		
		add rcx, 64 ; mi sposto di 16 posizioni in quanto leggo 4 valori di 4 byte alla volta
		cmp rcx, rax ; (ecx) < (eax) ? 
		jb ciclo

; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global terzaFase_U1
;RDI = n
;RSI = MAT
align 32
due_1: dq 2.0,2.0,2.0,2.0

align 32
zeropuntocinque_1: dq 0.5, 0.5, 0.5, 0.5

align 32
menouno_1: dq -1.0, -1.0,-1.0,-1.0

align 32
zero_1: dq 0.0

align 32
maschera1_1 : dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			dq  0111111111111111111111111111111111111111111111111111111111111111b
			
terzaFase_U1: 
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	mov rax, rdi ;valore di rdi=n
	mov rcx, rdi ;in rsi abbiamo il puntatore alla matrice
	imul rax,rcx ;n*n (dimensione della matrice)
	imul rax, 8  ;in rax ho n*n*8 (double) ovvero quanti byte occupati dalla matrice

	vmovapd ymm1,[due_1] ; in xmm1 ho i valori da sommare	
	vmovapd ymm2,[zeropuntocinque_1]
	
	vmovapd ymm3, [menouno_1]
	vmovapd ymm5,[maschera1_1]
	
			
	xor rcx,rcx            ; ecx = 0 azzero l'indice
 
ciclo_1:	vmovapd ymm0 , [rsi+rcx] ;prelevo 4 valori alla volta dalla matrice, il contatore lo incremento di 16 devo prendere 4 float di 4 byte 
	    vmovapd ymm4,[rsi+rcx+32]

		vandpd ymm5, ymm0
		vdivpd ymm5, ymm0

		vandpd ymm6,ymm4
		vdivpd ymm6,ymm4
		
		; si nega tutto e si aggiunge uno 
		vaddpd ymm0, ymm1 ;ad ogni elemento di ymm0 gli sommo 2 
		vdivpd ymm3, ymm0 ;ymm0 = ymm3/ymm0 -> -1/(elem+2) 
		vmovapd ymm0,ymm3 
		vaddpd ymm0, ymm2 ; -1/(elem+2) + 0.5
		vmovapd ymm3, [menouno_1]
		vmulpd ymm0, ymm5 ; sgn(elem) *(-1/(elem+2) + 0.5)  devo trovare il modo di mettere il segno in xmm7
	    vmovapd ymm5, [maschera1_1]
		vaddpd ymm0, ymm2 ; sgn(elem) *(-1/(elem+2) + 0.5) + 0.5

		
		
		vmovapd [rsi+rcx],ymm0 ; sposto il contenuto di xmm0 nella locazione di memoria da dove è stato estratto inizialmente 
		
		
		add rcx, 32 ; mi sposto di 16 posizioni in quanto leggo 4 valori di 4 byte alla volta
		cmp rcx, rax ; (ecx) < (eax) ? 
		jb ciclo_1

; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global quartaFase_U1
;RDI = mat	
;RSI = v		
;RDX = ris		
;RCX = n		
;R8 = nn
;r9 = inizio		

quartaFase_U1:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		rax, 0			; i = 0
.fori:	mov		rbx, 0			; j = 0
.forj:	; t = 0
		vxorpd	ymm4, ymm4	
		
		mov		r10, 0			; k = 0
.fork:	; prendo 1 elemento da mat = [RDI]
		;in r11 mi calcolo i+k*n
		mov		r11, rax
		imul	r11, rcx
		add		r11, r10
		vmovsd 	xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà

		;in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r8
		add		r11, rbx
		;prendo 4 elementi da v
		vmovapd	ymm2, [rsi+r11]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2

		

		add		r10, dim		; k++
		imul	r11, rcx, dim
		cmp		r10, r11		; (k < n) ?
		jl		.fork

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r8
		add		r11, rbx
		add		r11, r9
		add		r11, rdx
		vmovapd	[r11], ymm4
		

		add		rbx, dim*p		; j+=p
		imul 	r11, r8, dim
		cmp		rbx, r11				; (j < nn) ?
		jb		.forj

		add		rax, dim		; i ++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global quartaFase_U2
;RDI = mat	
;RSI = v		
;RDX = ris		
;RCX = n		
;R8 = nn
;r9 = inizio		

quartaFase_U2:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		rax, 0			; i = 0
.fori42:	mov		rbx, 0			; j = 0
.forj42:	; t = 0
		vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		
		mov		r10, 0			; k = 0
.fork42:	; prendo 1 elemento da mat = [RDI]
		;in r11 mi calcolo i+k*n
		mov		r11, rax
		imul	r11, rcx
		add		r11, r10
		vmovsd 	xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà

		;in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r8
		add		r11, rbx
		;prendo 4 elementi da v
		vmovapd	ymm2, [rsi+r11]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2

		; UNROLLING=2
		vmovapd	ymm2, [rsi+r11+32]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2
		

		add		r10, dim		; k++
		imul	r11, rcx, dim
		cmp		r10, r11		; (k < n) ?
		jl		.fork42

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r8
		add		r11, rbx
		add		r11, r9
		add		r11, rdx
		vmovapd	[r11], ymm4
		vmovapd	[r11+32], ymm5
		

		add		rbx, dim*p*2	; j+=p
		imul 	r11, r8, dim
		cmp		rbx, r11				; (j < nn) ?
		jb		.forj42

		add		rax, dim		; i ++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori42
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global quartaFase_U3
;RDI = mat	
;RSI = v		
;RDX = ris		
;RCX = n		
;R8 = nn
;r9 = inizio		

quartaFase_U3:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		rax, 0			; i = 0
.fori43:	mov		rbx, 0			; j = 0
.forj43:	; t = 0
		vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
		
		mov		r10, 0			; k = 0
.fork43:	; prendo 1 elemento da mat = [RDI]
		;in r11 mi calcolo i+k*n
		mov		r11, rax
		imul	r11, rcx
		add		r11, r10
		vmovsd 	xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà

		;in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r8
		add		r11, rbx
		;prendo 4 elementi da v
		vmovapd	ymm2, [rsi+r11]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2

		; UNROLLING=2
		vmovapd	ymm2, [rsi+r11+32]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2
		; UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2
		

		add		r10, dim		; k++
		imul	r11, rcx, dim
		cmp		r10, r11		; (k < n) ?
		jl		.fork43

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r8
		add		r11, rbx
		add		r11, r9
		add		r11, rdx
		vmovapd	[r11], ymm4
		vmovapd	[r11+32], ymm5
		vmovapd	[r11+64], ymm6

		add		rbx, dim*p*3	; j+=p
		imul 	r11, r8, dim
		cmp		rbx, r11				; (j < nn) ?
		jb		.forj43

		add		rax, dim		; i ++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori43
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante


global quartaFase_U4
;RDI = mat	
;RSI = v		
;RDX = ris		
;RCX = n		
;R8 = nn
;r9 = inizio		

quartaFase_U4:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push	rbp				; salva il Base Pointer
mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
pushaq					; salva i registri generali
		
	; ------------------------------------------------------------
	; PRODOTTO TRA MATRICI
	; ------------------------------------------------------------		
		mov		rax, 0			; i = 0
.fori44:	mov		rbx, 0			; j = 0
.forj44:	; t = 0
		vxorpd	ymm4, ymm4	
		vxorpd	ymm5, ymm5	
		vxorpd	ymm6, ymm6	
		vxorpd	ymm7, ymm7	
		
		mov		r10, 0			; k = 0
.fork44:	; prendo 1 elemento da mat = [RDI]
		;in r11 mi calcolo i+k*n
		mov		r11, rax
		imul	r11, rcx
		add		r11, r10
		vmovsd 	xmm0, [rdi+r11]   			;carica il valore double in xmm0
		vmovddup xmm1, xmm0                 ;duplica il valore in xmm1
		vinsertf128 ymm1, ymm1, xmm1, 1     ;inserisce xmm1 in ymm0 nella seconda metà

		;in r11 mi calcolo k*nn+j
		mov		r11, r10
		imul	r11, r8
		add		r11, rbx
		;prendo 4 elementi da v
		vmovapd	ymm2, [rsi+r11]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm4, ymm2

		; UNROLLING=2
		vmovapd	ymm2, [rsi+r11+32]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm5, ymm2
		; UNROLLING=3
		vmovapd	ymm2, [rsi+r11+64]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm6, ymm2
		; UNROLLING=4
		vmovapd	ymm2, [rsi+r11+96]	
		vmulpd	ymm2, ymm1
		vaddpd	ymm7, ymm2
		

		add		r10, dim		; k++
		imul	r11, rcx, dim
		cmp		r10, r11		; (k < n) ?
		jl		.fork44

		; copio t in ris[i][j] = ris[i*nn+j] = [ris+i*nn+j] = [edx+edi]
		;in r11 mi calcolo i*nn+j
		mov		r11, rax
		imul	r11, r8
		add		r11, rbx
		add		r11, r9
		add		r11, rdx
		vmovapd	[r11], ymm4
		vmovapd	[r11+32], ymm5
		vmovapd	[r11+64], ymm6
		vmovapd	[r11+96], ymm7
		

		add		rbx, dim*p*4	; j+=p
		imul 	r11, r8, dim
		cmp		rbx, r11				; (j < nn) ?
		jb		.forj44

		add		rax, dim		; i ++
		imul 	r11, rcx, dim
		cmp		rax, r11		; (i < n) ?
		jb		.fori44
			
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
popaq				; ripristina i registri generali
mov		rsp, rbp	; ripristina lo Stack Pointer
pop		rbp		; ripristina il Base Pointer
ret				; torna alla funzione C chiamante