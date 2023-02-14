
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
