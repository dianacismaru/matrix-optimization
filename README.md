**Copyright (c) 2024, Cismaru-Diana Iuliana (331CA / 2023 - 2024)**

# Optimizarea operatiilor cu matrici
### ASC - Tema 3

## Descriere
> Tema consta in implementarea a trei variante de rezolvare pentru expresia matriceala
> C = (A<sup>t</sup> x B + B x A) x B<sup>t</sup>, unde toate matricele sunt
> patratice, de dimensiune N x N, iar matricea A este **superior triunghiulara**.

In toate variantele, matricele sunt alocate liniar in memorie. Rularea a fost 
efectuata pe partitia `haswell`, folosind `sbatch`, iar compilarea a fost facuta 
folosind `gcc-8.5.0`.

## 1. Implementarea BLAS
  In aceasta parte am folosit functii din biblioteca BLAS (Basic Linear Algebra 
Subprograms) pentru a efectua diverse operatii cu matrici, dar si cu vectori.

  Functiile folosite sunt:
  - `dgemm` - pentru inmultirea a doua matrici uzuale
  - `dtrmm` - pentru inmultirea unei matrici cu o matrice triunghiulara
  - `dcopy` - pentru copierea continutului unui vector in altul
  - `daxpy` - pentru adunarea a doi vectori

Pentru inceput, am realizat inmultirea matricelor **A<sup>t</sup>** si **B**, pastrand rezultatul in
variabila *A<sup>t</sup>B*. Am folosit functia `dtrmm`, deoarece matricea **A** este triunghiulara.
Avand in vedere ca in calculul dorit matricea A apare sub forma ei transpusa, am 
folosit parametrul **CblasTrans**. De asemenea, am folosit parametrul **CblasUpper**
pentru a specifica ca matricea **A** este superior triunghiulara si parametrul
**CblasLeft** pentru a specifica ca matricea **A** reprezinta partea stanga a 
inmultirii. `dtrmm` face ca rezultatul inmultirii sa fie stocat in matricea din
dreapta, motiv pentru care am alocat o copie a matricei **B** in **A<sup>t</sup>B**, folosind
functia `dcopy`.

Analog, pentru realizarea inmultirii **B x A**, am folosit functia `dtrmm`, dar fara
a mai fi nevoie sa transpun vreo matrice, iar parametrul **CblasRight** indica ca
matricea **A** reprezinta partea dreapta a inmultirii. Rezultatul este stocat in **BA**.

Avand in vedere ca matricile rezultate **A<sup>t</sup>B** si **BA** sunt stocate liniar in memorie,
am putut utiliza functia `daxpy` pentru a le aduna. Rezultatul se pastreaza in al
doilea termen al adunarii, adica **BA**.

In final, am inmultit rezultatul anterior cu matricea **B** transpusa, folosind
functia `dgemm`. Rezultatul final este stocat in matricea **C**.

## 2. Implementarea NEOPTimizata
  Varianta neoptimizata presupune un calcul brut al fiecarei operatii, fara a se
folosii artificii avansate pentru a eficientiza accesul la memorie.

  Pentru inceput, transpunerea matricei triunghiulare **A** se face in functia
`transpose_upper_triangular`, in care evitam operatiile inutile prin alocarea
matricei **A<sup>t</sup>** prin `calloc` si prin copierea doar a elementelor de deasupra diagonalei
din matricea A.

  Functia de adunare a doua matrici `add_matrices` presupune suma elementelor
corespondente din cele doua matrici, stocand rezultatul intr-o noua matrice.

  Inmultirile matricilor se realizeaza in mod clasic, folosind 3 for-uri imbricate,
in care se calculeaza suma produselor elementelor de pe fiecare linie si coloana, 
dupa formula `c[i][j] = a[i][k] * b[k][j]`, cu mentiunea ca matricele sunt alocate
liniar in memorie, formula devenind `c[i * N + j] = a[i * N + k] * b[k * N + j]`.
Totusi, cele 3 inmultiri necesare au obtinut implementari distincte, fiind
diferentiate de tipul matricilor implicate. Astfel, avem urmatoarele cazuri:

* **A<sup>t</sup> x B** utilizeaza functia `multiply_lower_triangular`, deoarece, prin transpunere,
matricea **A<sup>t</sup>** devine una inferior triunghiulara. Astfel, se evita inmultirea cu 
elementele de sub diagonala principala a matricei A, care sunt 0, prin modificarea
limitelor for-ului interior, in care **k** va merge doar pana la **i** inclusiv

* **B x A** utilizeaza functia `multiply_upper_triangular`, analog cu cazul anterior,
dar **k** va merge doar pana la **j** inclusiv

* **(A<sup>t</sup> x B + B x A) x B<sup>t</sup>** utilizeaza functia `multiply_with_transpose`, care
utilizeaza varianta clasica, dar in care matricea **B** este transpusa pe loc:
deci, in loc sa acceseze elementele matricei **B** ca `B[k][j]`, se acceseaza ca
`B[j][k]`.

## 3. Implementarea OPTimizata
- descrierea implementarii

Ideea din spatele tuturor optimizarilor porneste de la secventa de cod prezentata
in laboraorul 9, prin care se eficientizeaza accesul la memorie prin intermediul
pointerilor. Astfel, se evita accesul vectorial prin dereferentiere:

```
for(i = 0; i < N; i++){
	double *orig_pa = &a[i][0];
	for(j = 0; j < N; j++){
		double *pa = orig_pa;
		double *pb = &b[0][j];
		register double sum = 0;
		for(k = 0; k < N; k++){
			sum += *pa * *pb;
			pa++;
			pb += N;
		}
		c[i][j] = sum;
	}
}
```

Implementarea mea adapteaza codul de mai sus la matrici ce sunt alocate liniar in 
memorie si mentine un pointer chiar si pentru matricea rezultat. 
Deoarece pointerii sunt folositi in mod repetat in calculele din for-uri, acestia
sunt declarati ca **register**, pentru a se trimite un semnal compilatorului sa ii
tina in registrele procesorului pentru un acces mult mai rapid. Totusi, pentru a nu
epuiza registrele de pe CPU, declararile au fost realizate la inceputul functiilor,
nu in interiorul buclelor. Aceasta modificare a adus imbunatatiri semnificative in
timpul de executie al programului.

De asemenea, in majoritatea functiilor implementate, pentru a evita repetarea
calculului `i * N`, am folosit un registru auxiliar `index`, la care am adunat
`N` la fiecare iteratie a for-ului exterior.

## Memoria
In toate cele 3 variante, am facut alocari cu `calloc`, pentru a asigura ca toate
elementele matricilor sunt initializate cu 0. De asemenea, am eliberat memoria
folosita pentru matrici cu `free`, pentru a evita memory leak-urile.

Pentru a verifica ca nu exista probleme de acces la memorie, s-au rulat cele 3
executabile cu urmatoarea comanda:
```
valgrind --tool=memcheck --leak-check=full ./tema3_varianta ../input/input_valgrind > ../memory/varianta.memory 2>&1
```
Rezultatele comenzii sunt disponibile in directorul `memory/`, evidentiindu-se
faptul ca nu exista probleme de acces la memorie, nici memory leak-uri.

## Cache-ul
Pentru a analiza informatii legate de accesul la cache, am folosit urmatoarea comanda
pe toate cele 3 executabile:
```
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema3_varianta ../input/input_valgrind > ../cache/varianta.cache 2>&1
```

#### Rezultatele obtinute cu `cachegrind`
- `I refs` - Intruction references
	* tema3_blas: 248,609,027
	* tema3_neopt: 3,626,867,963
	* tema3_opt_m: 1,765,290,607
- `I1 misses` - L1 instruction cache misses
	* tema3_blas: 23,240
	* tema3_neopt: 1,672
	* tema3_opt_m: 1,643
- `LLi misses` - L2 instruction cache misses
	* tema3_blas: 3,526
	* tema3_neopt: 1,592
	* tema3_opt_m: 1,568
- `D refs` - Data references
	* tema3_blas: 94,660,230
	* tema3_neopt: 1,944,113,346
	* tema3_opt_m: 722,999,257
- `D1 misses` - L1 data cache misses
	* tema3_blas: 1,754,629
	* tema3_neopt: 52,454,225
	* tema3_opt_m: 52,454,401
- `LLd misses` - L2 data cache misses
	* tema3_blas: 140,181
	* tema3_neopt: 133,278
	* tema3_opt_m: 133,273
- `LL refs` - Last level cache references
	* tema3_blas: 1,777,869
	* tema3_neopt: 52,455,897
	* tema3_opt_m: 52,456,044
- `LL misses` - Last level cache misses
	* tema3_blas: 143,707
	* tema3_neopt: 134,870
	* tema3_opt_m: 134,841
- `Branches` - Branch instructions
	* tema3_blas: 5,876,858
	* tema3_neopt: 133,573,337
	* tema3_opt_m: 133,573,272
- `Mispredicts` - Branch mispredictions
	* tema3_blas: 67,881
	* tema3_neopt: 503,137
	* tema3_opt_m: 503,108

#### Interpretarea valorilor
> Varianta **BLAS** este cea mai eficienta in termeni de referinte la instructiuni si
> accesari la date, obtinand cele mai mici valori pentru majoritatea metricilor, cu exceptia
> miss-urilor de cache
 
> Din **I refs** si **D refs** se observa ca varianta neoptimizata executa de departe cele 
> mai multe instructiuni si acceseaza cele mai multe date. Pe de alta parte, varianta
> optimizata aduce imbunatariri semnificative in acest sens, avand un numar de referinte
> la instructiuni si de accesari la date de aproape 3 ori mai mic, datorita faptului ca
> am utilizat pointeri si am folosit registre, nemaiavand nevoie de accesari inutile la 
> memorie. In timp ce aceste variante aduc numere de ordinul miliardelor, varianta BLAS are 
> valori de ordinul sutelor de milioane, ceea ce demonstreaza eficienta folosirii functiilor 
> BLAS.

> Pentru toate celelalte metrici, varianta optimizata are cam aceleasi valori cu varianta
> neoptimizata, invingatoare iesind, din nou, varianta BLAS.

## Analiza comparativa a performantei
grafice!
	- legenda, unitati de masura
	- timpi de rulare
400
800
1000
1200
1400

## Resurse:
* [BLAS Quick Reference Guide](https://www.netlib.org/blas/blasqr.pdf)
* [Laborator 9 - Tehnici de Optimizare de Cod - Inmultirea Matricelor](https://ocw.cs.pub.ro/courses/asc/laboratoare/09)
* [Tema 3 - ASC](https://ocw.cs.pub.ro/courses/asc/teme/tema3)