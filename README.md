# LLA - Little Linear Algebra

Biblioteca pequena em C para operações de álgebra linear com vetores e matrizes de números reais (`double`).

O objetivo é manter o código simples, sem dependências externas além da biblioteca padrão do C, e fácil de incorporar em outros projetos.

Por ser um projeto recreativo, a biblioteca prioriza clareza sobre desempenho bruto. A partir da versão 0.3.0 as matrizes são guardadas em um **único bloco contíguo de doubles** (formato column-major), uma alocação por matriz, e as operações acessam o bloco diretamente por ponteiros — a validação de índices acontece só na API pública (`matrix_get`/`matrix_set`). A seção [Limites de dimensão](#limites-de-dimensão) traz os números concretos: em resumo, matrizes de até algumas centenas de linhas são instantâneas, e as operações cúbicas (produto, LU, determinante, inversa) ficam viáveis até n ≈ 2000.

## Estrutura

```
.
├── Makefile
└── src
    ├── vector.h / vector.c   # tipo Vector e operações sobre vetores
    ├── matrix.h / matrix.c   # tipo Matrix e operações sobre matrizes
    ├── lu.h / lu.c           # decomposição LU com pivoteamento parcial
    ├── svd.h / svd.c         # decomposição em valores singulares
    ├── lsystem.h / lsystem.c # resolução de sistemas lineares (LU e mínimos quadrados via SVD)
    ├── examples               # exemplos de uso (ver examples/README.md)
    ├── main.c                # testes
    ├── bech.c                # benchmark de performance
    └── lim.c                 # limites de dimensão
```

## Compilação

Requer `gcc` (ou outro compilador C99 compatível) e `make`.

```sh
make build      # gera o executável ./al
make run        # executa os testes
make bench      # gera o executável ./bech
make bench-run  # compila e executa o benchmark
make lim        # gera o executável ./lim
make lim-run    # compila e executa a busca de limites
make libs       # gera liblla.a e liblla.so
make static     # gera apenas a biblioteca estática liblla.a
make shared     # gera apenas a biblioteca compartilhada liblla.so
make OPENMP=1 ...   # paraleliza matrix_mul com OpenMP (qualquer alvo)
make clean      # remove os executáveis e as bibliotecas
```

O build padrão usa `-O2 -g -Wall -Wextra -Werror` — qualquer aviso interrompe a compilação. Com `make OPENMP=1` (por exemplo, `make OPENMP=1 bench-run`), o produto matricial é paralelizado com OpenMP sobre as colunas do resultado: não muda o comportamento do build padrão, só adiciona `-fopenmp` e mais desempenho em máquinas multi-core.

## Biblioteca estática e compartilhada

O alvo `make libs` compila os módulos da biblioteca (`vector`, `matrix`, `lu`,
`svd` e `lsystem`) em `obj/` com `-fPIC` e gera os dois formatos:

- `liblla.a`: biblioteca estática;
- `liblla.so`: biblioteca compartilhada.

Uso:

```sh
# estática (o código fica embutido no executável)
cc -Isrc meu_app.c -llla -lm

# compartilhada (o executável carrega liblla.so em tempo de execução)
cc -Isrc meu_app.c -llla -lm -Wl,-rpath,/home/user/lla
```

O build usa `-O2 -g -Wall -Wextra -Werror`, então qualquer aviso interrompe a compilação.

## Benchmark

`src/bech.c` mede o desempenho das operações de matriz em dimensões
quadradas crescentes. Por padrão roda com as dimensões **5, 10, 20 e 50**,
podendo ser executado com outras dimensões:

```sh
./bech            # 5, 10, 20, 50
./bech 3 7 100    # dimensões customizadas (1 a 512)
```

Para cada dimensão o benchmark gera duas matrizes pseudoaleatórias
diagonalmente dominantes (portanto não singulares) e mede, em regime de
tempo constante, as operações `matrix_mul`, `matrix_add`, `matrix_sub`,
`matrix_mul_scalar`, `matrix_transpose`, a decomposição LU, o determinante e
a inversa. A saída traz, por operação, o número de iterações, o tempo total,
os milissegundos por operação, o throughput e uma estimativa em MFLOP/s,
seguida de uma tabela de escala (ms por operação em cada dimensão) e dos
checksums usados para verificar a consistência dos resultados.

Como a biblioteca aloca o resultado de cada operação, os tempos medidos
incluem o custo de alocação e liberação de memória.

## Limites de dimensão

`src/lim.c` procura a maior matriz quadrada que a biblioteca ainda consegue
operar. Para cada operação, a busca dobra a dimensão (8, 16, 32, ...) e usa
o expoente de complexidade medido entre duas execuções para prever o custo
do próximo passo, de modo a não travar em operações O(n³).

```sh
./lim        # orçamento padrão de 10 s por execução
./lim 2.5    # orçamento de 2,5 s por execução
```

O programa relata três limites:

- **memória**: até onde as alocações da biblioteca têm sucesso (a busca é
  limitada a 60% da memória física da máquina);
- **tempo**: até onde cada operação cabe no orçamento de tempo;
- **numérico**: se o resultado ainda bate com uma implementação de
  referência independente no maior tamanho suportado.

Também é reportado o tamanho a partir do qual o determinante deixa de caber
em um `double` (cresce exponencialmente com n) e passa a devolver `inf`.

### Números medidos

Medidos com `./lim 10` (orçamento de 10 s por execução) em um Xeon E5-2695 v4
de 16 núcleos com 30 GB de memória, gcc 15.2.0 e `-O2`, com matrizes
quadradas aleatórias e diagonalmente dominantes. Os valores dependem do
hardware e do orçamento de tempo, mas a ordem de grandeza é a mesma em
qualquer máquina moderna:

| Operação | Maior `n` | Tempo em `n` | Limite |
| --- | --- | --- | --- |
| `matrix_add` | 16384 | 1,8 s | memória |
| `matrix_sub` | 16384 | 1,8 s | memória |
| `matrix_mul_scalar` | 16384 | 1,7 s | memória |
| `matrix_transpose` | 8192 | 2,5 s | tempo |
| `matrix_mul` | 2048 | 8,6 s | tempo |
| `lu_create_decomposition` | 2048 | 2,3 s | tempo |
| `matrix_determinant` | 2048 | 2,3 s | tempo |
| `matrix_inverse` | 2048 | 6,7 s | tempo |

Em todos esses tamanhos o resultado foi conferido contra uma implementação de
referência independente e confere (erro relativo da ordem de `1e-14` para LU
e inversa, zero para as operações exatas). O limite é de tempo e de memória,
não de correção.

**Antes e depois (v0.2.0 → v0.3.0, mesma máquina).** Os três gargalos da
versão anterior eram: uma alocação por coluna (uma matriz `n x n` era
`n` alocações), validação de índice a cada acesso dentro dos laços e o
producto em laço `i -> j -> k` com acesso não contíguo. Com o bloco denso,
acesso direto por ponteiros e o produto reordenado (`j -> k -> i`, saxpy por
colunas), em `n = 200`:

| Operação | v0.2.0 | v0.3.0 | Aceleração |
| --- | --- | --- | --- |
| `matrix_mul` | 190 ms | 5,8 ms | 33× |
| `matrix_add` | 1,95 ms | 0,10 ms | 19× |
| `matrix_sub` | 1,85 ms | 0,10 ms | 18× |
| `matrix_mul_scalar` | 1,31 ms | 0,093 ms | 14× |
| `matrix_transpose` | 0,69 ms | 0,11 ms | 6,6× |
| `lu_create_decomposition` | 99 ms | 1,7 ms | 57× |
| `matrix_determinant` | 95 ms | 1,7 ms | 57× |
| `matrix_inverse` | 415 ms | 5,4 ms | 77× |

Além do ganho por operação, o produto matricial passou de ~80 para ~2800
MFLOP/s em `n = 512`, e os limites do `lim` subiram de `n = 512` (produto,
LU, determinante) e `n = 256` (inversa) para `n = 2048` em todos eles.

**Limite de memória.** A biblioteca aloca uma matriz `n x n` como **um único
bloco contíguo** de `n*n` doubles — uma alocação por matriz:

| Alocação | Maior `n` | Tamanho |
| --- | --- | --- |
| 1 matriz | 32768 | 8,00 GB |
| 3 matrizes (pior caso de `matrix_mul`) | 16384 | 6,00 GB |
| 6 matrizes (folga da sonda; a inversa usa 4) | 16384 | 12,00 GB |

Quando a memória acaba, a biblioteca falha de forma graciosa: as funções
devolvem `NULL` com `errno = EMATRIXALLOCATION`.

**Limite numérico.** O determinante cresce exponencialmente com `n` e deixa
de caber em um `double` (`|det| > 1.8e308`) a partir de **n = 162** para
matrizes aleatórias desse tipo. Acima disso `matrix_determinant` devolve
`inf`, que não é um erro reportado pela API, e sim um resultado inválido.

### Regra prática

| `n` | Comportamento |
| --- | --- |
| até 100 | instantâneo, incluindo `matrix_inverse` |
| 100 a 500 | confortável para as operações cúbicas (milissegundos) |
| 500 a 2000 | `matrix_mul`, LU, determinante e inversa levam de segundos a dezenas de segundos |
| acima de 2000 | inviável para as operações cúbicas; as operações lineares e quadráticas continuam bem até onde a memória deixar |

Para uso recreativo, a faixa de 10 a 500 linhas cobre com folga o que se faz
com álgebra linear na prática, e é onde a biblioteca se comporta melhor.

## Tipos

### Vector

```c
typedef struct vector_rn_t {
  double* x;        // coordenadas
  unsigned int dim; // dimensão
} Vector;
```

### Matrix

A matriz é armazenada em um único bloco contíguo de doubles, no formato
column-major: `a(i, j)` fica em `data[j * nlins + i]`.

```c
typedef struct matrix_rn_t {
  double*      data;  // bloco contíguo column-major
  unsigned int nlins; // número de linhas
  unsigned int ncols; // número de colunas
} Matrix;
```

## API

### Vetores (`vector.h`)

| Função | Descrição |
| --- | --- |
| `vector_create_empty(dim)` | Cria um vetor de dimensão `dim` com todas as coordenadas zeradas |
| `vector_create(dim, coords)` | Cria um vetor a partir de um array de coordenadas |
| `vector_destroy(v)` | Libera a memória do vetor |
| `vector_add(a, b)` | Retorna `a + b` |
| `vector_sub(a, b)` | Retorna `a - b` |
| `vector_mul(k, a)` | Retorna `k * a` |
| `vector_inner_product(a, b)` | Retorna o produto interno `a · b` |
| `vector_norm_p(a, p)` | Retorna a norma `p` de `a`: `(Σ |xᵢ|^p)^(1/p)` |
| `vector_norm_inf(a)` | Retorna a norma infinita de `a`: `max(|xᵢ|)` |
| `vector_norm(a)` | Retorna a norma euclidiana de `a`: `sqrt(Σ xᵢ²)` |
| `vector_product(dim, vectors)` | Calcula o produto vetorial generalizado de `dim - 1` vetores de dimensão `dim` |
| `vector_print(v, label)` / `vector_println(v, label)` | Imprime o vetor no formato `(x1, x2, ...)`, prefixado com `label = ` quando o rótulo é informado |

`vector_inner_product` devolve `NAN` (com `errno` definido) quando os vetores
têm dimensões diferentes ou quando algum deles é nulo.

As funções de norma devolvem `NAN` (com `errno` definido) quando recebem um
vetor nulo ou, no caso de `vector_norm_p`, quando `p` é zero
(`errno = EINVALIDP`). `vector_norm` é equivalente a `vector_norm_p` com
`p = 2`.

#### Produto vetorial generalizado

`vector_product(dim, vectors)` recebe `dim - 1` vetores, todos com dimensão
`dim`, e devolve um novo vetor de dimensão `dim`. O vetor resultante é
ortogonal a todos os vetores de entrada:

```
resultado · entrada[j] = 0, para todo j
```

A implementação monta uma matriz `M` cujas linhas são os vetores de entrada e
calcula a primeira linha de cofatores:

```
resultado[i] = (-1)^i · det(M sem a coluna i),  i = 0, ..., dim - 1
```

Esse é o produto vetorial generalizado: em R³, com dois vetores de entrada, o
resultado coincide com o produto vetorial tradicional `A x B`. A orientação do
resultado segue a convenção dos determinantes; trocar a ordem de duas entradas
inverte o sinal do resultado.

O vetor devolvido é alocado pela biblioteca e deve ser liberado com
`vector_destroy`. Se `dim < 2`, a função devolve `NULL` e define
`errno = EINVALIDDIM`. Entradas inválidas ou falhas de alocação também podem
fazer a função devolver `NULL` com `errno` definido.

Exemplo em R³:

```c
double a_coords[] = { 1.0, 0.0, 0.0 };
double b_coords[] = { 0.0, 1.0, 0.0 };
Vector entradas[] = {
  { .dim = 3, .x = a_coords },
  { .dim = 3, .x = b_coords },
};

Vector* resultado = vector_product(3, entradas);
if (resultado == NULL) {
  return 1;
}

vector_println(&entradas[0], "A");
vector_println(&entradas[1], "B");
vector_println(resultado, "A x B");

vector_destroy(resultado);
```

### Matrizes (`matrix.h`)

| Função | Descrição |
| --- | --- |
| `matrix_create_empty(m, n)` | Cria uma matriz `m x n` zerada |
| `matrix_create(m, n, columns)` | Cria uma matriz `m x n` copiando os vetores coluna informados |
| `matrix_create_identity(m, n)` | Cria a matriz identidade `m x n` |
| `matrix_copy(a)` | Retorna uma cópia profunda da matriz |
| `matrix_destroy(m)` | Libera a memória da matriz |
| `matrix_get(m, i, j)` | Retorna o elemento `a(i,j)` |
| `matrix_set(m, i, j, v)` | Define o elemento `a(i,j)` e retorna a própria matriz |
| `matrix_validate_index(m, i, j)` | Verifica se os índices estão dentro das dimensões da matriz |
| `matrix_add(a, b)` | Retorna `A + B` |
| `matrix_sub(a, b)` | Retorna `A - B` |
| `matrix_mul(a, b)` | Retorna o produto `A x B` |
| `matrix_mul_scalar(k, a)` | Retorna `k * A` |
| `matrix_transpose(a)` | Retorna `Aᵀ` |
| `matrix_mul_vector(a, x)` | Retorna `A · x`, com `x` de dimensão igual ao número de colunas de `A` |
| `matrix_inverse(a)` | Retorna `A⁻¹` (matriz quadrada) via decomposição LU |
| `matrix_determinant(a)` | Retorna `det(A)` (matriz quadrada) via decomposição LU |
| `matrix_trace(a)` | Retorna o traço de `A` (soma da diagonal principal) |
| `matrix_rank(a)` | Retorna o posto de `A` via SVD (valores singulares acima de `eps · max(m,n) · σ₁`) |
| `matrix_print(m, label)` / `matrix_println(m, label)` | Imprime a matriz com um rótulo |

`matrix_determinant` devolve `0.0` para matrizes singulares e `NAN` (com `errno`
definido) para matriz nula, matriz não quadrada ou falha de alocação. Para
matrizes grandes o determinante pode estourar a faixa finita do `double` e o
retorno é `inf` — ver [Limites de dimensão](#limites-de-dimensão).

`matrix_add_column_vector_at(m, i, xlen, x)` **copia** os `xlen` valores de
`x` para a coluna `i` (desde a v0.3.0 não assume mais a posse do buffer —
o chamador continua responsável por liberá-lo, e `xlen` precisa ser igual ao
número de linhas). `matrix_get_vector_column_at` também mudou: devolve uma
**cópia** da coluna como um novo `Vector`, a liberar com `vector_destroy`.

Toda função que retorna um novo `Vector*` ou `Matrix*` aloca memória, que deve ser liberada pelo chamador com `vector_destroy` ou `matrix_destroy`.

### Decomposição LU (`lu.h`)

| Função | Descrição |
| --- | --- |
| `lu_create_decomposition(a)` | Retorna a decomposição LU de `a` (matriz quadrada) ou `NULL` com `errno` |
| `lu_destroy(lu)` | Libera a memória da decomposição |

O resultado é um `LU` com as três matrizes e o vetor de permutação:

```c
typedef struct LU_t {
  Matrix* l;   // triangular inferior, diagonal unitária
  Matrix* u;   // triangular superior
  Matrix* p;   // permutação
  int*    piv; // piv[i] = coluna do 1 na linha i de P (a permutação em vetor)
} LU;
```

A decomposição usa **pivoteamento parcial**, então a relação válida é

```
P * A = L * U      (equivalentemente, A = Pᵀ * L * U)
```

Como o pivoteamento reordena as linhas de `A`, `L * U` sozinho **não** devolve
`A`: é preciso aplicar a permutação. Exemplo para

```
A = | 2 1 1 |      L = | 1    0     0   |     U = | 8  7    3   |
    | 4 3 1 |          | 0.25 1     0   |         | 0 -0.75 0.25|
    | 8 7 3 |          | 0.50 0.67  1   |         | 0  0   -0.67|
```

tem-se `P = | 0 0 1 |`, `P * A = L * U` e `A = Pᵀ * L * U`.

Matrizes singulares (sem solução única) retornam `NULL` com `errno =
EMATRIXSINGULAR`.

### Valores singulares (`svd.h`)

| Função | Descrição |
| --- | --- |
| `svd_golub_reinsch(a)` | Retorna a decomposição em valores singulares de `a`: `A = U · diag(σ) · Vᵀ` |
| `svd_householder_bidiag(a)` | Reduz `a` à forma bidiagonal superior: `A = U · B · Vᵀ` |
| `svd_qr_bidiag_iter(hbdr)` | Iteração QR sobre a bidiagonal, acumulando as rotações em `U` e `V` |
| `svd_pos_proc(bdr)` | Normaliza os sinais, ordena `σ` e monta `VT` |
| `svd_destroy(r)` | Libera a memória de um `SVDResult` |
| `svd_householder_bidiag_destroy(r)` | Libera a memória de um `HouseholderBidiagResult` |
| `svd_qr_bidiag_destroy(r)` | Libera a memória de um `QRBidiagResult` |
| `svd_givens_rotation(f, g, c, s)` | Calcula a rotação de Givens `(c, s)` que leva `(f, g)` a `(r, 0)` |

O resultado da decomposição é:

```c
typedef struct svd_result_t {
  Matrix* U;      // m linhas x n colunas (m >= n)
  Vector* sigma;  // os n valores singulares, em ordem decrescente
  Matrix* VT;     // n x n, com Vᵀ
} SVDResult;
```

Com `A` de `m` linhas e `n` colunas (`m >= n`):

- `U` tem `m` linhas e `n` colunas, com colunas ortonormais;
- `sigma` tem os `n` valores singulares, não negativos e em ordem
  decrescente;
- `VT` é `n x n`, com linhas ortonormais;
- `A = U · diag(sigma) · VT` (reconstrução exata até o erro de arredondamento).

Matrizes com mais colunas que linhas são transpostas internamente e o
resultado é devolvido já no formato acima, com `min(m, n)` valores singulares.
A matriz recebida **não** é modificada.

#### Algoritmo

A implementação segue o método de Golub-Reinsch, em três etapas:

1. **Bidiagonalização de Householder.** Cada passo aplica um refletor
   `H = I - 2 v vᵀ` (à esquerda, nas linhas) ou `H = I - 2 v vᵀ` (à direita,
   nas colunas), reduzindo `A` a `B = Uᵀ · A · V`, bidiagonal superior.
2. **Iteração QR com deslocamento de Wilkinson.** A superdiagonal de `B` é
   anulada por rotações de Givens; o "calombo" (bulge) que aparece a cada
   passo é perseguido da diagonal principal até o canto da matriz, com
   deflações a cada subdiagonal que se anula.
3. **Pós-processamento.** Os sinais são normalizados para que `σ >= 0` e os
   valores singulares são ordenados em ordem decrescente.

A redução custa `O(m · n²)` e a iteração QR custa `O(n · (m + n))` por
varredura; na prática a decomposição converge em poucas varreduras.

#### O produto `v vᵀ` nunca é materializado

A etapa de Householder precisa de `v vᵀ`, mas esse produto **não** é
construído. Como `v` é um vetor (e não uma matriz), a transposta `vᵀ` é
apenas um produto interno: o refletor é aplicado coluna a coluna (ou linha a
linha) usando `vᵀ · x` calculado diretamente, sem alocar matrizes
intermediárias.

Isso evita a alocação de uma matriz `m x n` por refletor (a versão anterior
do código criava `matrix_transpose` e `matrix_mul` só para obter `v vᵀ`) e
troca o custo `O(m · n²)` da multiplicação completa por `O(m · n)` do
refletor aplicado diretamente. Em outras palavras: **não existe um "vetor
transposto" na biblioteca — a transposição de vetor é implícita no produto
interno**, que é a forma eficiente de calcular `vᵀ · x`.

#### Exemplo

```c
#include "matrix.h"
#include "svd.h"

int main(void) {
  double c1[] = { 1.0, 3.0 };
  double c2[] = { 2.0, 4.0 };
  Vector colunas[] = {
    { .dim = 2, .x = c1 },
    { .dim = 2, .x = c2 },
  };

  Matrix* a = matrix_create(2, 2, colunas);
  SVDResult* r = svd_golub_reinsch(a);
  if (r == NULL) {
    return 1;
  }

  matrix_println(r->U, "U");
  vector_println(r->sigma, "sigma");
  matrix_println(r->VT, "VT");

  svd_destroy(r);
  matrix_destroy(a);
  return 0;
}
```

#### Etapas intermediárias

Além de `svd_golub_reinsch`, as três etapas são públicas e podem ser usadas
separadamente. `svd_qr_bidiag_iter` assume a posse de `d`, `U` e `V` do
`HouseholderBidiagResult` recebido, e `svd_pos_proc` assume a posse de `d` e
`U` do `QRBidiagResult`; nos dois casos o chamador libera a estrutura
intermediária com a função `*_destroy` correspondente:

```c
HouseholderBidiagResult* hbdr = svd_householder_bidiag(a);
QRBidiagResult* bdr = svd_qr_bidiag_iter(hbdr);
svd_householder_bidiag_destroy(hbdr);

SVDResult* r = svd_pos_proc(bdr);
svd_qr_bidiag_destroy(bdr);
```

### Sistemas lineares (`lsystem.h`)

| Função | Descrição |
| --- | --- |
| `lsystem_solve(a, b)` | Resolve o sistema linear `A · x = b` |

A escolha do método é automática:

- matriz **quadrada e regular**: decomposição LU com pivoteamento parcial
  (substituições direta e inversa);
- matriz **singular** ou **retangular** (excesso ou pouco definida):
  mínimos quadrados via SVD com corte de valores singulares abaixo de
  `eps · max(m, n) · σ₁` — a solução devolvida é a de **norma mínima**.

```c
double c1[] = { 3.0, 2.0 };
double c2[] = { -1.0, 4.0 };
Vector colunas[] = {
  { .dim = 2, .x = c1 },
  { .dim = 2, .x = c2 },
};

Matrix* a = matrix_create(2, 2, colunas);
double bx[] = { 4.0, 26.0 };
Vector* b = vector_create(2, bx);

Vector* x = lsystem_solve(a, b); // x = (3, 5)
```

## Exemplo

```c
#include "matrix.h"

int main(void) {
  double c1[] = { 1.0, 3.0 };
  double c2[] = { 2.0, 4.0 };
  Vector colunas[] = {
    { .dim = 2, .x = c1 },
    { .dim = 2, .x = c2 },
  };

  // A = | 1  2 |
  //     | 3  4 |
  Matrix* a = matrix_create(2, 2, colunas);
  Matrix* b = matrix_mul_scalar(2.0, a);
  Matrix* c = matrix_mul(a, b);

  matrix_println(a, "A");
  matrix_println(b, "B = 2A");
  matrix_println(c, "C = A x B");

  matrix_destroy(c);
  matrix_destroy(b);
  matrix_destroy(a);
  return 0;
}
```

## Tratamento de erros

Em caso de falha, as funções retornam `NULL` (ou `NAN`, no caso de `matrix_get`,
`vector_inner_product` e das funções de norma) e definem `errno` com um dos
códigos abaixo:

| Código | Valor | Significado |
| --- | --- | --- |
| `EINVALIDDIM` | 1000 | Dimensão de vetor inválida ou incompatível |
| `EVECTORALLOCATION` | 1001 | Falha ao alocar vetor |
| `EVECTORNULL` | 1002 | Vetor nulo passado como argumento |
| `EINVALIDP` | 1003 | Ordem `p` inválida para a norma de vetor |
| `EMATRIXALLOCATION` | 2000 | Falha ao alocar matriz |
| `EMATRIXINVALIDADDDIM` | 2001 | Dimensões incompatíveis para soma/subtração |
| `EMATRIXINVALIDLINEINDEX` | 2002 | Índice de linha fora dos limites |
| `EMATRIXINVALIDCOLINDEX` | 2003 | Índice de coluna fora dos limites |
| `EMATRIXINVALIDMUL` | 2004 | Dimensões incompatíveis para multiplicação |
| `EMATRIXNULL` | 2999 | Matriz nula passada como argumento |
| `EMATRIXINVALIDDIM` | 2005 | Matriz não quadrada ou de dimensão inválida |
| `EMATRIXSINGULAR` | 2006 | Matriz singular, sem decomposição LU única |
| `EMATRIXINVALIDDETERMINANT` | 2007 | Falha ao calcular o determinante |
| `ESVDNULL` | 4000 | Argumento nulo na decomposição em valores singulares |
| `ESVDINVALIDDIM` | 4001 | Matriz sem linhas ou sem colunas |
| `ESVDALLOCATION` | 4002 | Falha ao alocar memória na decomposição |
| `ESVDNOTCONVERGED` | 4003 | A iteração QR não convergiu dentro do limite de varreduras |
| `ELSYSTEMNULL` | 3000 | Argumento nulo na resolução do sistema linear |
| `ELSYSTEMINVALIDDIM` | 3001 | Dimensões incompatíveis entre `A` e `b`, ou nulas |

## Próximos passos

- Separar os testes do executável principal
- Aproveitar a faixa ampliada (n até ~2000) como base das próximas
  bibliotecas: estatística (`lstat`), probabilidade (`lprob`) e
  machine learning (`lml`)

