# LLA - Little Linear Algebra

Biblioteca pequena em C para operações de álgebra linear com vetores e matrizes de números reais (`double`).

O objetivo é manter o código simples, sem dependências externas além da biblioteca padrão do C, e fácil de incorporar em outros projetos.

Por ser um projeto recreativo, a biblioteca prioriza clareza sobre desempenho bruto. Na prática isso significa que as matrizes são guardadas como vetores coluna, com uma alocação por coluna, e que todo acesso valida os índices — o que torna as operações O(n³) (produto, LU, determinante, inversa) bem mais lentas do que uma implementação otimizada. A seção [Limites de dimensão](#limites-de-dimensão) traz os números concretos: em resumo, matrizes de até algumas centenas de linhas são tranquilamente usáveis, e a partir de ~1000 as operações cúbicas começam a ficar pesadas.

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
make clean      # remove os executáveis e as bibliotecas
```

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

O build usa `-Wall -Wextra -Werror`, então qualquer aviso interrompe a compilação.

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

Medidos com `./lim 10` (orçamento de 10 s por execução) em um Mac com 8 GB de
memória, com matrizes quadradas aleatórias e diagonalmente dominantes. Os
valores dependem do hardware e do orçamento de tempo, mas a ordem de grandeza
é a mesma em qualquer máquina moderna:

| Operação | Maior `n` | Tempo em `n` | Limite |
| --- | --- | --- | --- |
| `matrix_transpose` | 8192 | 1,4 s | memória |
| `matrix_add` | 8192 | 6,4 s | memória |
| `matrix_sub` | 8192 | 8,0 s | memória |
| `matrix_mul_scalar` | 8192 | 3,4 s | memória |
| `matrix_mul` | 512 | 5,1 s | tempo |
| `lu_create_decomposition` | 512 | 2,8 s | tempo |
| `matrix_determinant` | 512 | 2,8 s | tempo |
| `matrix_inverse` | 256 | 1,5 s | tempo |

Em todos esses tamanhos o resultado foi conferido contra uma implementação de
referência independente e confere (erro relativo da ordem de `1e-15` para LU
e inversa, zero para as operações exatas). O limite é de tempo e de memória,
não de correção.

**Limite de memória.** A biblioteca aloca uma matriz `n x n` como `n` vetores
coluna, ou seja, `n` alocações separadas:

| Alocação | Maior `n` | Tamanho |
| --- | --- | --- |
| 1 matriz | 16384 | 2,00 GB |
| 3 matrizes (pior caso de `matrix_mul`) | 8192 | 1,50 GB |
| 6 matrizes (pior caso de `matrix_inverse`) | 8192 | 3,00 GB |

Quando a memória acaba, a biblioteca falha de forma graciosa: as funções
devolvem `NULL` com `errno = EMATRIXALLOCATION`.

**Limite numérico.** O determinante cresce exponencialmente com `n` e deixa
de caber em um `double` (`|det| > 1.8e308`) a partir de **n = 162** para
matrizes aleatórias desse tipo. Acima disso `matrix_determinant` devolve
`inf`, que não é um erro reportado pela API, e sim um resultado inválido.

### Regra prática

| `n` | Comportamento |
| --- | --- |
| até 50 | instantâneo, incluindo `matrix_inverse` |
| 50 a 500 | confortável para as operações cúbicas (segundos) |
| 500 a 1000 | `matrix_mul`, LU, determinante e inversa levam de segundos a dezenas de segundos |
| acima de 1000 | inviável para as operações cúbicas; as operações lineares e quadráticas continuam bem |

Para uso recreativo, a faixa de 10 a 200 linhas cobre com folga o que se faz
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

A matriz é armazenada como um conjunto de vetores coluna.

```c
typedef struct matrix_rn_t {
  Vector* columns;    // vetores coluna
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

Toda função que retorna um novo `Vector*` ou `Matrix*` aloca memória, que deve ser liberada pelo chamador com `vector_destroy` ou `matrix_destroy`.

### Decomposição LU (`lu.h`)

| Função | Descrição |
| --- | --- |
| `lu_create_decomposition(a)` | Retorna a decomposição LU de `a` (matriz quadrada) ou `NULL` com `errno` |
| `lu_destroy(lu)` | Libera a memória da decomposição |

O resultado é um `LU` com três matrizes:

```c
typedef struct LU_t {
  Matrix* l; // triangular inferior, diagonal unitária
  Matrix* u; // triangular superior
  Matrix* p; // permutação
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
- Reduzir o overhead das operações internas (validação de índice a cada
  acesso e alocação por coluna), o que ampliaria os limites de dimensão
  descritos em [Limites de dimensão](#limites-de-dimensão)

