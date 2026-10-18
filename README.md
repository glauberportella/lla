# LLA - Little Linear Algebra

Biblioteca pequena em C para operações de álgebra linear com vetores e matrizes de números reais (`double`).

O objetivo é manter o código simples, sem dependências externas além da biblioteca padrão do C, e fácil de incorporar em outros projetos.

## Estrutura

```
.
├── Makefile
└── src
    ├── vector.h / vector.c   # tipo Vector e operações sobre vetores
    ├── matrix.h / matrix.c   # tipo Matrix e operações sobre matrizes
    └── main.c                # testes
```

## Compilação

Requer `gcc` (ou outro compilador C99 compatível) e `make`.

```sh
make build   # gera o executável ./al
make run     # executa os testes
```

O build usa `-Wall -Wextra -Werror`, então qualquer aviso interrompe a compilação.

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
| `vector_print(v)` / `vector_println(v)` | Imprime o vetor no formato `(x1, x2, ...)` |

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
| `matrix_print(m, label)` / `matrix_println(m, label)` | Imprime a matriz com um rótulo |

Toda função que retorna um novo `Vector*` ou `Matrix*` aloca memória, que deve ser liberada pelo chamador com `vector_destroy` ou `matrix_destroy`.

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

Em caso de falha, as funções retornam `NULL` (ou `NAN`, no caso de `matrix_get`) e definem `errno` com um dos códigos abaixo:

| Código | Valor | Significado |
| --- | --- | --- |
| `EINVALIDDIM` | 1000 | Dimensão de vetor inválida ou incompatível |
| `EVECTORALLOCATION` | 1001 | Falha ao alocar vetor |
| `EVECTORNULL` | 1002 | Vetor nulo passado como argumento |
| `EMATRIXALLOCATION` | 2000 | Falha ao alocar matriz |
| `EMATRIXINVALIDADDDIM` | 2001 | Dimensões incompatíveis para soma/subtração |
| `EMATRIXINVALIDLINEINDEX` | 2002 | Índice de linha fora dos limites |
| `EMATRIXINVALIDCOLINDEX` | 2003 | Índice de coluna fora dos limites |
| `EMATRIXINVALIDMUL` | 2004 | Dimensões incompatíveis para multiplicação |
| `EMATRIXNULL` | 2999 | Matriz nula passada como argumento |

## Próximos passos

- Produto interno, norma e produto vetorial
- Transposta
- Determinante e matriz inversa
- Resolução de sistemas lineares (eliminação de Gauss, decomposição LU)
- Build como biblioteca estática/compartilhada
- Separar os testes do executável principal
