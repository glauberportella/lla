# Exemplos

Sete exemplos independentes, cada um mostrando um uso real da lib com um
contexto de problema (fisica, jogos, nutricao, previsao e processamento de
sinais). Todos compilam direto contra os fontes em `../src`: nao e preciso
gerar `liblla` antes.

## Como executar

```sh
cd examples
make          # compila os sete exemplos
make run      # compila e executa todos
make clean
```

Ou, para compilar e executar um exemplo especifico:

```sh
make 04_decomposicao_lu && ./04_decomposicao_lu
```

## O que cada exemplo mostra

| exemplo | tema | features da lib |
|---------|------|-----------------|
| `01_vetores.c` | forcas, trabalho e torque | `vector_create(_empty)`, `vector_copy`, `vector_set`, `vector_add`, `vector_sub`, `vector_mul`, `vector_inner_product`, `vector_norm*`, `vector_product`, `vector_print(ln)` |
| `02_matrizes.c` | transformacoes 2D de um sprite | `matrix_create(_empty/_identity)`, `matrix_get/set`, `matrix_validate_index`, `matrix_add_column_vector_at`, `matrix_remove_column_vector_at`, `matrix_get_line_vectors`, `matrix_add`, `matrix_sub`, `matrix_mul`, `matrix_mul_scalar`, `matrix_mul_vector`, `matrix_transpose`, `matrix_trace`, `matrix_print(ln)` |
| `03_determinante_inversa.c` | estabilidade de estrutura | `matrix_determinant`, `matrix_inverse`, diagnostico de singularidade e tratamento de erro (`errno` + retorno `NULL`) |
| `04_decomposicao_lu.c` | circuito eletrico | `lu_create_decomposition`, verificacao `P*A = L*U`, uso manual de `L` e `U` com substituicao forward/back e cross-check com `lsystem_solve` |
| `05_sistema_linear.c` | planejamento de dieta | `lsystem_solve` em sistema quadrado 3x3 (caminho LU) |
| `06_minimos_quadrados.c` | previsao de tempo de tarefa | `lsystem_solve` em sistema retangular 6x2 (caminho SVD/pseudoinversa) |
| `07_svd.c` | leituras de um sensor | `svd_golub_reinsch` e as 3 etapas isoladas (`svd_householder_bidiag`, `svd_qr_bidiag_iter`, `svd_pos_proc`), destrutores de cada resultado, `svd_givens_rotation`, `lla_signal` |

## Boas praticas

- **Memoria**: `matrix_create` copia os dados recebidos, mas o resultado de
  `matrix_get_line_vectors` tem de ser liberado pelo chamador (os `x` de cada
  linha e o array em si). Veja `02_matrizes.c`.
- **Posse**: `matrix_add_column_vector_at` assume a posse do buffer passado,
  por isso o buffer e alocado com `malloc`. As funcoes de SVD em etapas assumem
  a-posse do resultado intermediario, que deve sempre ir para o destrutor
  correspondente. Veja `07_svd.c`.
- **Erros**: toda funcao que devolve ponteiro retorna `NULL` com `errno`
  definido; os codigos de erro estao nos cabecalhos (`src/*.h`).
