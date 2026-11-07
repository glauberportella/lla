# Examples

Seven independent examples, each showing a real use of the library in a
problem context (physics, games, nutrition, forecasting and signal
processing). All of them compile directly against the sources in `../src`:
there is no need to build `liblla` first.

## How to run

```sh
cd examples
make          # builds the seven examples
make run      # builds and runs all of them
make clean
```

Or, to build and run a specific example:

```sh
make 04_decomposicao_lu && ./04_decomposicao_lu
```

## What each example shows

| example | topic | library features |
|---------|-------|------------------|
| `01_vetores.c` | forces, work and torque | `vector_create(_empty)`, `vector_copy`, `vector_set`, `vector_add`, `vector_sub`, `vector_mul`, `vector_inner_product`, `vector_norm*`, `vector_product`, `vector_print(ln)` |
| `02_matrizes.c` | 2D sprite transformations | `matrix_create(_empty/_identity)`, `matrix_get/set`, `matrix_validate_index`, `matrix_add_column_vector_at`, `matrix_remove_column_vector_at`, `matrix_get_line_vectors`, `matrix_add`, `matrix_sub`, `matrix_mul`, `matrix_mul_scalar`, `matrix_mul_vector`, `matrix_transpose`, `matrix_trace`, `matrix_print(ln)` |
| `03_determinante_inversa.c` | structural stability | `matrix_determinant`, `matrix_inverse`, singularity diagnostics and error handling (`errno` + `NULL` return) |
| `04_decomposicao_lu.c` | electrical circuit | `lu_create_decomposition`, `P*A = L*U` verification, manual use of `L` and `U` with forward/back substitution and cross-check with `lsystem_solve` |
| `05_sistema_linear.c` | diet planning | `lsystem_solve` on a square 3x3 system (LU path) |
| `06_minimos_quadrados.c` | task time forecasting | `lsystem_solve` on a rectangular 6x2 system (SVD/pseudoinverse path) |
| `07_svd.c` | sensor readings | `svd_golub_reinsch` and the 3 isolated stages (`svd_householder_bidiag`, `svd_qr_bidiag_iter`, `svd_pos_proc`), destructors for each result, `svd_givens_rotation`, `lla_signal` |

## Best practices

- **Memory**: `matrix_create` copies the received data, but the result of
  `matrix_get_line_vectors` must be freed by the caller (the `x` of each row
  and the array itself). See `02_matrizes.c`.
- **Ownership**: `matrix_add_column_vector_at` copies the values from the
  passed buffer, which remains the caller's responsibility (the example frees
  it with `free`). The staged SVD functions take ownership of the intermediate
  result, which must always be passed to the corresponding destructor. See
  `07_svd.c`.
- **Errors**: every function that returns a pointer returns `NULL` with
  `errno` set; the error codes are in the headers (`src/*.h`).
