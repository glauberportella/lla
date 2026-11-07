# LLA - Little Linear Algebra

Version history: [CHANGELOG.md](CHANGELOG.md).

A small C library for linear algebra with vectors and matrices of real
numbers (`double`).

The goal is to keep the code simple, with no dependencies beyond the C
standard library, and easy to embed in other projects.

Because this is a recreational project, the library favors clarity over raw
performance. Since version 0.3.0 matrices are stored in a **single contiguous
block of doubles** (column-major layout), one allocation per matrix, and the
operations access that block directly through pointers — index validation
happens only in the public API (`matrix_get`/`matrix_set`). The
[Dimension limits](#dimension-limits) section has the concrete numbers: in
short, matrices with up to a few hundred rows are instantaneous, and the cubic
operations (product, LU, determinant, inverse) remain practical up to
n ≈ 2000.

## Layout

```
.
├── Makefile
├── README.md
├── CHANGELOG.md
├── examples                 # usage examples (see examples/README.md)
└── src
    ├── vector.h / vector.c   # Vector type and vector operations
    ├── matrix.h / matrix.c   # Matrix type and matrix operations
    ├── lu.h / lu.c           # LU decomposition with partial pivoting
    ├── svd.h / svd.c         # singular value decomposition
    ├── lsystem.h / lsystem.c # linear system solving (LU and SVD least squares)
    ├── main.c                # tests
    ├── bech.c                # performance benchmark
    └── lim.c                 # dimension limits
```

## Building

Requires `gcc` (or another C99-compatible compiler) and `make`.

```sh
make build      # builds the ./al executable
make run        # runs the tests
make bench      # builds the ./bech executable
make bench-run  # builds and runs the benchmark
make lim        # builds the ./lim executable
make lim-run    # builds and runs the dimension limit search
make libs       # builds liblla.a and liblla.so
make static     # builds only the static library liblla.a
make shared     # builds only the shared library liblla.so
make OPENMP=1 ...   # parallelizes matrix_mul with OpenMP (any target)
make clean      # removes executables and libraries
```

The default build uses `-O2 -g -Wall -Wextra -Werror` — any warning stops the
compilation. With `make OPENMP=1` (for example, `make OPENMP=1 bench-run`),
the matrix product is parallelized with OpenMP over the columns of the result:
the default build behavior does not change, it only adds `-fopenmp` and more
performance on multi-core machines.

## Static and shared library

The `make libs` target compiles the library modules (`vector`, `matrix`, `lu`,
`svd` and `lsystem`) into `obj/` with `-fPIC` and produces both formats:

- `liblla.a`: static library;
- `liblla.so`: shared library.

Usage:

```sh
# static (the code is embedded into the executable)
cc -Isrc my_app.c -llla -lm

# shared (the executable loads liblla.so at runtime)
cc -Isrc my_app.c -llla -lm -Wl,-rpath,/home/user/lla
```

The build uses `-O2 -g -Wall -Wextra -Werror`, so any warning stops the
compilation.

## Benchmark

`src/bech.c` measures the performance of the matrix operations on increasing
square dimensions. By default it runs with dimensions **5, 10, 20 and 50**,
and it can be run with other dimensions:

```sh
./bech            # 5, 10, 20, 50
./bech 3 7 100    # custom dimensions (1 to 512)
```

For each dimension the benchmark generates two pseudo-random, diagonally
dominant matrices (therefore non-singular) and measures, in a steady-state
regime, the operations `matrix_mul`, `matrix_add`, `matrix_sub`,
`matrix_mul_scalar`, `matrix_transpose`, the LU decomposition, the determinant
and the inverse. The output shows, for each operation, the number of
iterations, the total time, the milliseconds per operation, the throughput and
an estimate in MFLOP/s, followed by a scaling table (ms per operation at each
dimension) and the checksums used to verify the consistency of the results.

Since the library allocates the result of every operation, the measured times
include the cost of allocating and freeing memory.

## Dimension limits

`src/lim.c` searches for the largest square matrix the library can still
operate on. For each operation, the search doubles the dimension (8, 16, 32,
...) and uses the complexity exponent measured between two runs to predict
the cost of the next step, so that it does not hang on O(n³) operations.

```sh
./lim        # default budget of 10 s per run
./lim 2.5    # budget of 2.5 s per run
```

The program reports three limits:

- **memory**: how far the library allocations succeed (the search is limited
  to 60% of the machine's physical memory);
- **time**: how far each operation fits within the time budget;
- **numerical**: whether the result still matches an independent reference
  implementation at the largest supported size.

It also reports the size from which the determinant no longer fits in a
`double` (it grows exponentially with n) and starts returning `inf`.

### Measured numbers

Measured with `./lim 10` (10 s budget per run) on a Xeon E5-2695 v4 with 16
cores and 30 GB of memory, gcc 15.2.0 and `-O2`, with random, diagonally
dominant square matrices. The values depend on the hardware and on the time
budget, but the order of magnitude is the same on any modern machine:

| Operation | Largest `n` | Time at `n` | Limit |
| --- | --- | --- | --- |
| `matrix_add` | 16384 | 1.8 s | memory |
| `matrix_sub` | 16384 | 1.8 s | memory |
| `matrix_mul_scalar` | 16384 | 1.7 s | memory |
| `matrix_transpose` | 8192 | 2.5 s | time |
| `matrix_mul` | 2048 | 8.6 s | time |
| `lu_create_decomposition` | 2048 | 2.3 s | time |
| `matrix_determinant` | 2048 | 2.3 s | time |
| `matrix_inverse` | 2048 | 6.7 s | time |

At all of those sizes the result was checked against an independent reference
implementation and matches (relative error on the order of `1e-14` for LU and
the inverse, zero for the exact operations). The limit is time and memory, not
correctness.

**Before and after (v0.2.0 → v0.3.0, same machine).** The three bottlenecks of
the previous version were: one allocation per column (an `n x n` matrix was
`n` allocations), index validation on every access inside the loops, and the
`i -> j -> k` loop product with non-contiguous access. With the dense block,
direct pointer access and the reordered product (`j -> k -> i`, saxpy over the
columns), at `n = 200`:

| Operation | v0.2.0 | v0.3.0 | Speedup |
| --- | --- | --- | --- |
| `matrix_mul` | 190 ms | 5.8 ms | 33× |
| `matrix_add` | 1.95 ms | 0.10 ms | 19× |
| `matrix_sub` | 1.85 ms | 0.10 ms | 18× |
| `matrix_mul_scalar` | 1.31 ms | 0.093 ms | 14× |
| `matrix_transpose` | 0.69 ms | 0.11 ms | 6.6× |
| `lu_create_decomposition` | 99 ms | 1.7 ms | 57× |
| `matrix_determinant` | 95 ms | 1.7 ms | 57× |
| `matrix_inverse` | 415 ms | 5.4 ms | 77× |

Besides the per-operation gain, the matrix product went from ~80 to ~2800
MFLOP/s at `n = 512`, and the `lim` limits rose from `n = 512` (product, LU,
determinant) and `n = 256` (inverse) to `n = 2048` for all of them.

**Memory limit.** The library allocates an `n x n` matrix as a **single
contiguous block** of `n*n` doubles — one allocation per matrix:

| Allocation | Largest `n` | Size |
| --- | --- | --- |
| 1 matrix | 32768 | 8.00 GB |
| 3 matrices (worst case of `matrix_mul`) | 16384 | 6.00 GB |
| 6 matrices (probe headroom; the inverse uses 4) | 16384 | 12.00 GB |

When memory runs out, the library fails gracefully: the functions return
`NULL` with `errno = EMATRIXALLOCATION`.

**Numerical limit.** The determinant grows exponentially with `n` and no
longer fits in a `double` (`|det| > 1.8e308`) starting at **n = 162** for
random matrices of this kind. Above that, `matrix_determinant` returns `inf`,
which is not an error reported by the API, but an invalid result.

### Rule of thumb

| `n` | Behavior |
| --- | --- |
| up to 100 | instantaneous, including `matrix_inverse` |
| 100 to 500 | comfortable for the cubic operations (milliseconds) |
| 500 to 2000 | `matrix_mul`, LU, determinant and inverse take from seconds to tens of seconds |
| above 2000 | impractical for the cubic operations; linear and quadratic operations remain fine until memory runs out |

For recreational use, the 10 to 500 row range comfortably covers what is done
with linear algebra in practice, and it is where the library performs best.

## Types

Indices, dimensions and counts use `size_t`, the standard C library type for
sizes (`sizeof`, `strlen`, `malloc`) and array indices. The only parameter
that remains `unsigned int` is the order `p` of the norms
(`vector_norm_p`), which is not a size.

### Vector

```c
typedef struct vector_rn_t {
  double* x;   // coordinates
  size_t  dim; // dimension
} Vector;
```

### Matrix

The matrix is stored in a single contiguous block of doubles, in column-major
layout: `a(i, j)` lives at `data[j * nlins + i]`.

```c
typedef struct matrix_rn_t {
  double* data;  // contiguous column-major block
  size_t  nlins; // number of rows
  size_t  ncols; // number of columns
} Matrix;
```

## API

### Vectors (`vector.h`)

| Function | Description |
| --- | --- |
| `vector_create_empty(dim)` | Creates a vector of dimension `dim` with all coordinates set to zero |
| `vector_create(dim, coords)` | Creates a vector from an array of coordinates |
| `vector_destroy(v)` | Frees the vector memory |
| `vector_add(a, b)` | Returns `a + b` |
| `vector_sub(a, b)` | Returns `a - b` |
| `vector_mul(k, a)` | Returns `k * a` |
| `vector_inner_product(a, b)` | Returns the inner product `a · b` |
| `vector_norm_p(a, p)` | Returns the `p`-norm of `a`: `(Σ |xᵢ|^p)^(1/p)` |
| `vector_norm_inf(a)` | Returns the infinity norm of `a`: `max(|xᵢ|)` |
| `vector_norm(a)` | Returns the Euclidean norm of `a`: `sqrt(Σ xᵢ²)` |
| `vector_product(dim, vectors)` | Computes the generalized cross product of `dim - 1` vectors of dimension `dim` |
| `vector_print(v, label)` / `vector_println(v, label)` | Prints the vector as `(x1, x2, ...)`, prefixed with `label = ` when a label is given |

`vector_inner_product` returns `NAN` (with `errno` set) when the vectors have
different dimensions or when either one is null.

The norm functions return `NAN` (with `errno` set) when they receive a null
vector or, in the case of `vector_norm_p`, when `p` is zero
(`errno = EINVALIDP`). `vector_norm` is equivalent to `vector_norm_p` with
`p = 2`.

#### Generalized cross product

`vector_product(dim, vectors)` receives `dim - 1` vectors, all with dimension
`dim`, and returns a new vector of dimension `dim`. The resulting vector is
orthogonal to all input vectors:

```
result · input[j] = 0, for every j
```

The implementation builds a matrix `M` whose rows are the input vectors and
computes the first row of cofactors:

```
result[i] = (-1)^i · det(M without column i),  i = 0, ..., dim - 1
```

This is the generalized cross product: in R³, with two input vectors, the
result matches the traditional cross product `A x B`. The orientation of the
result follows the determinant convention; swapping the order of two inputs
flips the sign of the result.

The returned vector is allocated by the library and must be freed with
`vector_destroy`. If `dim < 2`, the function returns `NULL` and sets
`errno = EINVALIDDIM`. Invalid inputs or allocation failures can also make
the function return `NULL` with `errno` set.

Example in R³:

```c
double a_coords[] = { 1.0, 0.0, 0.0 };
double b_coords[] = { 0.0, 1.0, 0.0 };
Vector inputs[] = {
  { .dim = 3, .x = a_coords },
  { .dim = 3, .x = b_coords },
};

Vector* result = vector_product(3, inputs);
if (result == NULL) {
  return 1;
}

vector_println(&inputs[0], "A");
vector_println(&inputs[1], "B");
vector_println(result, "A x B");

vector_destroy(result);
```

### Matrices (`matrix.h`)

| Function | Description |
| --- | --- |
| `matrix_create_empty(m, n)` | Creates a zeroed `m x n` matrix |
| `matrix_create(m, n, columns)` | Creates an `m x n` matrix by copying the given column vectors |
| `matrix_create_identity(m, n)` | Creates the `m x n` identity matrix |
| `matrix_copy(a)` | Returns a deep copy of the matrix |
| `matrix_destroy(m)` | Frees the matrix memory |
| `matrix_get(m, i, j)` | Returns the element `a(i,j)` |
| `matrix_set(m, i, j, v)` | Sets the element `a(i,j)` and returns the matrix itself |
| `matrix_validate_index(m, i, j)` | Checks whether the indices are within the matrix dimensions |
| `matrix_add(a, b)` | Returns `A + B` |
| `matrix_sub(a, b)` | Returns `A - B` |
| `matrix_mul(a, b)` | Returns the product `A x B` |
| `matrix_mul_scalar(k, a)` | Returns `k * A` |
| `matrix_transpose(a)` | Returns `Aᵀ` |
| `matrix_mul_vector(a, x)` | Returns `A · x`, with `x` of dimension equal to the number of columns of `A` |
| `matrix_inverse(a)` | Returns `A⁻¹` (square matrix) via LU decomposition |
| `matrix_determinant(a)` | Returns `det(A)` (square matrix) via LU decomposition |
| `matrix_trace(a)` | Returns the trace of `A` (sum of the main diagonal) |
| `matrix_rank(a)` | Returns the rank of `A` via SVD (singular values above `eps · max(m,n) · σ₁`); on error it returns `SIZE_MAX` with `errno` set |
| `matrix_print(m, label)` / `matrix_println(m, label)` | Prints the matrix with a label |

`matrix_determinant` returns `0.0` for singular matrices and `NAN` (with
`errno` set) for a null matrix, a non-square matrix or an allocation failure.
For large matrices the determinant can overflow the finite range of a `double`
and the return value is `inf` — see
[Dimension limits](#dimension-limits).

`matrix_add_column_vector_at(m, i, xlen, x)` **copies** the `xlen` values of
`x` into column `i` (since v0.3.0 it no longer takes ownership of the buffer —
the caller remains responsible for freeing it, and `xlen` must equal the
number of rows). `matrix_get_vector_column_at` also changed: it returns a
**copy** of the column as a new `Vector`, to be freed with `vector_destroy`.

Every function that returns a new `Vector*` or `Matrix*` allocates memory,
which must be freed by the caller with `vector_destroy` or `matrix_destroy`.

### LU decomposition (`lu.h`)

| Function | Description |
| --- | --- |
| `lu_create_decomposition(a)` | Returns the LU decomposition of `a` (square matrix) or `NULL` with `errno` |
| `lu_destroy(lu)` | Frees the decomposition memory |

The result is an `LU` with the three matrices and the permutation vector:

```c
typedef struct LU_t {
  Matrix* l;   // lower triangular, unit diagonal
  Matrix* u;   // upper triangular
  Matrix* p;   // permutation
  size_t* piv; // piv[i] = column of the 1 in row i of P (the permutation as a vector)
} LU;
```

The decomposition uses **partial pivoting**, so the valid relation is

```
P * A = L * U      (equivalently, A = Pᵀ * L * U)
```

Because pivoting reorders the rows of `A`, `L * U` alone does **not** return
`A`: the permutation must be applied. For example, for

```
A = | 2 1 1 |      L = | 1    0     0   |     U = | 8  7    3   |
    | 4 3 1 |          | 0.25 1     0   |         | 0 -0.75 0.25|
    | 8 7 3 |          | 0.50 0.67  1   |         | 0  0   -0.67|
```

we have `P = | 0 0 1 |`, `P * A = L * U` and `A = Pᵀ * L * U`.

Singular matrices (no unique solution) return `NULL` with
`errno = EMATRIXSINGULAR`.

### Singular values (`svd.h`)

| Function | Description |
| --- | --- |
| `svd_golub_reinsch(a)` | Returns the singular value decomposition of `a`: `A = U · diag(σ) · Vᵀ` |
| `svd_householder_bidiag(a)` | Reduces `a` to upper bidiagonal form: `A = U · B · Vᵀ` |
| `svd_qr_bidiag_iter(hbdr)` | QR iteration on the bidiagonal matrix, accumulating the rotations in `U` and `V` |
| `svd_pos_proc(bdr)` | Normalizes the signs, sorts `σ` and builds `VT` |
| `svd_destroy(r)` | Frees an `SVDResult` |
| `svd_householder_bidiag_destroy(r)` | Frees a `HouseholderBidiagResult` |
| `svd_qr_bidiag_destroy(r)` | Frees a `QRBidiagResult` |
| `svd_givens_rotation(f, g, c, s)` | Computes the Givens rotation `(c, s)` that takes `(f, g)` to `(r, 0)` |

The decomposition result is:

```c
typedef struct svd_result_t {
  Matrix* U;      // m rows x n columns (m >= n)
  Vector* sigma;  // the n singular values, in descending order
  Matrix* VT;     // n x n, with Vᵀ
} SVDResult;
```

With `A` having `m` rows and `n` columns (`m >= n`):

- `U` has `m` rows and `n` columns, with orthonormal columns;
- `sigma` has the `n` singular values, non-negative and in descending order;
- `VT` is `n x n`, with orthonormal rows;
- `A = U · diag(sigma) · VT` (exact reconstruction up to rounding error).

Matrices with more columns than rows are transposed internally and the result
is returned in the format above, with `min(m, n)` singular values. The input
matrix is **not** modified.

#### Algorithm

The implementation follows the Golub-Reinsch method, in three stages:

1. **Householder bidiagonalization.** Each step applies a reflector
   `H = I - 2 v vᵀ` (on the left, over the rows) or `H = I - 2 v vᵀ` (on the
   right, over the columns), reducing `A` to `B = Uᵀ · A · V`, upper
   bidiagonal.
2. **QR iteration with Wilkinson shift.** The superdiagonal of `B` is zeroed
   by Givens rotations; the bulge that appears at each step is chased from the
   main diagonal to the corner of the matrix, with deflations whenever a
   subdiagonal becomes zero.
3. **Post-processing.** The signs are normalized so that `σ >= 0` and the
   singular values are sorted in descending order.

The reduction costs `O(m · n²)` and the QR iteration costs `O(n · (m + n))` per
sweep; in practice the decomposition converges in a few sweeps.

#### The product `v vᵀ` is never materialized

The Householder stage needs `v vᵀ`, but that product is **not** built. Since
`v` is a vector (not a matrix), the transpose `vᵀ` is just an inner product:
the reflector is applied column by column (or row by row) using `vᵀ · x`
computed directly, without allocating intermediate matrices.

This avoids allocating an `m x n` matrix per reflector (the previous version
of the code created `matrix_transpose` and `matrix_mul` just to obtain
`v vᵀ`) and replaces the `O(m · n²)` cost of the full multiplication with the
`O(m · n)` cost of the reflector applied directly. In other words: **there is
no "transposed vector" in the library — vector transposition is implicit in
the inner product**, which is the efficient way to compute `vᵀ · x`.

#### Example

```c
#include "matrix.h"
#include "svd.h"

int main(void) {
  double c1[] = { 1.0, 3.0 };
  double c2[] = { 2.0, 4.0 };
  Vector columns[] = {
    { .dim = 2, .x = c1 },
    { .dim = 2, .x = c2 },
  };

  Matrix* a = matrix_create(2, 2, columns);
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

#### Intermediate stages

Besides `svd_golub_reinsch`, the three stages are public and can be used
separately. `svd_qr_bidiag_iter` takes ownership of `d`, `U` and `V` from the
received `HouseholderBidiagResult`, and `svd_pos_proc` takes ownership of `d`
and `U` from the `QRBidiagResult`; in both cases the caller frees the
intermediate structure with the corresponding `*_destroy` function:

```c
HouseholderBidiagResult* hbdr = svd_householder_bidiag(a);
QRBidiagResult* bdr = svd_qr_bidiag_iter(hbdr);
svd_householder_bidiag_destroy(hbdr);

SVDResult* r = svd_pos_proc(bdr);
svd_qr_bidiag_destroy(bdr);
```

### Linear systems (`lsystem.h`)

| Function | Description |
| --- | --- |
| `lsystem_solve(a, b)` | Solves the linear system `A · x = b` |

The method is chosen automatically:

- **square and non-singular** matrix: LU decomposition with partial pivoting
  (forward and backward substitution);
- **singular** or **rectangular** matrix (overdetermined or underdetermined):
  least squares via SVD with a cutoff for singular values below
  `eps · max(m, n) · σ₁` — the returned solution is the one with **minimum
  norm**.

```c
double c1[] = { 3.0, 2.0 };
double c2[] = { -1.0, 4.0 };
Vector columns[] = {
  { .dim = 2, .x = c1 },
  { .dim = 2, .x = c2 },
};

Matrix* a = matrix_create(2, 2, columns);
double bx[] = { 4.0, 26.0 };
Vector* b = vector_create(2, bx);

Vector* x = lsystem_solve(a, b); // x = (3, 5)
```

## Example

```c
#include "matrix.h"

int main(void) {
  double c1[] = { 1.0, 3.0 };
  double c2[] = { 2.0, 4.0 };
  Vector columns[] = {
    { .dim = 2, .x = c1 },
    { .dim = 2, .x = c2 },
  };

  // A = | 1  2 |
  //     | 3  4 |
  Matrix* a = matrix_create(2, 2, columns);
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

## Error handling

On failure, the functions return `NULL` (or `NAN`, in the case of
`matrix_get`, `vector_inner_product` and the norm functions) and set `errno`
to one of the codes below:

| Code | Value | Meaning |
| --- | --- | --- |
| `EINVALIDDIM` | 1000 | Invalid or incompatible vector dimension |
| `EVECTORALLOCATION` | 1001 | Failed to allocate a vector |
| `EVECTORNULL` | 1002 | Null vector passed as an argument |
| `EINVALIDP` | 1003 | Invalid order `p` for the vector norm |
| `EMATRIXALLOCATION` | 2000 | Failed to allocate a matrix |
| `EMATRIXINVALIDADDDIM` | 2001 | Incompatible dimensions for addition/subtraction |
| `EMATRIXINVALIDLINEINDEX` | 2002 | Row index out of bounds |
| `EMATRIXINVALIDCOLINDEX` | 2003 | Column index out of bounds |
| `EMATRIXINVALIDMUL` | 2004 | Incompatible dimensions for multiplication |
| `EMATRIXNULL` | 2999 | Null matrix passed as an argument |
| `EMATRIXINVALIDDIM` | 2005 | Non-square matrix or invalid dimension |
| `EMATRIXSINGULAR` | 2006 | Singular matrix, no unique LU decomposition |
| `EMATRIXINVALIDDETERMINANT` | 2007 | Failed to compute the determinant |
| `ESVDNULL` | 4000 | Null argument in the singular value decomposition |
| `ESVDINVALIDDIM` | 4001 | Matrix with no rows or no columns |
| `ESVDALLOCATION` | 4002 | Memory allocation failure in the decomposition |
| `ESVDNOTCONVERGED` | 4003 | QR iteration did not converge within the sweep limit |
| `ELSYSTEMNULL` | 3000 | Null argument in the linear system solve |
| `ELSYSTEMINVALIDDIM` | 3001 | Incompatible dimensions between `A` and `b`, or zero dimensions |

## Next steps

- Separate the tests from the main executable
- Use the expanded range (n up to ~2000) as the basis for the next libraries:
  statistics (`lstat`), probability (`lprob`) and machine learning (`lml`)
