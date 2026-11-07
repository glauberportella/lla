# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.4.0] - 2026-11-07

### Added

- Single English documentation set (`README.md`, `examples/README.md` and
  this changelog) for broader visibility.
- `CHANGELOG.md` documenting the project history.

### Changed

- Dimensions, indices, sizes and counts migrated from `unsigned int` to
  `size_t` across the public API and the implementation:
  - `Vector.dim`;
  - `Matrix.nlins` and `Matrix.ncols`;
  - creation, access, index validation, column manipulation, cross product
    and matrix operations;
  - `LU.piv`, which changed from `int*` to `size_t*`.
- `matrix_rank` now returns `size_t` and returns `SIZE_MAX` with `errno` set
  on error.
- Size and index printing updated to `%zu` in the benchmark, `lim`, tests and
  examples.
- `unsigned int` kept only where the value is not a size: the order `p` of
  `vector_norm_p` and the 32-bit output of the pseudo-random generators.

### Fixed

- Duplicated `b == NULL` check in `svd_householder_bidiag`.
- Potential underflow in `matrix_remove_column_vector_at` when `ncols == 0`.
- `int`-based reverse loop in the LU example's back substitution.

### Migration notes

- This release breaks API/ABI compatibility: struct fields, several function
  signatures and the return type/error sentinel of `matrix_rank` changed.

## [0.3.0] - 2026-10-28

### Added

- Optional OpenMP build (`make OPENMP=1`) to parallelize `matrix_mul` over the
  columns of the result.
- `matrix_copy`, returning a deep copy of the matrix.

### Changed

- Matrices are now stored in a single contiguous block of doubles, in
  column-major layout, with one allocation per matrix.
- Internal operations access the block directly through pointers; index
  validation is restricted to the public API (`matrix_get`/`matrix_set`).
- `matrix_mul` reordered to the `j -> k -> i` loop (saxpy over columns).
- `matrix_inverse` rewritten using forward and backward substitution over the
  LU decomposition, without building the identity or `n x 1` matrices.
- `matrix_add_column_vector_at` now copies the values and no longer takes
  ownership of the input buffer.
- `matrix_get_vector_column_at` now returns a copy of the column.
- Benchmark and dimension limit search updated with the new performance
  numbers.

## [0.2.0] - 2026-10-26

### Added

- Seven independent examples covering the API: vectors, matrices, determinant
  and inverse, LU decomposition, linear systems, least squares and SVD.
- Dedicated `Makefile` to build and run the examples.

### Changed

- README expanded with the examples documentation.

## [0.1.0] - 2026-10-26

### Added

- Singular value decomposition (Householder + QR iteration).
- `lsystem_solve`, with LU solving and SVD-based least squares.
- `matrix_rank` and `matrix_trace`.
- Static (`liblla.a`) and shared (`liblla.so`) library builds, with `-fPIC`
  objects in `obj/`.

## [0.0.6] - 2026-10-24

### Added

- `vector_product`, the generalized cross product.
- Tests and documentation for the cross product.

## [0.0.5] - 2026-10-21

### Added

- `src/bech.c`, a matrix operation benchmark.
- `src/lim.c`, a search for practical dimension limits.
- Benchmark, dimension limit and rule-of-thumb sections in the README.

## [0.0.4] - 2026-10-21

### Added

- Vector inner product.
- Vector norms (`p`, infinity and Euclidean).
- Optional labels when printing vectors.
- Tests for the new operations.

## [0.0.3] - 2026-10-20

### Added

- `matrix_determinant` via LU decomposition.

## [0.0.2] - 2026-10-20

### Added

- `matrix_inverse`.

## [0.0.1] - 2026-10-20

### Added

- First version of the library: `Vector`, `Matrix`, basic operations,
  printing, LU decomposition and tests.
