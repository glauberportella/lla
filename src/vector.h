#ifndef _VECTOR_H_
#define _VECTOR_H_

#define EINVALIDDIM       1000
#define EVECTORALLOCATION 1001
#define EVECTORNULL       1002
#define EINVALIDP         1003

typedef struct vector_rn_t {
  double* x;
  unsigned int dim;
} Vector;

Vector* vector_create_empty(unsigned int dim);
Vector* vector_create(unsigned int dim, double coords[]);
void vector_destroy(Vector* v);

Vector* vector_add(const Vector* a, const Vector* b);
Vector* vector_sub(const Vector* a, const Vector* b);
Vector* vector_mul(double scalar, const Vector* a);
double vector_inner_product(const Vector* a, const Vector* b);
double vector_norm_p(const Vector* a, unsigned int p);
double vector_norm_inf(const Vector* a);
double vector_norm(const Vector* a);

Vector* vector_product(unsigned int dim, const Vector vectors[]);

void vector_print(const Vector* v, const char* vector_label);
void vector_println(const Vector* v, const char* vector_label);

#endif
