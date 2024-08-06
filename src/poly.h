#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

#define DEBUG_MODE 0
#define DEG_THRESHOLD 32
#define MOD UINT32_MAX
#define MAX(a,b) (((a)>(b))?(a):(b))

struct poly_u
{
    // degree
    size_t deg;
    // coeffs
    uint32_t* c;
};

typedef struct poly_u poly_u_t;

/* allocates space for a polynomial of unsigned integer coefficients */
poly_u_t* alloc_poly(size_t deg);

/* deallocates p by freeing the space occupied by p on the heap */
void free_poly_u(poly_u_t *p);

void display_poly_u(poly_u_t* p, bool inc);

bool equals(poly_u_t* r1, poly_u_t* r2);

poly_u_t* mulpu(poly_u_t* p, poly_u_t* q);

poly_u_t* mulpukr(poly_u* p, poly_u_t* q);