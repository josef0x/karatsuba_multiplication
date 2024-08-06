#include "poly.h"
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

poly_u_t* alloc_poly(size_t deg)
{
    assert(deg < UINT32_MAX - 1);

    poly_u_t* p = malloc(sizeof(poly_u_t));

    if (p == NULL)
    {
        printf("[ERROR] malloc() failed!\n");
        return NULL;
    }

    p->deg = deg;
    p->c = malloc((deg + 1) * sizeof(uint32_t));
    
    if (p->c == NULL)
    {
        printf("[ERROR] malloc() failed!\n");
        return NULL;
    }

    // intialize the polynomial to zero
    //memset(p->c, 0, deg+1);
    for (size_t i = 0 ; i <= p->deg ; i++)
    {
        p->c[i] = 0;
    }

    return p;
}

void free_poly_u(poly_u_t *p)
{
    free(p->c);
    free(p);
}

void display_poly_u(poly_u_t* p, bool inc)
{
    if (p == NULL)
    {
        printf("[ERROR] invalid polynomial\n");
        return;
    }

    if (inc)
    {
        // display the coefficients in increasing order
        for (size_t i = 0 ; i <= p->deg ; i++)
        {
            printf("%u ", p->c[i]);
        }
    }
    else
    {
        // display the coefficients in decreasing order
        for (size_t i = 0 ; i <= p->deg ; i++)
        {
            printf("%u ", p->c[p->deg-i]);
        }
    }
    printf("\n");
}

bool equals(poly_u_t* r1, poly_u_t* r2)
{
    if (r1->deg != r2->deg)
    {
        return false;
    }

    for (size_t i = 0 ; i <= r1->deg ; i++)
    {
        if (r1->c[i] != r2->c[i])
        {
            printf("{%lu}\n", i);
            return false;
        }
    }
    return true;
}

poly_u_t* mulpu(poly_u_t* p, poly_u_t* q)
{
    if (p == NULL || q == NULL) return NULL;
    
    uint32_t result_deg = p->deg + q->deg;
    poly_u_t* result = alloc_poly(result_deg);
    if (result == NULL) return NULL;

    for (uint32_t i = 0; i <= p->deg; i++) {
        for (uint32_t j = 0; j <= q->deg; j++) {
            result->c[i + j] += p->c[i] * q->c[j];
        }
    }

    return result;
}

poly_u_t* mulpukr(poly_u_t* p, poly_u_t* q)
{

    if (p == NULL || q == NULL || p->deg != q->deg) return NULL;
    
    if (p->deg <= 100)
    {
        return mulpu(p, q);
    }
    
    uint32_t n = p->deg;
    uint32_t half_n = n/2;
    
    poly_u_t* p0 = alloc_poly(half_n - 1);
    poly_u_t* p1 = alloc_poly(n - half_n);
    poly_u_t* q0 = alloc_poly(half_n - 1);
    poly_u_t* q1 = alloc_poly(n - half_n);
    
    // Split p and q into p0, p1, q0, q1
    for (uint32_t i = 0; i < half_n; i++)
    {
        p0->c[i] = p->c[i];
        q0->c[i] = q->c[i];
    }
    for (uint32_t i = half_n; i <= n; i++)
    {
        p1->c[i - half_n] = p->c[i];
        q1->c[i - half_n] = q->c[i];
    }
    
    poly_u_t* A = mulpukr(p0, q0);
    poly_u_t* C = mulpukr(p1, q1);
    
    poly_u_t* p0_plus_p1 = alloc_poly(MAX(n - half_n, half_n - 1));
    poly_u_t* q0_plus_q1 = alloc_poly(MAX(n - half_n, half_n - 1));

    for (uint32_t i = 0; i <= MAX(n - half_n, half_n - 1); i++)
    {
        // p0_plus_p1->c[i] = p0->c[i] + p1->c[i];
        // q0_plus_q1->c[i] = q0->c[i] + q1->c[i];

        if (i <= p0->deg)
        {
            p0_plus_p1->c[i] += p0->c[i];
        }
        if (i <= p1->deg)
        {
            p0_plus_p1->c[i] += p1->c[i];
        }
        if (i <= q0->deg)
        {
            q0_plus_q1->c[i] += q0->c[i];
        }
        if (i <= q1->deg)
        {
            q0_plus_q1->c[i] += q1->c[i];
        }
    }

    poly_u_t* B = mulpukr(p0_plus_p1, q0_plus_q1);

    size_t i = B->deg;
    while (B->c[i--] == 0)
    {
        B->deg = i;
    }

    poly_u_t* result = alloc_poly(2*n);

    // Combine A, B, C to form the result
    for (uint32_t i = 0; i <= A->deg; i++)
    {
        result->c[i] += A->c[i];
    }

    for (uint32_t i = 0; i <= B->deg; i++)
    {
        result->c[i + half_n] += B->c[i];

        if (i <= A->deg)
            result->c[i + half_n] -= A->c[i];

        if (i <= C->deg)
            result->c[i + half_n] -= C->c[i];

    }
    for (uint32_t i = 0; i <= C->deg; i++) {
        result->c[i + 2*half_n] += C->c[i];
    }

    // Free temporary polynomials
    free_poly_u(p0); free_poly_u(p1); free_poly_u(q0); free_poly_u(q1);
    free_poly_u(A); free_poly_u(B); free_poly_u(C);
    free_poly_u(p0_plus_p1); free_poly_u(q0_plus_q1);
    
    return result;
}