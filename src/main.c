#include "poly.h"
#include <bits/types/struct_timeval.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>

// #define DEG 500000

int main(int argc, char** argv)
{
    
    if (argc != 2)
    {
        printf("ERROR: Wrong number of args, usage: ./main degree \
        \nexample : ./main 50000\n");
        return 1;
    }

    struct timeval start_tv;
    struct timeval end_tv;

    srand(time(0));

    // degree of both polynomials
    int deg = atoi(argv[1]);

    poly_u_t* p = alloc_poly(deg);
    poly_u_t* q = alloc_poly(deg);


    for (size_t i = 0 ; i <= p->deg ; i++)
    {
        p->c[i] = 1 + (rand() % MOD);
    }

    for (size_t i = 0 ; i <= q->deg ; i++)
    {
        q->c[i] = 1 + (rand() % MOD);
    }

    // the result to be stored in r
    poly_u_t* r1 = alloc_poly(2*deg);
    poly_u_t* r2 = alloc_poly(2*deg);


    uint32_t n = deg/2;

    // TODO: change this later to ONLY 3 parameters:
    // 1 poly to store result + 2 polys to multiply
    // mulpu(r1, p, q, 0, p->deg, 0, q->deg);

    // printf("p(x) = ");
    // display_poly_u(p, true);
    // printf("q(x) = ");
    // display_poly_u(q, true);
    // printf("\n");

    gettimeofday(&start_tv, NULL);
    r1 = mulpu(p, q);
    // r1 = mulpu_(p, q);
    gettimeofday(&end_tv, NULL);
    // printf("p(x)q(x) = ");
    // display_poly_u(r1, true);
    printf("Time taken (naive): %ld seconds\n\n", end_tv.tv_sec - start_tv.tv_sec);

    gettimeofday(&start_tv, NULL);
    // mulpuk1(r2, p, q, deg);
    r2 = mulpukr(p, q);
    gettimeofday(&end_tv, NULL);
    // printf("p(x)q(x) = ");
    // display_poly_u(r2, true);
    printf("Time taken (mulpukr): %ld seconds\n\n", end_tv.tv_sec - start_tv.tv_sec);

    // Check correctness
    if (equals(r1, r2))
    {
        printf("[Success]: Consistent results\n");
    }
    else {
        printf("[Failure]: The results are not consistent!\n");
    }

    free_poly_u(p);
    free_poly_u(q);

    free_poly_u(r1);
    free_poly_u(r2);

    return 0;
}