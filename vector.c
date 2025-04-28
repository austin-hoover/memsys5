/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Project:   MemSys5  Quantified Maximum Entropy System.
 * 
 * Filename:  vector.c
 * 
 * Purpose:   Vector library.
 * 
 * History:   21 Feb 1992   Array-processor code
 *            01 May 1999   Normal deviates by Box-Muller/Knuth-subtractive
 * 
 *            Copyright (c) 1992-1999 Maximum Entropy Data Consultants Ltd.
 *-----------------------------------------------------------------------------
 */
#include <math.h>
#include "vector.h"

#define SQRT(x) (float)sqrt(x)
#define EXP(x)  (float)exp(x)
#define LOG(x)  (float)log(x)

#define ZERO    ((float)0.0)
#define ONE     ((float)1.0)
#define TWO     ((float)2.0)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *                               Vector Procedures
 *-----------------------------------------------------------------------------
 */

void VecFill(
float* u,
float  a,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        u[i] = a;
}

void VecMov(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = u[i];
}

void VecSign(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = (u[i] >= ZERO) ? ONE : -ONE;
}

void VecAdd(
float* u,
float* v,
float* w,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        w[i] = u[i] + v[i];
}

void VecSub(
float* u,
float* v,
float* w,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        w[i] = u[i] - v[i];
}

void VecMul(
float* u,
float* v,
float* w,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        w[i] = u[i] * v[i];
}

void VecDiv(
float* u,
float* v,
float* w,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        w[i] = u[i] / v[i];
}

void VecRecip(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = ONE / u[i];
}

void VecLog(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = LOG(u[i]);
}

void VecSqrt(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = SQRT(u[i]);
}

void VecExp(
float* u,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = EXP(u[i]);
}

void VecSadd(
float* u,
float  a,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = u[i] + a;
}

void VecSmul(
float* u,
float  a,
float* v,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        v[i] = a * u[i];
}

void VecSmula(
float* u,
float  a,
float* v,
float* w,
int    n)
{
    int   i;
    for( i = 0; i < n; ++i )
        w[i] = a * u[i] + v[i];
}

void VecSum(
float* u,
float* a,
int    n)
/* Vector sum
 *                        A := Sum u[.]
 * Sum the product as binary tree, e.g. A=((a+b)+(c+d))+((e+f)+(g+h)),
 *  effectively guaranteeing full accuracy.
 * This routine is a slower but more accurate version of the
 *  following more easily vectorisable code ................
 *      A = 0.;
 *      for( i = 0; i < n; ++i )
 *          A += u[i];
 * .........................................................
 * Allow for addresses up to 32 bits. Z stores the branch sums
 */
{
    int   i, j;
    float z[32];
    for( j = 0; j < 32; ++j)
        z[j] = ZERO;
    for( i = 0; i < n; ++i )
    {
        *a = u[i];
        j = 0;
        while (z[j] != ZERO)
        {
            *a += z[j];
            z[j] = ZERO;
            ++j;
        }
    z[j] = *a;
    }
/* Accumulate all the branches, in increasing order, */
    *a = z[0];
    for( j = 1; j < 32; ++j)
        *a += z[j];
}

void VecDot(
float* u,
float* v,
float* a,
int    n)
/* Vector dot
 *                        A := u[.].v[.]
 * Sum the product as binary tree, e.g. A=((a+b)+(c+d))+((e+f)+(g+h)),
 *  effectively guaranteeing full accuracy.
 * This routine is a slower but more accurate version of the
 *  following more easily vectorisable code ................
 *      A = 0.;
 *      for( i = 0; i < n; ++i )
 *          A += u[i] * v[i];
 * .........................................................
 * Allow for addresses up to 32 bits. Z stores the branch sums
 */
{
    int   i, j;
    float z[32];
    for( j = 0; j < 32; ++j)
        z[j] = ZERO;
    for( i = 0; i < n; ++i )
    {
        *a = u[i] * v[i];
        j = 0;
        while (z[j] != ZERO)
        {
            *a += z[j];
            z[j] = ZERO;
            ++j;
        }
        z[j] = *a;
    }
/* Accumulate all the branches, in increasing order, */
    *a = z[0];
    for( j = 1; j < 32; ++j)
        *a += z[j];
}

int  VecNonlin(     /* Negative return value transmits error code from umemx */
float* u,
float* v,
float* w,
int    n,
int    umemx(float,float*,float*))
{
    int   i;
    int   CALLvalue;
    for( i = 0; i < n; ++i )
        if( (CALLvalue = umemx(u[i], &v[i], &w[i])) < 0 )
            return CALLvalue;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  VecRand0
 * 
 * Purpose:   Initialise random generator state in krand[0...111].
 * 
 * Notes:     Knuth's subtractive method
 * 
 * History:   JS             06 May 1995     First release
 *                           19 Oct 1995     MassInf version
 *                           01 May 1999     MaxEnt version
 *-----------------------------------------------------------------------------
 */
void  VecRand0(
unsigned* krand,       /*   O  Random generator state in [0...111] */
int       iseed)       /* I    Seed */
{
    unsigned    i, j, k, m;

    m = (unsigned)(161803398 - iseed);
    krand[54] = m;
    krand[109] = 0;
    k = 1;
    for( i = 1; i < 55; ++i )
    {
        j = (21 * i - 1) % 55;
        krand[j] = k;
        k = m - k;
        m = krand[j];
        krand[54+i] = i;
    }
    for( k = 0; k < 4; ++k )
        for( i = 0; i < 55; ++i )
            krand[i] -= krand[(i + 31) % 55];
    krand[110] = 0;
    krand[111] = 31;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  VecRand
 * 
 * Purpose:   Normal random number generator. 
 *            u[0],u[1],...,u[n-1]  :=  samples from N(0,1)
 * 
 * Notes: (1) Requires previous call to VecRand0 to initialise the generator
 *            stored in krand[0..111]
 * 
 * History:   JS         22 Jan 1994   Box-Muller/Knuth-subtractive
 *                       19 Oct 1995   MassInf version
 *                       01 May 1999   Vector version
 *-----------------------------------------------------------------------------
 */
void VecRand(
unsigned* krand,   /* I O  Random generator state in [0...111] */
float*    u,       /* I O  Vector arithmetic workspace */
int       n)       /* I    Length of vector */
{
    static const double Z = (double)((unsigned)(-1) >> 1);
    unsigned   j, k;
    int   i, m;
    float a, r, g1, g2;
    float* v;

    m = (n + 1) / 2;
    v = u + n / 2;
    for( i = 0; i < m; ++i )
    {
/* Generate 2 normal deviates (g1*a) and (g2*a) */
        do
        {
            j = krand[110];    krand[110] = krand[55+j];
            k = krand[111];    krand[111] = krand[55+k];
            krand[j] -= krand[k];
            g1 = (float) (krand[j] / Z - ONE);                               /* (-1,1) */
            j = krand[110];    krand[110] = krand[55+j];
            k = krand[111];    krand[111] = krand[55+k];
            krand[j] -= krand[k];
            g2 = (float) (krand[j] / Z - ONE);                               /* (-1,1) */
            r = g1 * g1 + g2 * g2;
        } while( r >= ONE );
        a = (float) (sqrt(-TWO * log(r) / r));
        u[i] = g1 * a;
        v[i] = g2 * a;
    }
}
