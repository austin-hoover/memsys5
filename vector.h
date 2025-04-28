/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Project:   MemSys5  Quantified Maximum Entropy System.
 * 
 * Filename:  vector.h
 * 
 * Purpose:   Header for vector.c
 * 
 * History:   21 Feb 1992, 01 May 1999
 * 
 *            Copyright (c) 1992-1999 Maximum Entropy Data Consultants Ltd.
 *-----------------------------------------------------------------------------
 */
#ifndef VECTORH
#define VECTORH

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * 
 *                            Vector procedures
 *                            -----------------
 * 
 *   VecFill  (u,a,n)          u[] = a
 *   VecMov   (u,v,n)          v[] = u[]
 *   VecSign  (u,v,n)          v[] = sign(u[])
 *   VecAdd   (u,v,w,n)        w[] = u[] + v[]
 *   VecSub   (u,v,w,n)        w[] = u[] - v[]
 *   VecMul   (u,v,w,n)        w[] = u[] * v[]
 *   VecDiv   (u,v,w,n)        w[] = u[] / v[]
 *   VecRecip (u,v,n)          v[] = 1 / u[]
 *   VecLog   (u,v,n)          v[] = log(u[])
 *   VecExp   (u,v,n)          v[] = exp(u[])
 *   VecSqrt  (u,v,n)          v[] = sqrt(u[])
 *   VecSadd  (u,a,v,n)        v[] = u[] + a
 *   VecSmul  (u,a,v,n)        v[] = u[] * a
 *   VecSmula (u,a,v,w,n)      w[] = u[] * a + v[]
 *   VecSum   (u,a,n)          a = sum(u[])
 *   VecDot   (u,v,a,n)        a = u[].v[]
 *   VecNonlin(u,v,w,n,umemx)  v[] = umemx(u[]) ,  w[] = umemx'(u[])
 * 
 *   VecRand0 (krand,iseed)    Initialise random generator into  krand[0...111]
 *   VecRand  (krand,u,n)      u[] = Random normal samples using krand[0...111]
 * 
 *-----------------------------------------------------------------------------
 */

extern void VecFill  ( float*, float, int );
extern void VecMov   ( float*, float*, int );
extern void VecSwap  ( float*, float*, int );
extern void VecSign  ( float*, float*, int );
extern void VecAdd   ( float*, float*, float*, int );
extern void VecSub   ( float*, float*, float*, int );
extern void VecMul   ( float*, float*, float*, int );
extern void VecDiv   ( float*, float*, float*, int );
extern void VecRecip ( float*, float*, int );
extern void VecLog   ( float*, float*, int );
extern void VecSqrt  ( float*, float*, int );
extern void VecExp   ( float*, float*, int );
extern void VecSadd  ( float*, float, float*, int );
extern void VecSmul  ( float*, float, float*, int );
extern void VecSmula ( float*, float, float*, float*, int );
extern void VecSum   ( float*, float*, int );
extern void VecDot   ( float*, float*, float*, int );
extern int  VecNonlin( float*, float*, float*, int, int(float,float*,float*) );
extern void VecRand0 ( unsigned*, int );
extern void VecRand  ( unsigned*, float*, int );

#endif
