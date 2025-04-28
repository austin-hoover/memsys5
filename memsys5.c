/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Project:   MemSys5  Quantified Maximum Entropy System.
 *
 * Filename:  memsys5.c
 *
 * Purpose:   Kernel functions.
 *
 * History:   Conceptualised and programmed by
 *            Stephen Gull, John Skilling, and Mark Charter
 *
 *            Version 1.41                01 May 1999
 *
 *            Copyright (c) 1990-1999
 *            Maximum Entropy Data Consultants Ltd
 *            115c Milton Road, Cambridge CB4 1XE, England
 *-----------------------------------------------------------------------------
 */

#define FLOWCHART 0   /* 0 = off, 1 = print major procedure names */
#define CONTROL   0   /* 0 = off, 1 = extra numerical diagnostics */

/***************/
/* Environment */
/***************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "memsys5.h"
#include "vector.h"

/**************************/
/* Prototype declarations */
/**************************/

static int  MemOp   ( const MEMSYS*, int, int, short );
static int  MemTr   ( const MEMSYS*, int, int, short );
static int  MemQuant( const MEMSYS*, float, float, short, float*, float*,
                                     short*, short* );
static int  MemCgh  ( const MEMSYS*, float, float, short, float*, float*,
                                     short*, int*, float*, float*, float* );
static int  MemCgd  ( const MEMSYS*, float, float, float, float, short, float*,
                                     float*, float*, float*, float*,
                                     int*, int*, short*, short* );
static int  MemDet  (                float*, float*, float*, int, int, int );
static int  MemLa0  ( const MEMSYS*, float, float, float, float, float,
                                     float*, short* );
static int  MemLa   ( const MEMSYS*, float*, float, float, float, float, float,
                                     short, short, short*, short );
static int  MemLaw  ( const MEMSYS*, float, float, float, short );
static int  MemLar  ( const MEMSYS*, float, float*, float* );
static int  MemLb   (                float, float, float, float, float,
                                     float*, short* );
static int  MemZero ( const MEMSYS*, int );
static int  MemSmul ( const MEMSYS*, int, float, int );
static int  MemVmul ( const MEMSYS*, int, int, int );
static int  MemDot  ( const MEMSYS*, int, int, float* );
static int  MemCopy ( const MEMSYS*, int, int );
static int  MemVrand( const MEMSYS*, int, int, int );
static int  MemSma  ( const MEMSYS*, int, float, int, int );
static int  MemSmd  ( const MEMSYS*, int, float, int, float* );
static int  MemEnt  ( const MEMSYS*, float*, float*, float* );
static int  MemEnt1 ( const MEMSYS*, float, float*, float*, float* );
static int  MemEnt2 ( const MEMSYS*, float, float*, float*, float* );
static int  MemEnt3 ( const MEMSYS*, float, float*, float*, float* );
static int  MemEnt4 ( const MEMSYS*, float, float*, float*, float* );
static int  MemEnt51( const MEMSYS*, float, float*, float* );
static int  MemEnt52( const MEMSYS*, float, float,
                                     float*, float*, float* );
static int  MemProj ( const MEMSYS*, int, int );
static int  MemChi  ( const MEMSYS*, float*, float*, float* );
static int  MemPoi  ( const MEMSYS*, float*, float*, float* );
static int  MemTest ( const MEMSYS*, float* );

/************************/
/* Constants and macros */
/************************/

/* Memory handling */
#define CALLOC(p,n,t) {if((n)>0) if(!(p=(t*)calloc((size_t)(n),sizeof(t)))) return E_MALLOC;}
#define FREE(p)       {if(p) (void)free( (void*)p ); p=NULL;}

/* Error exit handling */
#define CALL(x)   {if( (CALLvalue = (x)) < 0 ) return CALLvalue;}

/* Diagnostic handling */
#define FPRINTF   (void)fprintf
#define SPRINTF   (void)sprintf

/* Arithmetic */
#define ABS(x)    ((x) >=  0  ? (x) :-(x))
#define MIN(x,y)  ((x) <= (y) ? (x) : (y))
#define MAX(x,y)  ((x) >= (y) ? (x) : (y))
#define SQRT(x)   ((float)sqrt(x))
#define EXP(x)    ((float)exp(x))
#define LOG(x)    ((float)log(x))
#define COS(x)    ((float)cos(x))
#define ATAN(x)   ((float)atan(x))

#define ZERO      ((float)0.0)
#define HALF      ((float)0.5)
#define ONE       ((float)1.0)
#define TWO       ((float)2.0)
#define THREE     ((float)3.0)
#define FOUR      ((float)4.0)
#define EIGHT     ((float)8.0)
#define NINE      ((float)9.0)
#define TWELVE    ((float)12.0)
#define TWENTY7   ((float)27.0)
#define FIFTY4    ((float)54.0)
#define TINY      ((float)1.0e-30)
#define EPS       ((float)1.0e-20)
#define SMALL     ((float)1.0e-6)
#define BIG       ((float)1.0e20)

#define TRUE      1
#define FALSE     0

/* Sizes */
#define LMAX     30   /* Maximum number of conjugate gradient iterations */
#define NSIZE     8   /* Maximum size of (alpha, omega) table */


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: Mem5
 *
 * Purpose:  Perform one iterate of data-space maximum entropy, update
 *           the Lagrange multipliers and find the new distribution,
 *           with statistical information such as
 *                  Evidence = log Prob(Data)     (base e)
 * Methods:
 *           iEntropy = 1 is "Standard entropy"
 *           iEntropy = 2 is "Positive/negative entropy"
 *           iEntropy = 3 is "Fermi-Dirac entropy" 
 *           iEntropy = 4 is "Quadratic regularisation"
 *           iEntropy = 5 is "Standard entropy: fixed sum"
 *  
 *           iBayes   = 1 is "Classic with known noise"
 *           iBayes   = 2 is "Classic with unknown noise"
 *           iBayes   = 3 is "Alpha given"
 *           iBayes   = 4 is "Historic Chisquared=Ndata"
 * 
 *           iGauss   = 1 is "Gaussian likelihood"
 *           iGauss   = 2 is "Poisson likelihood"
 *
 *           pfnUmemx = NULL is "linear data"
 *                    else is "nonlinear data"
 *
 * Areas:
 *   Enter with:
 *     Default   -->  <3> = default model                         (if def <= 0)
 *     Data      --> <21> = data D
 *     Acc       --> <22> = accuracies = 1/sigma      (for Gaussian if acc > 0)
 *                       or Poisson accuracies      (for Poisson if memrun > 1)
 *     Lagrange  --> <23> = Lagrange multipliers w              (if memrun > 1)
 *     Nonlinear --> <29> = Nonlinearities       (if nonlinear, and memrun > 1)
 * 
 *   Exit with:
 *     Result    <--  <1> = new MaxEnt distribution f
 *     Acc       <-- <22> = Poisson accuracies                     (if Poisson)
 *     Lagrange  <-- <23> = new Lagrange multipliers
 *     Residuals <-- <24> = normalised residuals = (D-Transform(old f))/sigma
 *     Mock      <-- <27> = new mock data = Transform(new f)    (if memrun = 4)
 *     Nonlinear <-- <29> = Nonlinearities                       (if nonlinear)
 * 
 *   Full list:
 *        <1>       O       New MaxEnt distribution f
 *        <2>      (O)      MaxEnt workspace
 *        <3>     I         Default model m (if def <= 0)
 *       <21>     I         Data D
 *       <22>     I(O)      Accuracies = 1/sigma (if acc <= 0)
 *                            (Poisson = sqrt(D+1)/Transform(f))
 *       <23>     I O       Lagrange multipliers w  (not input if memrun = 1)
 *       <24>       O       Normalised residuals = (D-Transform(old f))/sigma
 *       <25>      (O)      Workspace
 *       <26>      (O)      Workspace
 *       <27>       O       Transform(new f) if memrun = 4
 *       <28>      (O)      Workspace
 *       <29>    (I O)      Rescaled accuracies (if nonlinear)
 *  
 * History:           19 Dec 1990     First release
 *                    04 Feb 1991     MemProb
 *                    13 Feb 1991     MemChi, MemPoi zeros in <29>
 *                    17 Feb 1991     n <= 1 return in MemDet fixed
 *                    20 Feb 1991     iEntropy=5 normalised option
 *                    08 Mar 1991     Revised routine MemLaw
 *                    02 May 1991     Revised output formats
 *                    22 May 1991     x**third -> exp(log(x)/three)
 *                    03 Oct 1991     log[10] for probabilities (now log[e])
 *                    09 Oct 1991     Negative rate switch
 *                    27 Nov 1991     Mock data in <27> for memrun=4
 *                    13 Jan 1992     Version 1.25 in 'C'
 *                    21 Feb 1992     "ulong" vector addressing
 *                    02 Jul 1992     #define "*H" non-recursive header files
 *                    08 Aug 1992     Bugfix in setting var (historic)
 *                    27 Apr 1999     v1.40 full revision for flat memory
 *                    01 May 1999     v1.41 ICF withdrawn
 * 
 * Returns:  Subsidiary error if present,
 *           else
 *            (64*code6+32*code5+16*code4+8*code3+4*code2+2*code1+code0)
 *           where
 *            code6 = 0 is "conjugate gradient returns GOOD to within UTOL"
 *            code6 = 1 is "conjgrad failure to maximise accurately"
 *            code5 = 0 is "conjugate gradient returns dw,H,|dr| within UTOL"
 *            code5 = 1 is "conjgrad failure to maximise accurately"
 *            code4 = 0 is "test <= |rate| / ( 1 + |rate|)"
 *            code4 = 1 is "test >  |rate| / ( 1 + |rate|) , unacceptable"
 *            code3 = 0 is "Alpha control is satisfied by current alpha"
 *            code3 = 1 is "Alpha control changes alpha"
 *            code2 = 0 is "Beta control does not impose distance penalty"
 *            code2 = 1 is "Beta control imposes distance penalty"
 *            code1 = 0 is "stopping criterion omega = 1 reached within UTOL"
 *            code1 = 1 is "stopping criterion on omega not reached"
 *            code0 = 0 is "bubble overlaps true bubble(alpha) to within UTOL"
 *            code0 = 1 is "bubble not overlapping true bubble(alpha)"
 *           but return 0 if Omega > 1 on setup.
 *-----------------------------------------------------------------------------
 */
int  Mem5(
const MEMSYS* psMemsys, 
int     memrun)            /* I    1 is  start,  > 1 is continue */
{
    int       ncell;       /* I    # result cells */
    int       ndata;       /* I    # data */

    int       mbayes;      /* I    Stopping switch */
    int       mentrp;      /* I    Entropy switch */
    int       mgauss;      /* I    Likelihood switch */
    int       nrand;       /* I    # random vecs */
    int       iseed;       /* I    Regulariser seed */
    float     aim;         /* I    Constraint target */
    float     rate;        /* I    Dist & beta control */
    float     def;         /* I    Default (if >0)    */
    float     acc;         /* I    Accuracies (if >0) */
    float     utol;        /* I    Tolerance, O(0.1) */

    float     s;           /*   O  Entropy */
    float     test;        /*   O  Gradient misfit */
    float     chisq;       /*   O  Chisquared */
    float     scale;       /*   O  Scaling sigma(D) */
    float     plow;        /*   O  < evidence */
    float     phigh;       /*   O  > evidence */
    float     pdev;        /*   O  Std dev evidence */
    float     glow;        /*   O  < good */
    float     ghigh;       /*   O  > good */
    float     gdev;        /*   O  Std dev good */
    float     omega;       /*   O  Termination */

    float     alpha;       /* I O  Regularising parm */
    float     hhigh;       /* I O  Bubble mismatch */
    int       ntrans;      /* I O  # of transforms */
    int       istat;       /* I O  Mem5 return code */
    unsigned* krand;       /* I O  Randomiser state */

    int     (*Umemx)(float,float*,float*); /* Nonlinearity */
    int     (*Uinfo)(FILE*,const char*);   /* Diagnostics */
    FILE*     LogFile;                     /* Diagnostics */

    static  float gam[LMAX+1];
    static  float del[LMAX+1];
    static  float off[LMAX+1];
    static  float vec[LMAX+1];
    static  float val[LMAX+1];
    float   agrads;
    float   alfs;
    float   alf2s;
    float   alfgd;
    float   alhood;
    float   data;
    float   ddev;
    float   detl;
    float   detlhi;
    float   detllo;
    float   detdev;
    float   dist;
    float   dhigh;
    float   dlow;
    float   good;
    float   goodhi;
    float   goodlo;
    float   goodev;
    float   grads;
    float   gsnew;
    float   h;
    float   hlow;
    float   probl;
    float   snew;
    float   smnew;
    float   summet;
    float   tol2;
    float   units;
    float   var;
    float   weight;
    float   w;
    float   x;
    float   y;

    int     i;
    int     j;
    int     m;
    int     n;
    short   lg;
    short   bcode0;
    short   bcode1;
    short   bcode2;
    short   bcode3;
    short   bcode4;
    short   bcode5;
    short   bcode6;
    char    info[80];
    int     CALLvalue = 0;

    mbayes  = psMemsys->psInput->iBayes;
    mentrp  = psMemsys->psInput->iEntropy;
    mgauss  = psMemsys->psInput->iGauss;
    nrand   = psMemsys->psInput->iNrand;
    iseed   = psMemsys->psInput->iIseed;
    aim     = psMemsys->psInput->fAim;
    rate    = psMemsys->psInput->fRate;
    def     = psMemsys->psInput->fDef;
    acc     = psMemsys->psInput->fAcc;
    utol    = psMemsys->psInput->fUtol;
    krand   = psMemsys->psOutput->iRand;
    Umemx   = psMemsys->psProcs->pfnUmemx;
    Uinfo   = psMemsys->psProcs->pfnUinfo;
    LogFile = psMemsys->fpLog;
    ncell   = psMemsys->psVpoint->iNcell;
    ndata   = psMemsys->psVpoint->iNdata;

#if FLOWCHART
    printf("  MemSys5 Version 1.41  01 May 1999   Method%2d%2d%2d%2d\n",
              Umemx ? 2 : 1, mgauss, mentrp, mbayes);
#endif
    if( memrun < 1 || 4 < memrun )
    {
        FPRINTF(stderr, " Illegal memrun value\n");
        return E_MEM_INPUT;
    }
    if( mbayes < 1 || 4 < mbayes ) 
    {
        FPRINTF(stderr, " Illegal Bayes value\n");
        return E_MEM_INPUT;
    }
    if( mentrp < 1 || 5 < mentrp )
    {
        FPRINTF(stderr, " Illegal Entropy value\n");
        return E_MEM_INPUT;
    }
    if( mgauss < 1 || 2 < mgauss )
    {
        FPRINTF(stderr, " Illegal Gauss value\n");
        return E_MEM_INPUT;
    }
    if( mgauss == 2 && mbayes == 2 )
    {
        FPRINTF(stderr, " Method incompatible with Poisson\n");
        return E_MEM_INPUT;
    }
    if( (mbayes == 1 || mbayes == 2) && nrand <= 0 ) 
    {
        FPRINTF(stderr, " Illegal Nrand value\n");
        return E_MEM_INPUT;
    }
    if( aim < ZERO )
    {
        FPRINTF(stderr, " Illegal Aim value\n");
        return E_MEM_INPUT;
    }
    if( rate == ZERO )
    {
        FPRINTF(stderr, " Illegal Rate value\n");
        return E_MEM_INPUT;
    }
    if( utol <= ZERO || ONE <= utol ) 
    {
        FPRINTF(stderr, " Illegal Utol value\n");
        return E_MEM_INPUT;
    }

/* Setup */
    for( i = 0; i <= 40; ++i )
    {
        psMemsys->psOutput->afSt [i] = NULL;
        psMemsys->psOutput->aiLen[i] = 0;
    }
    psMemsys->psOutput->afSt [1]  = psMemsys->psVpoint->afResult;
    psMemsys->psOutput->aiLen[1]  = ncell;
    CALLOC( psMemsys->psOutput->afSt[2], ncell, float )
    psMemsys->psOutput->aiLen[2]  = ncell;
    if( def <= ZERO )
    {
        psMemsys->psOutput->afSt [3]  = psMemsys->psVpoint->afDefault;
        psMemsys->psOutput->aiLen[3]  = ncell;
    }
    psMemsys->psOutput->afSt [21] = psMemsys->psVpoint->afData;
    psMemsys->psOutput->aiLen[21] = ndata;
    if( mgauss == 2 || acc <= ZERO )
    {
        psMemsys->psOutput->afSt [22] = psMemsys->psVpoint->afAcc;
        psMemsys->psOutput->aiLen[22] = ndata;
    }
    psMemsys->psOutput->afSt [23] = psMemsys->psVpoint->afLagrange;
    psMemsys->psOutput->aiLen[23] = ndata;
    psMemsys->psOutput->afSt [24] = psMemsys->psVpoint->afResiduals;
    psMemsys->psOutput->aiLen[24] = ndata;
    CALLOC( psMemsys->psOutput->afSt[25], ndata, float )
    psMemsys->psOutput->aiLen[25] = ndata;
    CALLOC( psMemsys->psOutput->afSt[26], ndata, float )
    psMemsys->psOutput->aiLen[26] = ndata;
    if( memrun == 4 )
    {
        psMemsys->psOutput->afSt [27] = psMemsys->psVpoint->afMock;
        psMemsys->psOutput->aiLen[27] = ndata;
    }
    else
    {
        CALLOC( psMemsys->psOutput->afSt[27], ndata, float )
        psMemsys->psOutput->aiLen[27] = ndata;
    }
    CALLOC( psMemsys->psOutput->afSt[28], ndata, float )
    psMemsys->psOutput->aiLen[28] = ndata;
    if( Umemx )
    {
        psMemsys->psOutput->afSt [29] = psMemsys->psVpoint->afNonlinear;
        psMemsys->psOutput->aiLen[29] = ndata;
    }

/* Proceed... */
    if( memrun == 1 ) 
    {
        istat = ntrans = 0;
        CALL( MemZero(psMemsys, 23) )
        if( mgauss == 2 ) 
            CALL( MemZero(psMemsys, 22) )
        if( Umemx ) 
            CALL( MemZero(psMemsys, 29) )
        CALL( MemZero(psMemsys, 1) )
        CALL( MemEnt (psMemsys, &s, &summet, &grads) )
        s = ZERO;
        grads = ZERO;
    }
    else 
    {
        alpha  = psMemsys->psOutput->fAlpha;
        hhigh  = psMemsys->psOutput->fHhigh;
        ntrans = psMemsys->psOutput->iNtrans;
        istat  = psMemsys->psOutput->iIstat;
        CALL( MemTr (psMemsys, 23, 1, FALSE) )
        ++ntrans;
        CALL( MemEnt(psMemsys, &s, &summet, &grads) )
    }
    CALL( MemOp  (psMemsys, 2, 25, FALSE) )
    ++ntrans;
    CALL( MemCopy(psMemsys, 25, 27) )
    if( mgauss == 1 ) 
        CALL( MemChi(psMemsys, &alhood, &data, &units) )
    else if( mgauss == 2 )
        CALL( MemPoi(psMemsys, &alhood, &data, &units) )
    units -= (float) (data * HALF * log(EIGHT * ATAN(ONE)));
    if( memrun == 1 ) 
    {
        CALL( MemTr (psMemsys, 24, 2, TRUE) )
        ++ntrans;
        CALL( MemDot(psMemsys, 2, 2, &alf2s) )
        agrads = SQRT(alf2s);
        alf2s = - alf2s / TWO;
        alfs = ZERO;
        test = ZERO;
    }
    else 
    {
        agrads = alpha * grads;
        alfs = alpha * s;
        CALL( MemTest(psMemsys, &test) )
    }
    bcode4 = (test <= ABS(rate) / (ABS(rate) + ONE));
    if( mbayes == 2 )
        scale = SQRT((alhood - alfs) * TWO / data);
    else 
        scale = ONE;
    chisq = scale * scale;
    chisq = TWO * alhood / MAX(chisq, EPS);
    SPRINTF(info,
          "  Entropy === %12.4e    Test ===%7.4f    Chisq ===%12.4e\n",
          s, test, chisq);
    CALL( Uinfo(LogFile, info) )
    if( mbayes == 2 ) 
    {
        SPRINTF(info,
          "                                                 Scale ===%12.4e\n",
          scale);
        CALL( Uinfo(LogFile, info) )
    }
/* Stopping criterion OMEGA */
    bcode0 = TRUE;
    bcode6 = TRUE;
    glow = ghigh = gdev = ZERO;
    plow = phigh = pdev = ZERO;
    dlow = dhigh = ddev = ZERO;
    var = alfgd = omega = ZERO;
    if( mbayes == 1 || mbayes == 2 )
    {
        VecRand0(krand, iseed);
        if( memrun == 1 )
        {
            for( i = 0; i < nrand; ++i )
            {
                CALL( MemVrand(psMemsys, 26, 26, -1) )
                CALL( MemTr   (psMemsys, 26, 2, TRUE) )
                ++ntrans;
                CALL( MemDot  (psMemsys, 2, 2, &alfgd) )
            }
            alfgd /= (float)nrand;
        }
        else 
        {
            weight = glow = ghigh = gdev = dlow = dhigh = ddev = ZERO;
            bcode6 = FALSE;
            tol2 = SQRT(FLT_EPSILON);
            for( i = 0; i < nrand; ++i ) 
            {
                CALL( MemVrand(psMemsys, 26, 26, -1) )
                goodlo = goodev = goodhi = detllo = detdev = detlhi = ZERO;
                CALL( MemCgd(psMemsys, ONE, ONE, alpha, utol, FALSE,
                      NULL, NULL, NULL, gam, del, &m, &n, NULL, &lg) )
                ntrans += CALLvalue;
                if( m > 0 )
                {
/* Lower limit and variance for Good and LogDet */
                    for( j = 0; j < m; ++j )
                    {
                        val[j] = gam[j] * del[j];
                        off[j + 1] = - gam[j + 1] * del[j];
                    }
                    val[m] = ZERO;
                    CALL( MemDet(val, vec, off, m + 1, 1, -1) )
                    for( j = 0; j < m; ++j )
                    {
                        x = val[j] / alpha;
                        y = vec[j] * vec[j] / alpha;
                        goodlo += y / (x + ONE);
                        goodev += y / (x + ONE) * (x / (x + ONE));
                        if( x > tol2 )
                        {
                            detllo += y * LOG(x + ONE) / x;
                            detdev += y * (LOG(x + ONE) * LOG(x + ONE)) / x;
                        }
                        else 
                        {
                            detllo += y;
                            detdev += y * x;
                        }
                    }
                    x = gam[0] * del[0] * gam[0];
                    x = x * x;
                    goodlo = x * goodlo;
                    goodev = x * TWO * goodev;
                    detllo = x * detllo;
                    detdev = x * TWO * detdev;
                }
/* Upper limit for Good and LogDet */
/* (if val[m] unavailable because m>n, worst case is val(m) big, ignore) */
                val[0] = gam[0] * del[0];
                for( j = 1; j <= n; ++j )
                {
                    val[j] = gam[j] * del[j];
                    off[j] = - gam[j] * del[j - 1];
                }
                CALL( MemDet(val, vec, off, n + 1, 1, 1) )
                for( j = 0; j <= n; ++j )
                {
                    x = val[j] / alpha;
                    goodhi += vec[j] * vec[j] * x / (x + ONE);
                    if( x > tol2 )
                        detlhi += vec[j] * vec[j] * LOG(x + ONE);
                    else
                        detlhi += vec[j] * vec[j] * x;
                }
                goodhi = gam[0] * goodhi * gam[0];
                detlhi = gam[0] * detlhi * gam[0];
/* Increment totals */
                w = lg ? ONE : SMALL;
                glow  += goodlo * w;
                ghigh += goodhi * w;
                gdev  += goodev * w;
                dlow  += detllo * w;
                dhigh += detlhi * w;
                ddev  += detdev * w;
                weight += w;
                bcode6 = bcode6 || lg;
            }
            glow  /= weight;
            ghigh /= weight;
            gdev = SQRT(gdev) / weight;
            dlow  /= weight;
            dhigh /= weight;
            ddev = SQRT(ddev) / weight;
            bcode0 = hhigh / (scale * scale) <= utol * glow;
        }
        plow  = units - data * LOG(scale) - dhigh / TWO +
                         alfs / (scale * scale) - chisq / TWO;
        phigh = units - data * LOG(scale) - dlow / TWO +
                         alfs / (scale * scale) - chisq / TWO;
        pdev  = ddev / TWO;
        probl = (plow + phigh) / TWO;
        good  = (glow + ghigh) / TWO;
        detl  = (dlow + dhigh) / TWO;
#if CONTROL
        SPRINTF(info," Good   =%14.6e= [%14.6e,%14.6e] +-%14.6e\n",
            good, glow, ghigh, gdev);
        CALL( Uinfo(LogFile, info) )
        SPRINTF(info," LogDet =%14.6e= [%14.6e,%14.6e] +-%14.6e\n",
            detl, dlow, dhigh, ddev);
        CALL( Uinfo(LogFile, info) )
#endif
        SPRINTF(info,
        "  LogProb ===%13.4e    Code ===   %1d       Good  ===%12.4e\n",
            probl, bcode6 ? 0 : 1, good);
        CALL( Uinfo(LogFile, info) )
        if( memrun == 1 )
        {
            omega = -TWO * alf2s;
            omega = alfgd * scale * scale * aim / MAX(omega, EPS);
        }
        else
        {
            omega = -TWO * alfs;
            omega = good * scale * scale * aim / MAX(omega, EPS);
        }
        var = ONE - test;
        var = test / MAX(var, EPS) + utol * utol / TWELVE;
    }
    else if( mbayes == 3 )
    {
        if( memrun == 1 )
            omega = ZERO;
        else 
            omega = aim / MAX(alpha, EPS);
        var = utol * utol / TWELVE;
    } 
    else if( mbayes == 4 )
    {
        omega = data * aim / MAX(chisq, EPS);
        var = ONE - test;
        var = test / MAX(var, EPS) + utol * utol / TWELVE;
    }
    bcode1 = (omega >= ONE - utol) && (omega <= ONE + utol);
    bcode2 = ((istat >> 2) & 1) == 0;
    istat = 0;
    if( ! bcode0 ) 
        istat += 1;
    if( ! bcode1 )
        istat += 2;
    if( ! bcode4 ) 
        istat += 16;
    if( ! bcode6 ) 
        istat += 64;

    if( memrun == 1 && omega >= ONE ) 
    {
        SPRINTF(info,"  Omega   === %10.6f\n",omega);
        CALL( Uinfo(LogFile, info) )
        istat = 0;
    } 
    else if( memrun == 4 )
    {
        SPRINTF(info,
            "  Omega   === %10.6f                         Alpha ===%12.4e\n",
            omega, alpha);
        CALL( Uinfo(LogFile, info) )
    }
    else 
    {
        if( memrun == 1 )
            CALL( MemLa0(psMemsys, rate, summet, agrads, omega, var,
                         &alpha, &bcode3) )
        else if( memrun == 2 )
            CALL( MemLa(psMemsys, &alpha, rate, summet, omega, var,
                        agrads, 0, bcode2, &bcode3, bcode4) )
        else if( memrun == 3 )
            CALL( MemLa(psMemsys, &alpha, rate, summet, omega, var,
                        agrads, 1, bcode2, &bcode3, bcode4) )
        CALL( MemSma(psMemsys, 23, - alpha, 24, 26) )
        CALL( MemCgd(psMemsys, rate, summet, alpha, utol, TRUE, &hlow, &hhigh,
                     &dist, NULL, NULL, NULL, NULL, &bcode2, &bcode5) )
        ntrans += CALLvalue;
        h = (hlow + hhigh) / TWO;
#if CONTROL
        SPRINTF(info, " H      =%14.6e= [%14.6e,%14.6e]\n", h, hlow, hhigh);
        CALL( Uinfo(LogFile, info) )
#endif
        SPRINTF(info,
            "  Omega   === %10.6f      dist ===%7.4f    Alpha ===%12.4e\n",
            omega, dist, alpha);
        CALL( Uinfo(LogFile, info) )
        if( ! bcode2 )
            istat += 4;
        if( ! bcode3 )
            istat += 8;
        if( ! bcode5 ) 
            istat += 32;
    }

/* The next executable lines are for 'diagnostic' purposes only, and  */
/*  can be disabled to save one transform per iterate................ */
    CALL( MemTr  (psMemsys, 23, 1, FALSE) )
    ++ntrans;
    CALL( MemEnt (psMemsys, &snew, &smnew, &gsnew) )
    CALL( MemCopy(psMemsys, 2, 1) )
/* ...........................................end of diagnostic code. */
    SPRINTF(info,"  Ntrans  === %6d  %s  %1d%1d%1d%1d%1d%1d\n",
      ntrans, "                           Code  ===",     (istat>>5)&1,
      (istat>>4)&1, (istat>>3)&1, (istat>>2)&1, (istat>>1)&1, istat&1);
    CALL( Uinfo(LogFile, info) )

    psMemsys->psOutput->fS      = s;
    psMemsys->psOutput->fTest   = test;
    psMemsys->psOutput->fChisq  = chisq;
    psMemsys->psOutput->fScale  = scale;
    psMemsys->psOutput->fPlow   = plow;
    psMemsys->psOutput->fPhigh  = phigh;
    psMemsys->psOutput->fPdev   = pdev;
    psMemsys->psOutput->fGlow   = glow;
    psMemsys->psOutput->fGhigh  = ghigh;
    psMemsys->psOutput->fGdev   = gdev;
    psMemsys->psOutput->fOmega  = omega;

    psMemsys->psOutput->fAlpha  = alpha;
    psMemsys->psOutput->fHhigh  = hhigh;
    psMemsys->psOutput->iNtrans = ntrans;
    psMemsys->psOutput->iIstat  = istat;

    FREE( psMemsys->psOutput->afSt[28] )
    if( memrun != 4 )
        FREE( psMemsys->psOutput->afSt[27] )
    FREE( psMemsys->psOutput->afSt[26] )
    FREE( psMemsys->psOutput->afSt[25] )
    FREE( psMemsys->psOutput->afSt[2] )
    return istat;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemMovie5
 * 
 * Purpose:  Obtain MaxEnt sample f+df from posterior probability bubble,
 *           optionally correlated with previous sample(s).
 * 
 * Notes:(1) If( ncorr >= 0 ) then, in usual use,
 *                set random normal vector as source v
 *           else, exceptionally,
 *                use <4> as supplied externally for source v
 *       (2) The return codes istat are
 *              0 is "offset vector was estimated accurately"
 *              1 is "conjugate gradient failed in some way"
 *              2 is "offset vector has wrong magnitude"
 *              3 is "conjugate gradient failed and offset has wrong magnitude"
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         20 Feb 1991     iEntropy = 5 added
 *      JS             01 May 1999     Detailed revision (ncorr redefined)
 *-----------------------------------------------------------------------------
 */
int  MemMovie5(
const MEMSYS* psMemsys, 
int     ncorr)   /* I    Sample overlap (none if ncorr = 0) */
/* Areas:
 *   Enter with:
 *     Default   -->  <3> = default model                         (if def <= 0)
 *     Acc       --> <22> = accuracies = 1/sigma  (unless Gaussian and acc > 0)
 *     Lagrange  --> <23> = Lagrange multipliers w
 *     Nonlinear --> <29> = Nonlinearities                       (if nonlinear)
 * 
 *   Exit with:
 *     Result    <--  <1> = new MaxEnt sample f+df
 *     Mask      <--  <4> = new MaxEnt sample offset df
 * 
 *   Full list:
 *       <1>            O              MaxEnt bubble sample = f+(B**-1/2)v 
 *       <2>           (O)                      workspace
 *       <3>          I                Default model m (if def <= 0)
 *       <4>         (I)               (Source vector v, if ncorr < 0)
 *                      O              MaxEnt offset        = (B**-1/2)v
 *      <22>         (I)               [acc]    (if linear and acc <= 0)
 *      <23>          I                Multipliers w
 *      <25>           (O)                      workspace
 *      <26>           (O)                      workspace
 *      <27>           (O)                      workspace
 *      <28>           (O)                      workspace
 *      <29>         (I)               Rescaled accuracies (if nonlinear)
 */
{
    int     ncell  = psMemsys->psVpoint->iNcell;
    int     ndata  = psMemsys->psVpoint->iNdata;
    int     mentrp = psMemsys->psInput->iEntropy;
    int     mgauss = psMemsys->psInput->iGauss;
    float   def    = psMemsys->psInput->fDef;
    float   acc    = psMemsys->psInput->fAcc;
    float   utol   = psMemsys->psInput->fUtol;
    float   scale  = psMemsys->psOutput->fScale;
    float   alpha  = psMemsys->psOutput->fAlpha;
    int     ntrans = psMemsys->psOutput->iNtrans;
    float   s;
    float   grads;
    float   summet;
    float   qhi;
    float   qlo;
    int     i;
    int     jstat;
    short   bcodec;
    short   bcodet;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("  MemMovie5\n");
#endif
    for( i = 0; i <= 40; ++i )
    {
        psMemsys->psOutput->afSt [i] = NULL;
        psMemsys->psOutput->aiLen[i] = 0;
    }
    psMemsys->psOutput->afSt [1]  = psMemsys->psVpoint->afResult;
    psMemsys->psOutput->aiLen[1]  = ncell;
    CALLOC( psMemsys->psOutput->afSt[2], ncell, float )
    psMemsys->psOutput->aiLen[2]  = ncell;
    if( def <= ZERO )
    {
        psMemsys->psOutput->afSt [3]  = psMemsys->psVpoint->afDefault;
        psMemsys->psOutput->aiLen[3]  = ncell;
    }
    psMemsys->psOutput->afSt [4]  = psMemsys->psVpoint->afMask;
    psMemsys->psOutput->aiLen[4]  = ncell;
    if( mgauss == 2 || acc <= ZERO )
    {
        psMemsys->psOutput->afSt [22] = psMemsys->psVpoint->afAcc;
        psMemsys->psOutput->aiLen[22] = ndata;
    }
    psMemsys->psOutput->afSt [23] = psMemsys->psVpoint->afLagrange;
    psMemsys->psOutput->aiLen[23] = ndata;
    CALLOC( psMemsys->psOutput->afSt[24], ndata, float )
    psMemsys->psOutput->aiLen[24] = ndata;
    CALLOC( psMemsys->psOutput->afSt[25], ndata, float )
    psMemsys->psOutput->aiLen[25] = ndata;
    CALLOC( psMemsys->psOutput->afSt[26], ndata, float )
    psMemsys->psOutput->aiLen[26] = ndata;
    CALLOC( psMemsys->psOutput->afSt[27], ndata, float )
    psMemsys->psOutput->aiLen[27] = ndata;
    CALLOC( psMemsys->psOutput->afSt[28], ndata, float )
    psMemsys->psOutput->aiLen[28] = ndata;
    if( psMemsys->psProcs->pfnUmemx )
    {
        psMemsys->psOutput->afSt [29] = psMemsys->psVpoint->afNonlinear;
        psMemsys->psOutput->aiLen[29] = ndata;
    }

/* <4> := random input vector, possibly correlated between samples */
    if( ncorr >= 0 )
        CALL( MemVrand(psMemsys, 4, 2, ncorr) )
/* <1> := sqrt([metric]) */
    CALL( MemTr (psMemsys, 23, 1, FALSE) )
    ++ntrans;
    CALL( MemEnt(psMemsys, &s, &summet, &grads) )
    if( mentrp == 5 ) 
        CALL( MemProj(psMemsys, 4, 1) )
/* <4> := B**(-1/2) <4>     ,    <2> = workspace */
    CALL( MemQuant(psMemsys, alpha, utol, TRUE, &qlo, &qhi, &bcodec, &bcodet) )
    ntrans += CALLvalue;
    jstat = 0;
    if( ! bcodec )
        jstat += 1;
    if( ! bcodet )
        jstat += 2;
/*   Rescale to include final sqrt([metric]) factor from <1> */
    CALL( MemSmul(psMemsys, 2, scale, 2) )
    CALL( MemVmul(psMemsys, 1, 2, 4) )
/* <2> := f */
    CALL( MemTr (psMemsys, 23, 1, FALSE) )
    ++ntrans;
    CALL( MemEnt(psMemsys, &s, &summet, &grads) )
/* <1> := <2> + <4>  =  f + (B**-1/2) random = MaxEnt bubble sample */
    CALL( MemSma(psMemsys, 4, ONE, 2, 1) )

    psMemsys->psOutput->iNtrans = ntrans;
    FREE( psMemsys->psOutput->afSt[28] )
    FREE( psMemsys->psOutput->afSt[27] )
    FREE( psMemsys->psOutput->afSt[26] )
    FREE( psMemsys->psOutput->afSt[25] )
    FREE( psMemsys->psOutput->afSt[24] )
    FREE( psMemsys->psOutput->afSt[2] )
    return jstat;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  MemMask5
 * 
 * Purpose:   Estimate  SUM  f * (MaxEnt mask)    with error bars
 *                       i    i               i
 * 
 * Notes: (1) A call to MemMask5 may be more expensive than a MaxEnt iterate.
 *        (2) Calculation of standard deviation can be very ill-conditioned,
 *            involving inverse of B without a factor A, and is prone to
 *            serious rounding error.  The value of the return code should be
 *            checked (zero indicates correct operation of the procedure).
 *            A remedy for serious error may be full double precision
 *            operation.
 *        (3) The return codes are
 *             0 is "standard deviation was estimated accurately"
 *             1 is "standard deviation may not have been estimated accurately"
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         20 Feb 1991     iEntropy = 5 added
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
int  MemMask5(
const MEMSYS* psMemsys, 
float*   ameanx,  /*   O  mean of required sum */
float*   errlox,  /*   O  lower limit on standard deviation */
float*   errhix)  /*   O  upper limit on standard deviation */
/* Areas:
 *   Enter with:
 *     Default   -->  <3> = default model                         (if def <= 0)
 *     Mask      -->  <4> = MaxEnt quantification mask
 *     Acc       --> <22> = accuracies = 1/sigma  (unless Gaussian and acc > 0)
 *     Lagrange  --> <23> = Lagrange multipliers w
 *     Nonlinear --> <29> = Nonlinearities                       (if nonlinear)
 * 
 *   Full list:
 *       <1>           (O)                      workspace
 *       <2>           (O)                      workspace
 *       <3>          I                Default model m (if def <= 0)
 *       <4>          I                Input MaxEnt mask
 *      <22>         (I)               [acc]    (if linear and acc <= 0)
 *      <23>          I                Multipliers w
 *      <25>           (O)                      workspace
 *      <26>           (O)                      workspace
 *      <28>           (O)                      workspace
 *      <29>         (I)               Rescaled accuracies (if nonlinear)
 */
{
    int     ncell  = psMemsys->psVpoint->iNcell;
    int     ndata  = psMemsys->psVpoint->iNdata;
    int     mentrp = psMemsys->psInput->iEntropy;
    int     mgauss = psMemsys->psInput->iGauss;
    float   def    = psMemsys->psInput->fDef;
    float   acc    = psMemsys->psInput->fAcc;
    float   utol   = psMemsys->psInput->fUtol;
    float   scale  = psMemsys->psOutput->fScale;
    float   alpha  = psMemsys->psOutput->fAlpha;
    int     ntrans = psMemsys->psOutput->iNtrans;
    float   amean;
    float   s;
    float   grads;
    float   errhi;
    float   varhi;
    float   errlo;
    float   varlo;
    float   summet;
    int     i;
    int     kstat;
    short   bcodec;
    short   bcodet;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("  MemMask5\n");
#endif
    for( i = 0; i <= 40; ++i )
    {
        psMemsys->psOutput->afSt [i] = NULL;
        psMemsys->psOutput->aiLen[i] = 0;
    }
    CALLOC( psMemsys->psOutput->afSt[1], ncell, float )
    psMemsys->psOutput->aiLen[1]  = ncell;
    CALLOC( psMemsys->psOutput->afSt[2], ncell, float )
    psMemsys->psOutput->aiLen[2]  = ncell;
    if( def <= ZERO )
    {
        psMemsys->psOutput->afSt [3]  = psMemsys->psVpoint->afDefault;
        psMemsys->psOutput->aiLen[3]  = ncell;
    }
    psMemsys->psOutput->afSt [4]  = psMemsys->psVpoint->afMask;
    psMemsys->psOutput->aiLen[4]  = ncell;
    if( mgauss == 2 || acc <= ZERO )
    {
        psMemsys->psOutput->afSt [22] = psMemsys->psVpoint->afAcc;
        psMemsys->psOutput->aiLen[22] = ndata;
    }
    psMemsys->psOutput->afSt [23] = psMemsys->psVpoint->afLagrange;
    psMemsys->psOutput->aiLen[23] = ndata;
    CALLOC( psMemsys->psOutput->afSt[24], ndata, float )
    psMemsys->psOutput->aiLen[24] = ndata;
    CALLOC( psMemsys->psOutput->afSt[25], ndata, float )
    psMemsys->psOutput->aiLen[25] = ndata;
    CALLOC( psMemsys->psOutput->afSt[26], ndata, float )
    psMemsys->psOutput->aiLen[26] = ndata;
    CALLOC( psMemsys->psOutput->afSt[27], ndata, float )
    psMemsys->psOutput->aiLen[27] = ndata;
    CALLOC( psMemsys->psOutput->afSt[28], ndata, float )
    psMemsys->psOutput->aiLen[28] = ndata;
    if( psMemsys->psProcs->pfnUmemx )
    {
        psMemsys->psOutput->afSt [29] = psMemsys->psVpoint->afNonlinear;
        psMemsys->psOutput->aiLen[29] = ndata;
    }

    CALL( MemTr   (psMemsys, 23, 1, FALSE) )
    ++ntrans;
    CALL( MemEnt  (psMemsys, &s, &summet, &grads) )
    CALL( MemDot  (psMemsys, 2, 4, &amean) )
    CALL( MemVmul (psMemsys, 1, 4, 4) )
    if( mentrp == 5 ) 
        CALL( MemProj(psMemsys, 4, 1) )
    CALL( MemQuant(psMemsys, alpha, utol, FALSE,
                   &varlo, &varhi, &bcodec, &bcodet) )
    ntrans += CALLvalue;
    kstat = 0;
    if( ! bcodec ) 
        ++kstat;
    errlo = SQRT(varlo) * scale;
    errhi = SQRT(varhi) * scale;
    *ameanx = amean;
    *errlox = errlo;
    *errhix = errhi;
    psMemsys->psOutput->iNtrans = ntrans;
    FREE( psMemsys->psOutput->afSt[28] )
    FREE( psMemsys->psOutput->afSt[27] )
    FREE( psMemsys->psOutput->afSt[26] )
    FREE( psMemsys->psOutput->afSt[25] )
    FREE( psMemsys->psOutput->afSt[24] )
    FREE( psMemsys->psOutput->afSt[2] )
    FREE( psMemsys->psOutput->afSt[1] )
    return kstat;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemTrq
 * 
 * Purpose:  Check for   Opus / Tropus   inconsistencies.
 *           Dimensionless return parameter err should be O(rounding error) 
 *           if Opus and Tropus are each other's transposes.
 * Method:
 *           if(  iseed > 0 ) then
 *             <1> := random   ,   <26> := random
 *           else
 *             use <1> and <26> as set by user.
 *  
 *                  ABS(  <1> . Tropus*<26>  -  <26> . Opus*<1> ) 
 *           err := ------------------------------------------------------- 
 *                  sqrt( |<1>| |Tropus*<26>| |<26>| |Opus*<1>| ) 
 *  
 * Notes:    
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
int  MemTrq(
const MEMSYS* psMemsys, 
int     jseed,   /* I    Seed for random vector (or 0) */
float*  errx)    /*   O  Fractional error */
/* Areas:
 *      <1>          (O)             seed random vector
 *      <2>          (O)             (transform workspace)
 *     <25>          (O)             (transform workspace)
 *     <26>          (O)             seed random vector
 */
{
    int       ncell                             = psMemsys->psVpoint->iNcell;
    int       ndata                             = psMemsys->psVpoint->iNdata;
    void*     psUser                            = psMemsys->psVpoint->psUser;
    float**   afSt                              = psMemsys->psOutput->afSt;
    unsigned* krand                             = psMemsys->psOutput->iRand;
    int    (*Opus)  (void*,const float*,float*) = psMemsys->psProcs->pfnOpus;
    int    (*Tropus)(void*,const float*,float*) = psMemsys->psProcs->pfnTropus;
    float    urru;
    float    vttv;
    float    u;
    float    v;
    float    uu;
    float    vv;
    float    vru;
    float    utv;
    float    err;
    int      i;
    int      CALLvalue = 0;
#if FLOWCHART
    printf("  MemTrq\n");
#endif
    for( i = 0; i <= 40; ++i )
    {
        psMemsys->psOutput->afSt [i] = NULL;
        psMemsys->psOutput->aiLen[i] = 0;
    }
    CALLOC( psMemsys->psOutput->afSt[1], ncell, float )
    psMemsys->psOutput->aiLen[1]  = ncell;
    CALLOC( psMemsys->psOutput->afSt[2], ncell, float )
    psMemsys->psOutput->aiLen[2]  = ncell;
    CALLOC( psMemsys->psOutput->afSt[25], ndata, float )
    psMemsys->psOutput->aiLen[25] = ndata;
    CALLOC( psMemsys->psOutput->afSt[26], ndata, float )
    psMemsys->psOutput->aiLen[26] = ndata;

/* <1> := random, <26> := random */
    VecRand0(krand, jseed);
    CALL( MemVrand(psMemsys,  1,  1, -1) )
    CALL( MemVrand(psMemsys, 26, 26, -1) )
/* <2> := Tropus <26> */
/*  vv := <26>.<26> */
    CALL( MemDot    (psMemsys, 26, 26, &vv) )
    CALL( MemCopy   (psMemsys, 26, 25) )
    CALL( Tropus(psUser, afSt[25], afSt[2]) )
/*  utv := <1>.<2> */
/* vttv := <2>.<2> */
    CALL( MemDot    (psMemsys,  1,  2, &utv) )
    CALL( MemDot    (psMemsys,  2,  2, &vttv) )
/* <25> := Opus <1> */
/*   uu := <1>.<1> */
    CALL( MemDot    (psMemsys,  1,  1, &uu) )
    CALL( MemCopy   (psMemsys,  1,  2) )
    CALL( Opus  (psUser, afSt[2], afSt[25]) )
/*  vru := <26>.<25> */
/* urru := <25>.<25> */
    CALL( MemDot    (psMemsys, 25, 26, &vru) )
    CALL( MemDot    (psMemsys, 25, 25, &urru) )
/* Fractional error */
    u = SQRT(uu) * SQRT(vv);
    v = SQRT(urru) * SQRT(vttv);
    err = (utv - vru) * (utv - vru) / (u * v);
    *errx = SQRT(err);
    FREE( psMemsys->psOutput->afSt[26] )
    FREE( psMemsys->psOutput->afSt[25] )
    FREE( psMemsys->psOutput->afSt[2] )
    FREE( psMemsys->psOutput->afSt[1] )
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemOp
 * 
 * Purpose:  Differential forward transform, with or without full scaling
 *           if( flag ) then
 *                        <k> :=  [acc] * Opus * sqrt([metric]) * <j>
 *           else
 *                        <k> :=          Opus * <j>
 * 
 * Notes:    (1)  Opus must NOT change frequently between iterates
 *           (2)  The actual Opus transform is performed from <2> to <25>
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         20 Feb 1991     iEntropy = 5 added
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemOp(
const MEMSYS* psMemsys, 
int     j,       /* I    Input MaxEnt area number */
int     k,       /* I    Output data area number */
short   flag)    /* I    Opus / MockData switch */
/* Areas:
 *      <j>         I                Input MaxEnt area
 *      <k>           O              Output data area
 *      <1>        (I)                 (metric)
 *      <2>          (O)               (transform workspace)
 *     <22>         I                [acc]    (if linear)
 *     <25>          (O)               (transform workspace)
 *     <29>         I                [acc]    (if nonlinear)
 */
{
    int      mentrp                            = psMemsys->psInput->iEntropy;
    int      mgauss                            = psMemsys->psInput->iGauss;
    float    acc                               = psMemsys->psInput->fAcc;
    void*    psUser                            = psMemsys->psVpoint->psUser;
    float**  afSt                              = psMemsys->psOutput->afSt;
    int    (*Opus) (void*,const float*,float*) = psMemsys->psProcs->pfnOpus;
    int    (*Umemx)(float,float*,float*)       = psMemsys->psProcs->pfnUmemx;
    int      CALLvalue = 0;

#if FLOWCHART
    printf("   MemOp   (%2d,%2d)\n", j, k);
#endif
    CALL( MemCopy(psMemsys, j, 2) )
    if( flag ) 
    {
        if( mentrp == 5 )
            CALL( MemProj(psMemsys, 2, 1) )
        CALL( MemVmul(psMemsys, 2, 1, 2) )
    }
    CALL( Opus(psUser, afSt[2], afSt[25]) )
    if( flag ) 
    {
        if( Umemx )
        {
            CALL( MemVmul(psMemsys, 25, 29, 25) )
        } 
        else
        {
            if( mgauss == 1 && acc > ZERO )
                CALL( MemSmul(psMemsys, 25, acc, 25) )
            else
                CALL( MemVmul(psMemsys, 25, 22, 25) )
        }
    }
    CALL( MemCopy(psMemsys, 25, k) )
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemTr
 * 
 * Purpose:  Differential backward transform, with or without full scaling 
 * 
 *           if( flag ) then
 *                        <j> := sqrt([metric]) * Tropus * [acc] * <k> 
 *           else
 *                        <j> :=                  Tropus * [acc] * <k>
 * 
 * Notes:    (1)  Tropus must NOT change frequently between iterates
 *           (2)  The actual Tropus transform is performed from <25> to <2>
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         20 Feb 1991     iEntropy = 5 added
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemTr(
const MEMSYS* psMemsys, 
int     k,       /* I    Input data area number */
int     j,       /* I    Output MaxEnt area number */
short   flag)    /* I    Tropus / SetDistribution switch */
{
/* Areas:
 *      <k>         I                Input data area
 *      <j>           O              Output MaxEnt area
 *      <1>        (I)                 (metric)
 *      <2>          (O)               (transform workspace)
 *     <22>         I                [acc]    (if linear)
 *     <25>          (O)               (transform workspace)
 *     <29>         I                [acc]    (if nonlinear)
 */
    int      mentrp                             = psMemsys->psInput->iEntropy;
    int      mgauss                             = psMemsys->psInput->iGauss;
    float    acc                                = psMemsys->psInput->fAcc;
    void*    psUser                             = psMemsys->psVpoint->psUser;
    float**  afSt                               = psMemsys->psOutput->afSt;
    int    (*Tropus)(void*,const float*,float*) = psMemsys->psProcs->pfnTropus;
    int    (*Umemx) (float,float*,float*)       = psMemsys->psProcs->pfnUmemx;
    int      CALLvalue = 0;

#if FLOWCHART
    printf("   MemTr   (%2d,%2d)\n", k, j);
#endif
    CALL( MemCopy(psMemsys, k, 25) )
    if( Umemx )
    {
        CALL( MemVmul(psMemsys, 25, 29, 25) )
    }
    else 
    {
        if( mgauss == 1 && acc > ZERO )
            CALL( MemSmul(psMemsys, 25, acc, 25) )
        else 
            CALL( MemVmul(psMemsys, 25, 22, 25) )
    } 
    CALL( Tropus(psUser, afSt[25], afSt[2]) )
    if( flag )
    {
        CALL( MemVmul(psMemsys, 2, 1, 2) )
        if( mentrp == 5 )
            CALL( MemProj(psMemsys, 2, 1) )
    }
    CALL( MemCopy(psMemsys, 2, j) )
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemQuant
 * 
 * Purpose:  
 *                ------------------------------------------
 *               |          T                        -1     |
 *           Set |   Q :=  b . ( Alpha*Identity + A )  . b  |
 *               |         -                             -  |
 *                ------------------------------------------
 *      where      A = Memtr.Memop ,
 *                 Memop = ( [acc] *  Opus  * sqrt([metric]) *)
 *                 Memtr = ( sqrt([metric]) * Tropus * [acc] *)
 *  
 *   IF ( bubble ) then additionally
 *                ------------------------------------------
 *               |                               -1/2       |
 *           set |   t :=  ( Alpha*Identity + A )    . b    |
 *               |   -                                 -    |
 *                ------------------------------------------
 *           When b is a random normal vector, the statistics of t are
 *                   T                          -1
 *              < t t > = ( Alpha*Identity + A )
 * 
 * Notes: (1) The status codes are
 *                bcodec   TRUE : Conjugate gradient returned OK
 *                         FALSE: Conjugate gradient not satisfactory
 *            and bcodet   TRUE : Bubble not called or bubble sample OK
 *                         FALSE: Bubble offset has detectably incorrect length
 *
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemQuant( /*   O  # of transforms */
const MEMSYS* psMemsys, 
float   alpha,        /* I    Regularising parm */
float   utol,         /* I    Finish tolerance */
short   bubble,       /* I    Switch for "bubble" option */
float*  qlo,          /*   O  Numerical lower limit on Q */
float*  qhi,          /*   O  Numerical upper limit on Q */
short*  bcodec,       /*   O  Status code for conj. gradients */
short*  bcodet)       /*   O  Status code for bubble sample */
/* Areas:
 *      <1>          I                sqrt([metric])
 *      <2>            O              t   (only if bubble)
 *      <4>          I                b = input vector
 *     <22>         (I)               [acc]    (if linear and acc <= 0)
 *     <25>           (O)                      workspace
 *     <26>           (O)                      workspace
 *     <27>           (O)                      workspace (only if bubble) 
 *     <28>           (O)                      workspace
 *     <29>         (I)               [acc]    (if nonlinear)
 */
{
    static  float gam[LMAX+1];
    static  float del[LMAX+1];
    static  float off[LMAX+1];
    static  float val[LMAX+1];
    static  float vec[(LMAX+1)*(LMAX+1)];
    float   temp;
    int     i;
    int     j;
    int     k;
    int     m;
    int     n;
    int     ntrans;
    int     CALLvalue = 0;

    ntrans = 0;
    CALL( MemCgh(psMemsys, alpha, utol, FALSE,
                 qlo, qhi, bcodec, &n, gam, del, off) )
    ntrans += CALLvalue;
    *bcodet = TRUE;
    if( bubble ) 
    {
/* Calculate coefficients OFF() of gradient vectors g for random sample */
        val[0] = gam[0] * del[0];
        for( i = 1; i <= n; ++i )
        {
            val[i] = gam[i] * del[i];
            off[i] = - gam[i] * del[i - 1];
        }
        m = n + 1;
        CALL( MemDet(val, vec, off, m, m, 1) )
        for( i = 0; i <= n; ++i )
            off[i] = ZERO;
        k = 0;
        for( j = 0; j <= n; ++j ) 
        {
            temp = vec[k] / SQRT(alpha + val[j]);
            for( i = 0; i <= n; ++i )
                off[i] += vec[k++] * temp * gam[0] / gam[i];
        }
/* Build up output vector */
        CALL( MemCgh(psMemsys, alpha, utol, TRUE,
                     qlo, qhi, bcodet, &n, gam, del, off) )
        ntrans += CALLvalue;
    }
    return ntrans;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemCgh
 * 
 * Purpose:  
 *    Object-space conjugate gradient
 *    -------------------------------
 *  
 *     Estimate q = [ qlo , qhi ]  and set up gamma and delta scalars
 *       by maximising
 *                           qb =  2 y.b  -  y.B.y    (scalars only)
 *                          qab = 2 z.A.b - z.A.B.z   (scalars only)
 *  
 *     Slave these to master conjugate gradient maximisation of
 *                           qa =  2 x.b  -  x.A.x    (scalars only)
 *      where      B = alpha*I + A    ,  A = Memtr.Memop ,
 *                 Memop = ( [acc] *  Opus  * sqrt([metric]) *)
 *                 Memtr = ( sqrt([metric]) * Tropus * [acc] *)
 *  
 *     if( flag ) then additionally construct the vector
 *                         N
 *             t(.)  :=   Sum coeff(i) * g(.,i)
 *                        i=0
 *         repeating conjugate gradient to regenerate the gradients g,
 *         (as simulated in dataspace).
 * 
 * Notes: (1) Conjugate gradient termination is either when the user's
 *            criterion is satisfied, by
 *            (1.+utol) * alpha * qb + qab   passing   (length of b)**2
 *            or when loop counter passes
 *            LMAX = maximum number of matrix applications
 *            or when gamma underflows the internal arithmetic tolerance.
 *        (2) A repeat run (with flag = TRUE) to calculate the vector output
 *            must reproduce the first run exactly, otherwise serious rounding
 *            errors can accumulate.  This is why the same code is re-used.
 *        (3) IF( ! flag ) the status code relates to the conjugate gradient:
 *              bcode = ( icode <= 2 )
 *            where the internal codes are
 *              icode = 0: User's termination criterion satisfied      (OK)
 *              icode = 1: Convergence to within arithmetic tolerance  (OK)
 *              icode = 2: Zero results from eps input vector b      (Warning)
 *              icode = 3: Iteration limit LMAX exceeded             (Warning)
 *              icode = 4: Probable error in transform routines       (Error)
 *
 * IF( flag ) the status code relates to the output vector length
 *     bcode = ( t.t is within a fraction utol of qlo )
 * this being a consistency test when called from MemQuant.
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemCgh( /*   O  # of transforms */
const MEMSYS* psMemsys, 
float   alpha,      /* I    Regularising parm */
float   utol,       /* I    Finish tolerance */
short   flag,       /* I    Output vector switch */
float*  qlo,        /*   O  Numerical lower limit on Q */
float*  qhi,        /*   O  Numerical upper limit on Q */
short*  bcode,      /*   O  Status code */
int*    n,          /*   O  gam[i], del[i], for i = 0,...,n */
float*  gam,        /*   O  Gamma from each iteration */
float*  del,        /*   O  Delta from each iteration */
float*  coeff)      /* I    Coefficients of g (only if flag) */
/* Areas:
 *      <1>         I                sqrt([metric])
 *      <2>           O              t = output (if flag), and workspace 
 *      <4>         I                b = input vector
 *     <22>        (I)               [acc]    (if linear)
 *     <25>          (O)             uu (workspace for transforms)
 *     <26>          (O)             gg = simulated gradient  (workspace) 
 *     <27>          (O)             tt = simulated output t  (if flag)
 *     <28>          (O)             hh = simulated conjugate (workspace) 
 *     <29>        (I)               [acc]    (if nonlinear)
 */
{
    int     icode;
    int     l;
    int     ntrans;
    float   delb;
    float   phib;
    float   epsb;
    float   temp;
    float   delab;
    float   gamma0;
    float   t0;
    float   qb;
    float   gamma;
    float   phiab;
    float   scale;
    float   delta;
    float   epsab;
    float   qab;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemCgh\n");
#endif
/* Initialise */
    l = ntrans = 0;
    qb = qab = t0 = ZERO;
    qab = ZERO;
    *n = 0;
    gam[0] = ZERO;
    del[0] = ZERO;
/* Scale arbitrarily (make input vector unit) to reduce overflow risk */
    CALL( MemDot(psMemsys, 4, 4, &temp) )
    if( temp < TINY ) 
    {
        gamma0 = ZERO;
        t0 = ZERO;
        scale = ONE;
        icode = 2;
        goto Exit;
    }
    scale = SQRT(temp);
    gamma0 = (temp / scale) / scale;
    gamma = gamma0;
    gam[l] = SQRT(gamma) * scale;
    temp = ONE / scale;
    CALL( MemSmul(psMemsys, 4, temp, 2) )
    CALL( MemZero(psMemsys, 26) )
    CALL( MemZero(psMemsys, 28) )
    if( flag )
    {
        CALL( MemZero(psMemsys, 27) )
        t0 = coeff[0];
    }
    epsb = ONE;
    phib = alpha / gamma;
    epsab = gamma;
    phiab = ZERO;
    delab = ZERO;
/* Loop */
    for( ; ; )
    {
        temp = ONE / gamma;
        CALL( MemOp (psMemsys, 2, 25, TRUE) )
        ++ntrans;
        CALL( MemSmd(psMemsys, 25, temp, 28, &delta) )
        *n = l;
        del[l++] = SQRT(delta) / scale;
        delb = phib + delta / epsb / epsb;
        qb += ONE / delb;
        if( alpha * qb * utol + alpha * qb + qab >= gamma0 )
        {
            icode = 0;
            break;
        }
        if( delta <= TINY )
        {
            icode = 4;
            break;
        }
        if( l > LMAX )
        {
            icode = 3;
            break;
        }
        temp = ONE / epsab - delab * epsab / gamma;
        phiab += alpha / (epsab * delta * epsab) + temp * gamma * temp;
        temp = - ONE / delta;
        CALL( MemSma(psMemsys, 28, temp, 26, 26) )
        CALL( MemTr (psMemsys, 26, 2, TRUE) )
        ++ntrans;
        temp = ONE / scale;
        CALL( MemSmd(psMemsys, 4, temp, 2, &gamma) )
        gam[l] = SQRT(gamma) * scale;
        delab = phiab + gamma / epsab / epsab;
        epsab = gamma / (epsab * delab);
        qab += ONE / delab;
        if( alpha * qb * utol + alpha * qb + qab >= gamma0 ) 
        {
            icode = 0;
            break;
        }
        if( gamma <= gamma0 * FLT_EPSILON )
        {
            icode = 1;
            break;
        }
        temp = ONE / epsb;
        epsb = delta / (epsb * delb);
        phib += alpha / (epsb * gamma * epsb) +
               (ONE / epsb - temp) * delta * (ONE / epsb - temp);
        if( flag )
        {
            CALL( MemSma(psMemsys, 26, coeff[l], 27, 27) )
            t0 += coeff[l];
        }
    }
/* Rescale and exit */
Exit:
    *qlo = scale * qb * scale;
    *qhi = scale * (gamma0 - qab) * scale / alpha;
    *bcode = icode <= 2;
    if( flag )
    {
        CALL( MemTr (psMemsys, 27, 2, TRUE) )

        ++ntrans;
        temp = t0 / scale;
        CALL( MemSmd(psMemsys, 4, temp, 2, &temp) )
        CALL( MemSmul(psMemsys, 2, scale, 2) )
        temp = scale * temp * scale - *qlo;
        temp = ABS(temp);
        *bcode = (temp <= utol * *qlo);
    }
    return ntrans;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemCgd
 * Purpose:  
 * 
 *   Data-space conjugate gradient
 *   -----------------------------
 *                    ------------------------------------
 *     Find          |                      -1            |
 *                   |    z = ( beta*I + A )  b           |
 *     by maximising |                                    |
 *                   |  qab = 2 z.A.b - z.A.(beta*I+A).z  |
 *                    ------------------------------------
 *             where  beta >= alpha  ensures (if flag is TRUE) that
 *                     dist = dimensionless length of z
 *                          = sqrt( z.A.z / summet )      <=   rate
 *     Return  dist  and  h = [ hlow , hhigh ] = qab/2 = cross-entropy
 *  
 *                        ---------------------------
 *     Also set up       |  gamma and delta scalars  |    defining A
 *                        ---------------------------
 *  
 *             where      A = Memop.Memtr
 *                    Memtr = ( sqrt([metric]) * Tropus * [acc] *)
 *                    Memop = ( [acc] *  Opus  * sqrt([metric]) *)
 *  
 *      In terms of normalised (and orthogonal) gradient vectors
 *        --------       --------  -------- 
 *       |        |     |        ||\       |
 *       |   ^    |     |        ||  1/gam |
 *       |   g    |  =  |   g    ||       \|
 *       |        |     |        | -------- 
 *       |        |     |        |
 *        --------       -------- 
 *      the matrix A becomes represented by the factorised tridiagonal form
 *        --------       --------  --------  --------  --------  --------
 *       | ^t   ^ |     |\       || 1      ||\       || 1 -1   ||\       |
 *       | g .A.g |  =  |  gam   ||-1  1   || del**2 ||    1 -1||  gam   |
 *       |        |     |       \||   -1  1||       \||       1||       \|
 *        --------       --------  --------  --------  --------  --------
 * 
 * Notes: (1) Maximise qab and  qb = 2 y.b - y.(beta*I+A).y  and set gam, del
 *            by slaving them to master conjugate gradient maximisation of
 *                           qa = 2 x.b - x.A.x    (scalars only)
 *        (2) Termination is either when qab has tolerably converged, by
 *            (1.+utol) * qab + beta * qb   passing   (length of b)**2
 *            or when loop counter passes
 *            LMAX = maximum number of matrix applications
 *            or when gamma underflows the internal arithmetic tolerance.
 *        (3) On exit, either m = n or m = n+1 .
 *        (4) Within MemSys, the subspace approximation (g,gam,del,m,n) to
 *            the matrix A is not actually needed if flag is TRUE
 *            Likewise, the vector z and its scalars are not actually needed
 *            if flag is FALSE
 *        (5) The status codes are
 *               bcodeb = ( beta == alpha , no distance penalty being invoked)
 *             and
 *               bcodec = ( icode <= 2 )
 *            where the internal codes are
 *               icode = 0: User's termination criterion satisfied      (OK)
 *               icode = 1: Convergence to within arithmetic tolerance  (OK)
 *               icode = 2: Zero results from eps input vector b      (Warning)
 *               icode = 3: Iteration limit LMAX exceeded             (Warning)
 *               icode = 4: Probable error in transform routines       (Error)
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemCgd( /*   O  # of transforms */
const MEMSYS* psMemsys, 
float   rate,       /* I    Dimensionless distance limit */
float   summet,     /* I    SUM([metric]) */
float   alpha,      /* I    Regularising parm */
float   utol,       /* I    Finish tolerance */
short   flag,       /* I    Update <23> with z ? */
float*  hlow,       /*   O  Lower limit on cross-entropy */
float*  hhigh,      /*   O  Upper limit on cross-entropy */
float*  dist,       /*   O  Dimensionless distance moved */
float*  gam,        /*   O  Gamma from each iteration ( ! flag ) */
float*  del,        /*   O  Delta from each iteration ( ! flag ) */
int*    m,          /*   O  gam(0..m) set on exit     ( ! flag ) */
int*    n,          /*   O  del(0..n) set on exit     ( ! flag ) */
short*  bcodeb,     /*   O  Status code for beta control */
short*  bcodec)     /*   O  Status code for conjugate gradients */
/* Areas:
 *      <1>          I                sqrt([metric])
 *      <2>           (O)             u (transform workspace)
 *     <22>         (I)               [acc]    (if linear)
 *     <23>          I O              <23> += z                 (if flag)
 *     <25>           (O)             u (transform workspace)
 *     <26>          I(O)             b input vector (& g gradient for x)
 *     <27>           (O)             w conjugate vector for z  (if flag)
 *     <28>           (O)             h conjugate vector for x
 *     <29>         (I)               [acc]    (if nonlinear)
 */
{
#if CONTROL
    FILE*   LogFile = psMemsys->fpLog;
    int     (*Uinfo)(FILE*,const char*) = psMemsys->psProcs->pfnUinfo;
    char    info[80];
#endif
    int     icode;
    int     l;
    int     ntrans;
    float   delb;
    float   beta;
    float   phib;
    float   epsb;
    float   temp;
    float   delab;
    float   d2;
    float   gamma0;
    float   qb;
    float   qab;
    float   gamma;
    float   phiab;
    float   scale;
    float   delta;
    float   epsab;
    float   zaz;
    float   zx;
    float   zy;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemCgd\n");
#endif
/* Initialise */
    l = ntrans = 0;
    if( m )
        *m = 0;
    if( n )
        *n = 0;
    epsb = ONE;
    qb = ZERO;
    if( flag )
        CALL( MemZero(psMemsys, 27) )
    qab = ZERO;
/* Scale arbitrarily (make input vector unit) to reduce overflow risk
 * Internal variables rescale proportionally to powers of scale:
 * u +1,               gamma -2, delta +2,     g -1,   h +1
 *              qb -2,   epsb 0,  phib +2,  delb +2
 * w +1, z -1, qab -2,  epsab-2, phiab +2, delab +2, zaz -2, zx +2, zy 0
 */
    CALL( MemDot(psMemsys, 26, 26, &temp) )
    if( temp < TINY )
    {
        icode = 2;
        scale = ONE;
        zaz = gamma0 = ZERO;
        if( gam )
            gam[0] = ZERO;
        if( del )
            del[0] = ZERO;
        goto Exit;
    }
    scale = SQRT(temp);
    gamma0 = temp / scale / scale;
    gamma = gamma0;
    if( gam )
        gam[0] = SQRT(gamma) * scale;
    temp = ONE / scale;
    CALL( MemSmul(psMemsys, 26, temp, 26) )
    epsab = gamma;
    phiab = ZERO;
    delab = ZERO;
    temp = ONE / gamma;
    CALL( MemSmul(psMemsys, 26, temp, 28) )
    CALL( MemTr  (psMemsys, 28, 2, TRUE) )
    ++ntrans;
    CALL( MemDot (psMemsys, 2, 2, &delta) )
    CALL( MemOp  (psMemsys, 2, 25, TRUE) )
    ++ntrans;
    if( del )
        del[0] = SQRT(delta) / scale;
    beta = alpha;
    if( bcodeb )
        *bcodeb = TRUE;
    if( flag )
    {
/*  Beta control: D2 = h.A.A.h */
        CALL( MemDot(psMemsys, 25, 25, &d2) )
        temp = scale * scale / (rate * summet * rate);
        CALL( MemLb(alpha, rate, temp, delta, d2, &beta, bcodeb) )
#if CONTROL
        SPRINTF(info,"      beta  %16.7e\n", beta);
        CALL( Uinfo(LogFile, info) )
#endif
    }
    zx = zy = zaz = ZERO;
    phib = beta / gamma;
/* Loop */
    for( ; ; )
    {
        ++l;
        delb = phib + delta / epsb / epsb;
        qb += ONE / delb;
        if( utol * qab + qab + beta * qb >= gamma0 ) 
        {
            icode = 0;
            break;
        }
        if( delta <= TINY ) 
        {
            icode = 4;
            break;
        }
        if( l > LMAX )
        {
            icode = 3;
            break;
        }
        if( l > 1 )
        {
            CALL( MemOp(psMemsys, 2, 25, TRUE) )
            ++ntrans;
        }
        temp = ONE / epsab - delab * epsab / gamma;
        phiab += beta / (epsab * delta * epsab) + temp * gamma * temp;
        temp = - ONE / delta;
        CALL( MemSmd(psMemsys, 25, temp, 26, &gamma) )
        delab = phiab + gamma / epsab / epsab;
        if( flag )
        {
            temp = ONE / (epsab * delta);
            CALL( MemSma(psMemsys, 28, temp, 27, 27) )
            temp = scale / delab;
            CALL( MemSma(psMemsys, 27, temp, 23, 23) )
        }
        else if( m )
        {
            *m = l;
            gam[*m] = SQRT(gamma) * scale;
        }
        zaz += zy / delab;
        zx  += ONE / (epsab * delta * epsab);
        zy  += zx / delab;
        zaz += zy / delab;
        epsab = gamma / (epsab * delab);
        qab += ONE / delab;
        if( utol * qab + qab + beta * qb >= gamma0 )
        {
            icode = 0;
            break;
        }
        if( gamma <= gamma0 * FLT_EPSILON ) 
        {
            icode = 1;
            break;
        }
        temp = ONE / gamma;
        CALL( MemSma(psMemsys, 26, temp, 28, 28) )
        temp = ONE / epsb;
        epsb = delta / (epsb * delb);
        phib += beta / (epsb * gamma * epsb) +
               (ONE / epsb - temp) * delta * (ONE / epsb - temp);
        CALL( MemTr (psMemsys, 28, 2, TRUE) )

        ++ntrans;
        CALL( MemDot(psMemsys, 2, 2, &delta) )
        if( !flag ) 
            if( n )
            {
                *n = l;
                del[*n] = SQRT(delta) / scale;
            }
    }
/* Rescale and exit */
Exit:
    if( bcodec )
        *bcodec = (icode <= 2);
    if( hlow )
        *hlow = scale * qab * scale / TWO;
    if( hhigh )
        *hhigh = scale * (gamma0 - alpha * qb) * scale / TWO;
    if( dist )
        *dist = scale * SQRT(zaz / summet);
    return ntrans;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemDet
 * 
 * Purpose:  Eigenstructure of symmetric tridiagonal matrix
 *                   t            t
 *            A = L L   or   A = L  L    where L is lower bidiagonal.
 * 
 * Notes: (1)  val must be strictly positive, except that val[n-1] may be 0.
 *             off must be strictly negative, except that off[n-1] may be 0.
 *             The internal signs in this routine are specifically chosen for
 *             this case, to ensure that the rotations preserve this form.
 *        (2)  The final run, after the input angle has been found by binary
 *             chop, must repeat the final simulation exactly, which is why
 *             the same code is re-used.
 *        (3)  DO NOT try to edit this - almost all changes will cause at
 *             least occasional damage!
 *        (4)  The eigenvalues return in decreasing order, unless an early
 *             one is tolerably small - see loop at instruction "Small".
 *        (5)  Running time proportional to  n * n * ( m + log(2) / FLT_EPS) )
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         17 Feb 1991     n <= 1 return fixed
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemDet(
float*  val,   /* I    Diagonal elements 0..n-1 of L
                    O  Eigenvalues of A                                   */
float*  vec,   /*   O  vec[k,.] is 1st part of evec k of A (m components) */
float*  off,   /* I    Subdiagonal elements 1..n-1 of L (off[0] unused)
                    O               (workspace)                           */
int     n,     /* I    Dimension                                          */
int     m,     /* I    Number of evector components needed                */
int     iflag) /* I    +ve is A = L Lt , -ve is A = Lt L                  */
{
    int     i;
    int     j;
    int     k;
    int     l;
    float   dcold;
    float   dsold;
    float   dc;
    float   ds;
    float   chi;
    float   ceq;
    float   clo;
    float   shi;
    float   seq;
    float   slo;
    short   exit;
    short   finish;
/* Register...... */
    register float c;
    register float s;
    register float r;
    register float x;
    register float d0;
    register float d1;
    register float e0;
    register float e1;
/* .....variables */

    off[0] = ZERO;
/* Set vectors vec to be rotated into eigenvector coords. */
    for( j = 0; j < m; ++j )
    {
        for( i = 0; i < n; ++i )
            vec[j + i * m] = ZERO;
        vec[j + j * m] = ONE;
    }
/* Exceptional return if n == 1 */
    if( n > 1 )
    {
/* Iterate each successive off[k] ( k = n-1,n-2,...) towards 0. */
/* Ignore any early off[l] ( l = 0,1,2,...) which are tolerably small. */
        k = n - 1;
        l = 0;
Small:
        if( off[l + 1] + FLT_EPSILON * (val[l] + val[l + 1]) > ZERO )
        {
            ++l;
            if( l < k - 1 )
                goto Small;
        }
/* Initialise binary chop, using both Cos and Sin for full accuracy. */
        do
        {
            finish = FALSE;
            clo = ONE;
            slo = ZERO;
            ceq = HALF;
            seq = - HALF;
            chi = ZERO;
            shi = - ONE;
            dc = ONE;
            ds = ONE;
/* Chop the initial rotation angle (Q) between 0 and -pi/2, */
/* until the last off-diagonal element vanishes. */
            do
            {
                exit = finish;
                r = SQRT(ceq * ceq + seq * seq);
                c = ceq / r;
                s = seq / r;
/*     Enter with lower bidiagonal matrix L and rotation Q (C=cos,S=sin). 
 *     Apply Q to L,
 *              -----------   -----------
 *             | C -S      | |d0         |
 *             | S  C      | |e1 d1      |   (e.g. dimension n=4)
 *      Q L =  |       1   | |   e2 d2   |   (d in VAL, e in OFF)
 *             |          1| |      e3 d3|
 *              -----------   -----------
 *     with such other right and left rotations as are needed to preserve 
 *     the lower bidiagonal form of A.  Successively (left to right):-
 *  -------   -------   -------   -------   -------   -------   -------
 * |d      | |d X    | |d      | |d      | |d      | |d      | |d      | 
 * |e d    | |e d    | |e d    | |e d X  | |e d    | |e d    | |e d    | 
 * |  e d  | |  e d  | |X e d  | |  e d  | |  e d  | |  e d X| |  e d  | 
 * |    e d| |    e d| |    e d| |    e d| |  X e d| |    e d| |    e d| 
 *  -------   -------   -------   -------   -------   -------   -------
 */
                i = l;
                d1 = val[l];
                e1 = off[l + 1];
Angle:
                d0 = d1;
                d1 = val[++i];
                r = c * e1 - s * d0;
                d0 = s * e1 + c * d0;
                e1 = r;
                x = s * d1;
                d1 = c * d1;
                if( exit )
                {
                    if( iflag >= 0 )
                    {
                        for( j = 0; j < m; ++j )
                        {
                            r = c * vec[j + i * m] -
                                s * vec[j + (i-1) * m];
                            vec[j + (i-1) * m] = s * vec[j + i * m] +
                                                 c * vec[j + (i-1) * m];
                            vec[j + i * m] = r;
                        }
                    }
                } 
/*     X is now overflow element on upper-right
 *                ----
 *     Apply  i-1|C -S|    on right to eliminate X
 *             i |S  C|
 *                ----
 */
                r = ABS(d0) + ABS(x);
                c = d0 / r;
                s = x / r;
                r *= SQRT(c * c + s * s);
                c = d0 / r;
                s = x / r;
                d0 = s * x + c * d0;
                r = c * d1 - s * e1;
                e1 = s * d1 + c * e1;
                d1 = r;
                if( exit )
                {
                    val[i - 1] = d0;
                    if( iflag < 0 )
                    {
                        for( j = 0; j < m; ++j )
                        {
                            r = c * vec[j + i * m] -
                                s * vec[j + (i - 1) * m];
                            vec[j + (i-1) * m] = s * vec[j + i * m] +
                                                 c * vec[j + (i-1) * m];
                            vec[j + i * m] = r;
                        }
                    }
                }
/* The zeros of the last off-diagonal terms off[i] at successive
 * dimensions i interleave.  Hence, if any off[i] would change sign,
 * then reduce the initial angle to approach the first zero of off[k].
 */
                if( e1 < ZERO )
                {
                    if( i < k )
                    {
                        e0 = e1;
                        e1 = off[i + 1];
                        x = s * e1;
                        e1 = c * e1;
/*     x is now overflow element on lower-left
 *                ----
 *     Apply   i | C S|    on left to eliminate x
 *            i+1|-S C|
 *                ----
 */
                        r = ABS(e0) + ABS(x);
                        c = e0 / r;
                        s = x / r;
                        r = - r * SQRT(c * c + s * s);
                        c = e0 / r;
                        s = x / r;
                        if( exit )
                            off[i] = r;
                        goto Angle;
                    }
                    clo = ceq;
                    slo = seq;
/* Exit if last off-diagonal element off[k] is tolerably small. */
                    finish = finish || e1 + FLT_EPSILON * (d0 + d1) > ZERO;
                } 
                else 
                {
                    chi = ceq;
                    shi = seq;
                }
                dcold = dc;
                dsold = ds;
                dc = clo - chi;
                ds = slo - shi;
                finish = finish || (dcold <= dc && dsold <= ds);
/* Repeat the final successful (low) iterate to get output */
                if( finish )
                {
                    ceq = clo;
                    seq = slo;
                } 
                else 
                {
                    ceq = chi + dc / TWO;
                    seq = slo - ds / TWO;
                }
            } while( ! exit );
            val[i] = d1;
/* Exit having re-used the rotation angle which kills the last        */
/* off-diagonal term off[k], effectively reducing the dimension by 1. */
/* e1 contains new off[k], which should be small negative diagnostic. */
            --k;
        } while( k > l );
    }
/* Square up eigenvalues on exit */
    for( i = 0; i < n; ++i )
        val[i] = val[i] * val[i];
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemLa0
 * 
 * Purpose:  Alpha control if alpha = infinity
 * 
 * Notes: 
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemLa0(
const MEMSYS* psMemsys, 
float   rate,   /* I    Dimensionless distance limit */
float   summet, /* I    SUM([metric]) */
float   agrads, /* I    alpha*gradS */
float   omega,  /* I    value of stopping criterion */
float   var,    /* I    Intrinsic variance of (omega) */
float*  alpha,  /*   O  Output alpha */
short*  bcode)  /*   O  Alpha control diagnostic */
{
#if CONTROL
    FILE*   LogFile = psMemsys->fpLog;
    int     (*Uinfo)(FILE*,const char*) = psMemsys->psProcs->pfnUinfo;
    char    info[80];
#endif
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemLa0\n");
#endif
    CALL( MemLaw(psMemsys, BIG, omega, var, 1) )
    if( omega < ONE ) 
    {
        *alpha = agrads / SQRT(rate * summet * rate);
        *bcode = FALSE;
    } 
    else 
    {
        *alpha = BIG;
        *bcode = TRUE;
    }
#if CONTROL
    SPRINTF(info,"     alpha  %16.7e\n", *alpha);
    CALL( Uinfo(LogFile, info) )
#endif
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemLa
 * 
 * Purpose:  Alpha control, iterating towards Omega rising to 1.00
 * 
 * Notes: 
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemLa(
const MEMSYS* psMemsys, 
float*  alpha,   /* I O  Regularising parm */
float   rate,    /* I    Dimensionless distance limit */
float   summet,  /* I    SUM([metric]) */
float   omega,   /* I    value of stopping criterion */
float   var,     /* I    Intrinsic variance of (omega) */
float   agrads,  /* I    alpha*gradS */
short   init,    /* I    Flag for initialisation of table */
short   bcodeb,  /* I    Beta control diagnostic */
short*  bcodea,  /*   O  Alpha control diagnostic */
short   bcodet)  /* I    Test diagnostic */
{
#if CONTROL
    FILE*   LogFile = psMemsys->fpLog;
    int     (*Uinfo)(FILE*,const char*) = psMemsys->psProcs->pfnUinfo;
    char    info[80];
#endif
    float   alf1;
    float   alf2;
    float   y1;
    float   y2;
    float   ynew;
    float   r;
    float   sigma;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   Control\n");
#endif
    CALL( MemLaw(psMemsys, *alpha, omega, var, init) )
#if CONTROL
    SPRINTF(info,"     rate   %16.7e\n", rate);
    CALL( Uinfo(LogFile, info) )
    SPRINTF(info,"     summet %16.7e\n", summet);
    CALL( Uinfo(LogFile, info) )
    SPRINTF(info,"     al*grS %16.7e\n", agrads);
    CALL( Uinfo(LogFile, info) )
    SPRINTF(info,"     alpha  %16.7e\n", *alpha);
    CALL( Uinfo(LogFile, info) )
#endif
    if( bcodeb && bcodet ) 
    {
        CALL( MemLar(psMemsys, *alpha, &ynew, &sigma) )
        if( ABS(ynew) < sigma)
            *bcodea = TRUE;
        else 
        {
            *bcodea = FALSE;
            r = *alpha * SQRT(rate * summet * rate) / agrads + ONE;
            alf1 = *alpha * r;
            alf2 = *alpha / r;
            alf1 = MAX(alf1, *alpha);
            alf2 = MIN(alf2, *alpha);
            CALL( MemLar(psMemsys, alf1, &y1, &sigma) )
            if( y1 > ZERO )
            {
                *alpha = alf1;
            }
            else 
            {
                CALL( MemLar(psMemsys, alf2, &y2, &sigma) )
                if( y2 <= ZERO ) 
                {
                    *alpha = alf2;
                } 
                else 
                {
                    r = ONE;
                    do
                    {
                        r /= TWO;
                        *alpha = (alf1 + alf2) / TWO;
                        CALL( MemLar(psMemsys, *alpha, &ynew, &sigma) )
                        if( ynew <= ZERO )
                        {
                            alf1 = *alpha;
                            y1 = ynew;
                        }
                        else
                        {
                            alf2 = *alpha;
                            y2 = ynew;
                        }
                    } while( r > FLT_EPSILON );
                    if( y1 - y2 != ZERO )
                    {
                          r = y1 / (y1 - y2);
                          *alpha = alf1 + (alf2 - alf1) * r;
                    }
                }
            }
        }
    }
    else 
    {
        *bcodea = TRUE;
    }
#if CONTROL
    SPRINTF(info,"     alpha  %16.7e\n", *alpha);
    CALL( Uinfo(LogFile, info) )
#endif
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemLaw
 * 
 * Purpose:  Insert new entry into table of omega(alpha) with variances,
 *           discarding least significant earlier entry if table is full.
 * 
 * Notes: (1) Deletion algorithm is ad hoc.
 * 
 * History:
 *      JS              7 Dec 1988     MEMSYS3 version
 *      MKC/JS         19 Dec 1990     First release
 *      JS              8 Mar 1991     Overflow risk reduced
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemLaw(
const MEMSYS* psMemsys, 
float   alpha,  /* I    New alpha */
float   omega,  /* I    New omega(alpha) */
float   var,    /* I    Intrinsic variance of log(omega) */
short   init)   /* I    Initialisation flag */
{
    int*    ntable = &(psMemsys->psOutput->iNtable);
    float*  xtable =   psMemsys->psOutput->fXtable;
    float*  ytable =   psMemsys->psOutput->fYtable;
    float*  vtable =   psMemsys->psOutput->fVtable;
    float   d;
    float   s;
    float   x;
    float   xnew;
    float   ynew;
    float   vnew;
    int     i;
    int     j;
#if CONTROL
    FILE*   LogFile = psMemsys->fpLog;
    int     (*Uinfo)(FILE*,const char*) = psMemsys->psProcs->pfnUinfo;
    char    info[80];
    int     CALLvalue = 0;
#endif

    xnew = LOG(alpha);
    ynew = LOG((MAX(omega, EPS)));
    vnew = MAX(var, EPS);
    if( init )
        *ntable = 0;
    for( i = 0; i < *ntable; ++i )
    {
        if( xnew == xtable[i] ) 
        {
            s = vnew / (vnew + vtable[i]);
            x = vtable[i] / (vnew + vtable[i]);
            ytable[i] = ynew * x + ytable[i] * s;
            vtable[i] = vnew * x;
            return 0;
        }
    }
    if( *ntable >= NSIZE )
    {
/* Delete worst alpha, preserving order */
        s = MAX(var, EPS);
        j = -1;
        for( i = 0; i < *ntable; ++i ) 
        {
            x = LOG(alpha) - xtable[i];
            d = vtable[i] + x * x * x * x;
            if( d > s )
            {
                s = d;
                j = i;
            }
        }
        if( j == -1 )
            return 0;
        -- *ntable;
        for( i = j; i < *ntable; ++i ) 
        {
            xtable[i] = xtable[i + 1];
            ytable[i] = ytable[i + 1];
            vtable[i] = vtable[i + 1];
        }
    }
/* Insert alpha */
    j = 0;
    for( i = 0; i < *ntable; ++i )
        if( xnew > xtable[i] ) 
            j = i + 1;
    for( i = *ntable - 1; i >= j; --i ) 
    {
        xtable[i + 1] = xtable[i];
        ytable[i + 1] = ytable[i];
        vtable[i + 1] = vtable[i];
    }
    xtable[j] = xnew;
    ytable[j] = ynew;
    vtable[j] = vnew;
    ++ *ntable;
#if CONTROL
    SPRINTF(info,"   Control on  log alpha       log omega       var omega\n");
    CALL( Uinfo(LogFile, info) )
    for( i = 0; i < *ntable; ++i )
    {
        SPRINTF(info," %4d %19.7e %15.7e %15.7e\n",
                      i, xtable[i], ytable[i], vtable[i]);
        CALL( Uinfo(LogFile, info) )
    }
#endif
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemLar
 * 
 * Purpose:  Interpolate table to estimate log(omega(alpha)) and its error
 * 
 * Notes: 
 * 
 * History:
 *      JS              7 Dec 1988     MEMSYS3 version
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemLar(
const MEMSYS* psMemsys, 
float   alf,    /* I    Trial alpha */
float*  yval,   /*   O  Estimate of log(omega(alpha)) */
float*  sigma)  /*   O  Standard deviation of estimate */
{
    int     ntable = psMemsys->psOutput->iNtable;
    float*  xtable = psMemsys->psOutput->fXtable;
    float*  ytable = psMemsys->psOutput->fYtable;
    float*  vtable = psMemsys->psOutput->fVtable;
    float   xbar;
    float   ybar;
    float   x;
    float   y;
    float   w;
    float   wx;
    float   wy;
    float   wxx;
    float   wxy;
    float   wt;
    int     i;

    w = wx = wxx = wy = wxy = ZERO;
    if( ntable == 1 ) 
    {
        *yval = ytable[0];
        *sigma = ZERO;
    } 
    else 
    {
        for( i = 0; i < ntable; ++i ) 
        {
            x  = LOG(alf) - xtable[i];
            y  = ytable[i];
            wt = vtable[i] + x * x * x * x;
            wt = ONE / MAX(wt, FLT_EPSILON);
            w  += wt;
            wx += wt * x;
            wy += wt * y;
        }
        xbar = wx / w;
        ybar = wy / w;
        for( i = 0; i < ntable; ++i )
        {
            x = LOG(alf) - xtable[i];
            y = ytable[i];
            wt = vtable[i] + x * x * x * x;
            wt = ONE / MAX(wt, FLT_EPSILON);
            x -= xbar;
            y -= ybar;
            wxx += wt * x * x;
            wxy += wt * x * y;
        }
        *yval = ybar - xbar * (wxy / MAX(wxx, FLT_EPSILON));
        *sigma = ONE / SQRT(w);
    }
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemLb
 * 
 * Purpose:  Beta control.
 *  
 *    Given g.g , g.A.g , g.A.A.g , optionally choose beta >= alpha so that 
 *                -1    -1
 *             g.B  .A.B  .g  <  summet * rate * rate   where  B = beta*I + A
 * 
 * Notes: 
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS              9 Oct 1991     Negative rate switch
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemLb(
float   alpha,  /* I    Regularising parm */
float   rate,   /* I    Dimensionless distance limit, O(1) */
float   d0,     /* I    g.g / (summet * rate * rate) */
float   d1,     /* I    g A g / g.g */
float   d2,     /* I    g A A g / g.g */
float*  beta,   /*   O  Revised alpha */
short*  bcode)  /*   O  Return code (FALSE = distance penalty invoked) */
{
    float   a0;
    float   a1;
    float   a2;
    float   a3;
    float   bmin;
    float   q;
    float   r;
    float   scale;
    float   x;
    float   discr;
    float   theta;

#if FLOWCHART
    printf("   MemLb\n");
#endif
    scale = SQRT(d2);
    scale = ONE / MAX(d0,scale);
    a0 = d0 * scale;
    a1 = d1 * scale;
    a2 = scale * d2 * scale;
    if( a1 * FOUR <= a0 ) 
        bmin = SQRT(a0 * a1) - a1;
    else
    {
        a3 = a0 * (a1 * a1 - a2) / FOUR;
        a2 -= a0 * a1;
        a1 *= TWO;
        q = (a1 * a1 - a2 * THREE) / NINE;
        r = (a1 * a1 * a1 * TWO - a1 * NINE * a2 + a3 * TWENTY7) / FIFTY4;
        x = q * q * q - r * r;
        if( x >= ZERO )
        {
            theta = SQRT(x);
            theta = MAX(theta, EPS);
            theta = ATAN(r / theta);
            bmin = (theta + TWO * ATAN(ONE)) / THREE;
            bmin = SQRT((MAX(q, ZERO))) * TWO * COS(bmin) - a1 / THREE;
        }
        else
        {
            discr = EXP(LOG(SQRT(-x) + ABS(r)) / THREE);
            discr = ( r >= ZERO )
                   ?   ABS( discr + q / discr )
                   : - ABS( discr + q / discr );
            bmin = - discr - a1 / THREE;
        }
    }
    bmin /= scale;
    if( bmin > alpha && rate > ZERO )
    {
        *beta = bmin;
        *bcode = FALSE;
    }
    else 
    {
        *beta = alpha;
        *bcode = TRUE;
    }
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemZero
 * 
 * Purpose:              <m> := zero
 *-----------------------------------------------------------------------------
 */
static int  MemZero(
const MEMSYS* psMemsys, 
int     m)       /* I    Output area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemZero (%2d)\n", m);
#endif
    VecFill(st[m], ZERO, len[m]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemSmul
 * 
 * Purpose:      <n> := a * <m>
 *-----------------------------------------------------------------------------
 */
static int  MemSmul(
const MEMSYS* psMemsys, 
int     m,       /* I    Input area number */
float   a,       /* I    Scalar */
int     n)       /* I    Output area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemSmul (%2d,%2d)\n", m, n);
#endif
    VecSmul(st[m], a, st[n], len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemVmul
 * 
 * Purpose:      <n> := <l> * <m>
 *-----------------------------------------------------------------------------
 */
static int  MemVmul(
const MEMSYS* psMemsys, 
int     l,       /* I    Input area number */
int     m,       /* I    Input area number */
int     n)       /* I    Output area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemVmul (%2d,%2d,%2d)\n", l, m, n);
#endif
    VecMul(st[l], st[m], st[n], len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemSma
 * 
 * Purpose:        <n> :=  <l> * x + <m>
 *-----------------------------------------------------------------------------
 */
static int  MemSma(
const MEMSYS* psMemsys, 
int     l,       /* I    Input area number */
float   x,       /* I    Scalar */
int     m,       /* I    Input area number */
int     n)       /* I    Output area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemSma  (%2d,%2d,%2d)\n", l, m, n);
#endif
    VecSmula(st[l], x, st[m], st[n], len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemSmd
 * 
 * Purpose:      <m> :=  <l> * x + <m>
 *              prod :=  <m>.<m>
 *-----------------------------------------------------------------------------
 */
static int  MemSmd(
const MEMSYS* psMemsys, 
int     l,       /* I    Input area number */
float   x,       /* I    Scalar */
int     m,       /* I    Input/Output area number */
float*  prod)    /*   O  Dot product */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemSmd  (%2d,%2d)\n", l, m);
#endif
    VecSmula(st[l], x, st[m], st[m], len[m]);
    VecDot  (st[m], st[m], prod, len[m]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemDot
 * 
 * Purpose:       prod :=  <m>.<n>
 *-----------------------------------------------------------------------------
 */
static int  MemDot(
const MEMSYS* psMemsys, 
int     m,       /* I    Input area number */
int     n,       /* I    Input area number */
float*  prod)    /*   O  Scalar product */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemDot  (%2d,%2d)\n", m, n);
#endif
    VecDot(st[m], st[n], prod, len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemCopy
 * 
 * Purpose:  Copy area m to n   ,   <n> := <m>
 *-----------------------------------------------------------------------------
 */
static int  MemCopy(
const MEMSYS* psMemsys, 
int     m,       /* I    Input area number */
int     n)       /* I    Output area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

#if FLOWCHART
    printf("   MemCopy (%2d,%2d)\n", m, n);
#endif
    if( m != n ) 
        VecMov(st[m], st[n], len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemProj
 * 
 * Purpose:  Project area <m> orthogonally to <n>
 *-----------------------------------------------------------------------------
 */
static int  MemProj(
const MEMSYS* psMemsys, 
int     m,       /* I O    Updated area number */
int     n)       /* I      Orthogonalising area number */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   a;
    float   b;

#if FLOWCHART
    printf("   MemProj (%2d,%2d)\n", m, n);
#endif
    VecDot(st[n], st[n], &a, len[n]);
    VecDot(st[m], st[n], &b, len[n]);
    if( a > ZERO )
        VecSmula(st[n], -b / a, st[m], st[m], len[n]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemTest
 * 
 * Purpose:  test = 1 - cos(angle(<23>,<24>))
 *-----------------------------------------------------------------------------
 */
static int  MemTest(
const MEMSYS* psMemsys, 
float*    test)    /*   O   1 - cos(angle( w, D-Transform(h) )) */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   dd;
    float   wd;
    float   ww;

#if FLOWCHART
    printf("   MemTest\n");
#endif
    VecDot  (st[23], st[23], &ww, len[21]);
    VecDot  (st[23], st[24], &wd, len[21]);
    VecDot  (st[24], st[24], &dd, len[21]);
    *test = SQRT(ww) * SQRT(dd);
    *test = ONE - wd / MAX(*test, EPS);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt
 * 
 * Purpose:  Generate map, metric, entropy from v = Tropus(w)
 * 
 * Notes: (1)  This routine contains the entire definitions of entropy and
 *             related quantities.  Entering with <1>=v and a functional form
 *             for S(f), the variational equation dS/df=-v defines the
 *             relation between f and v, enabling both f and S to be set.
 *             Then the metric is  [metric] = -1/(d2S/df2) , and we also need
 *             gradS = sqrt( dS/df.[metric].dS/df )  for alpha control.
 *        (2)  The entropy S must be uncorrelated and strictly convex, so that
 *             [metric] is diagonal and positive definite.  Also, S must have
 *             a global maximum of zero at v=0.
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      MKC/JS         20 Feb 1991     iEntropy = 5 added
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemEnt(
const MEMSYS* psMemsys, 
float*   s,       /*   O  Entropy */
float*   summet,  /*   O  SUM([metric]) */
float*   grads)   /*   O  gradS */
{
    int     mentrp = psMemsys->psInput->iEntropy;
    float   def    = psMemsys->psInput->fDef;
    float   a;
    float   b;
    float   sa;
    float   sb;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemEnt\n");
#endif
    if( mentrp == 1 ) 
        CALL( MemEnt1(psMemsys, def, s, grads, summet) )
    else if( mentrp == 2 )
        CALL( MemEnt2(psMemsys, def, s, grads, summet) )
    else if( mentrp == 3 )
        CALL( MemEnt3(psMemsys, def, s, grads, summet) )
    else if( mentrp == 4 ) 
        CALL( MemEnt4(psMemsys, def, s, grads, summet) )
    else if( mentrp == 5 ) 
    {
        CALL( MemEnt51(psMemsys, def, &sa, &sb) )
        if( sb > ZERO )
        {
            a = sa / sb;
            b = (a > ZERO) ? LOG(a) : ZERO;
        }
        else
        {
            a = ONE;
            b = ZERO;
        }
        CALL( MemEnt52(psMemsys, a, b, s, grads, summet) )
        *s -= sa;
    }
    *grads = SQRT(*grads);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt1
 *-----------------------------------------------------------------------------
 */
static int  MemEnt1(
const MEMSYS* psMemsys, 
float   def,
float*  s,
float*  gs,
float*  sum)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   a;
    float   c;

/* One block of standard entropy */
    if( def > ZERO ) 
    {
        VecFill(st[2], def, len[1]);
        VecSum (st[2], &a, len[1]);
        VecExp (st[1], st[2], len[1]);
        VecSmul(st[2], def, st[2], len[1]);
    }
    else
    {
        VecExp (st[1], st[2], len[1]);
        VecSum (st[3], &a, len[1]);
        VecMul (st[2], st[3], st[2], len[1]);
    } 
    VecDot (st[2], st[1], &c, len[1]);
    VecMul (st[1], st[1], st[1], len[1]);
    VecDot (st[2], st[1], gs, len[1]);
    VecSum (st[2], sum, len[1]);
    VecSqrt(st[2], st[1], len[1]);
    *s = *sum - a - c;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt2
 *-----------------------------------------------------------------------------
 */
static int  MemEnt2(
const MEMSYS* psMemsys, 
float   def,
float*  s,
float*  gs,
float*  sum)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   a;
    float   b1;
    float   b2;
    float   c1;
    float   c2;
    float   d1;
    float   d2;

/* One block of positive/negative entropy */
    if( def > ZERO ) 
    {
        VecFill(st[2], def, len[1]);
        VecSum (st[2], &a, len[1]);
        VecExp (st[1], st[2], len[1]);
        VecSmul(st[2], def, st[2], len[1]);
    }
    else 
    {
        VecExp (st[1], st[2], len[1]);
        VecSum (st[3], &a, len[1]);
        VecMul (st[2], st[3], st[2], len[1]);
    }
    VecSum (st[2], &b1, len[1]);
    VecMul (st[2], st[1], st[2], len[1]);
    VecSum (st[2], &c1, len[1]);
    VecDot (st[2], st[1], &d1, len[1]);
    VecSmul(st[1], -ONE, st[1], len[1]);
    VecExp (st[1], st[2], len[1]);
    if( def > ZERO ) 
        VecSmul(st[2], def, st[2], len[1]);
    else
        VecMul(st[2], st[3], st[2], len[1]);
    VecSum (st[2], &b2, len[1]);
    VecMul (st[2], st[1], st[2], len[1]);
    VecSum (st[2], &c2, len[1]);
    VecDot (st[2], st[1], &d2, len[1]);
    VecMov (st[1], st[2], len[1]);
    VecExp (st[2], st[1], len[1]);
    VecSmul(st[2], -ONE, st[2], len[1]);
    VecExp (st[2], st[2], len[1]);
    VecSub (st[2], st[1], st[2], len[1]);
    VecAdd (st[1], st[1], st[1], len[1]);
    VecAdd (st[1], st[2], st[1], len[1]);
    if( def > ZERO ) 
    {
        VecSmul(st[2], def, st[2], len[1]);
        VecSmul(st[1], def, st[1], len[1]);
    } 
    else
    {
        VecMul(st[2], st[3], st[2], len[1]);
        VecMul(st[1], st[3], st[1], len[1]);
    }
    VecSum (st[1], sum, len[1]);
    VecSqrt(st[1], st[1], len[1]);
    *s = b1 + b2 - a - a - c1 - c2;
    *gs = d1 + d2;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt3
 *-----------------------------------------------------------------------------
 */
static int  MemEnt3(
const MEMSYS* psMemsys, 
float   def,
float*  s,
float*  gs,
float*  sum)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   a;
    float   b;
    float   c;
    float   d;
    float   e;
    float   f;

/* One block of Fermi-Dirac entropy */
    if( def > ZERO ) 
    {
        a = LOG(def + ONE);
        VecFill(st[2], a, len[1]);
        VecSum (st[2], &a, len[1]);
        VecExp (st[1], st[2], len[1]);
        VecSmul(st[2], def, st[2], len[1]);
    } 
    else
    {
        VecSadd(st[3], ONE, st[2], len[1]);
        VecLog (st[2], st[2], len[1]);
        VecSum (st[2], &a, len[1]);
        VecExp (st[1], st[2], len[1]);
        VecMul (st[2], st[3], st[2], len[1]);
    }
    VecSadd (st[2], ONE, st[2], len[1]);
    VecRecip(st[2], st[2], len[1]);
    VecSum  (st[1], &b, len[1]);
    VecDot  (st[2], st[1], &c, len[1]);
    VecMul  (st[1], st[1], st[1], len[1]);
    VecMul  (st[2], st[1], st[1], len[1]);
    VecSum  (st[1], &d, len[1]);
    VecDot  (st[2], st[1], &e, len[1]);
    VecMov  (st[2], st[1], len[1]);
    VecLog  (st[1], st[2], len[1]);
    VecSum  (st[2], &f, len[1]);
    VecSadd (st[1], -ONE, st[2], len[1]);
    VecSmul (st[2], -ONE, st[2], len[1]);
    VecMul  (st[2], st[1], st[1], len[1]);
    VecSum  (st[1], sum, len[1]);
    VecSqrt (st[1], st[1], len[1]);
    *s = c - a - b - f;
    *gs = d - e;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt4
 *-----------------------------------------------------------------------------
 */
static int  MemEnt4(
const MEMSYS* psMemsys, 
float   def,
float*  s,
float*  gs,
float*  sum)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

/* One block of quadratic 'entropy' */
    if( def > ZERO )
    {
        VecSmul(st[1], def, st[2], len[1]);
        VecDot (st[2], st[1], gs, len[1]);
        VecFill(st[1], def, len[1]);
    }
    else
    {
        VecMul(st[1], st[3], st[2], len[1]);
        VecDot(st[2], st[1], gs, len[1]);
        VecMov(st[3], st[1], len[1]);
    }
    VecSum (st[1], sum, len[1]);
    VecSqrt(st[1], st[1], len[1]);
    *s = - (*gs) / TWO;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt51
 *-----------------------------------------------------------------------------
 */
static int  MemEnt51(
const MEMSYS* psMemsys, 
float   def,
float*  a,
float*  b)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;

/* One block of <2> = <3> exp<1>, with A=Sum<3>, B=SUM<2> */
    if( def > ZERO )
    {
        VecFill(st[2], def, len[1]);
        VecSum (st[2], a, len[1]);
        VecExp (st[1], st[2], len[1]);
        VecSmul(st[2], def, st[2], len[1]);
    } 
    else 
    {
        VecExp(st[1], st[2], len[1]);
        VecSum(st[3], a, len[1]);
        VecMul(st[2], st[3], st[2], len[1]);
    }
    VecSum(st[2], b, len[1]);
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemEnt52
 *-----------------------------------------------------------------------------
 */
static int  MemEnt52(
const MEMSYS* psMemsys, 
float   a,
float   b,
float*  s,
float*  gs,
float*  sum)
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    float   c;

/* One block of normalised entropy */
    VecSmul(st[2], a, st[2], len[1]);
    VecSadd(st[1], b, st[1], len[1]);
    VecDot (st[2], st[1], &c, len[1]);
    VecMul (st[1], st[1], st[1], len[1]);
    VecDot (st[2], st[1], gs, len[1]);
    VecSum (st[2], sum, len[1]);
    VecSqrt(st[2], st[1], len[1]);
    *s = *sum - c;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemChi
 * 
 * Purpose:  Generate residuals and likelihood from mock Gaussian data, for
 *                pr(Data|Mock) = exp( units - alhood ) ,
 *  
 *               <24>  := ( <21> - <25> ) * <22>
 *                      = ( Data - Mockdata ) /sigma = Scaled residuals
 *              alhood := <24>.<24> / 2          ( >= 0)
 *  
 *           Also generate dimensionalities
 *              data  := Number of nonzero accuracies
 *              units := Sum log ( nonzero accuracies )
 *  
 *           If nonlinear, rescale <23> to keep the map constant.
 * 
 * Notes: (1) This routine contains the entire definition of likelihood and
 *            related quantities, with zero likelihood minimum at F=D.
 * 
 * History:
 *      JS              7 Dec 1988     First MEMSYS3 version
 *      MKC/JS         19 Dec 1990     First release
 *      MKC            13 Feb 1991     Division by zero in <29> fixed
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemChi(
const MEMSYS* psMemsys, 
float*   alhood,  /*   O  L = - log likelihood */
float*   data,    /*   O  Number of data */
float*   units)   /*   O  Sum log accuracy */
{
    float   acc = psMemsys->psInput->fAcc;
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    int     (*Umemx)(float,float*,float*) = psMemsys->psProcs->pfnUmemx;
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemChi\n");
#endif
    if( Umemx )
    {
        VecMul   (st[23], st[29], st[23], len[21]);
        CALL( VecNonlin(st[25], st[24], st[29], len[21], Umemx) )
        VecMov   (st[24], st[25], len[21]);
        if( acc > ZERO )
            VecSmul(st[29], acc, st[29], len[21]);
        else 
            VecMul (st[29], st[22], st[29], len[21]);
        VecSadd (st[29], EPS, st[26], len[21]);
        VecDiv  (st[23], st[26], st[23], len[21]);
    }
    if( acc > ZERO )
    {
        VecSub  (st[21], st[25], st[24], len[21]);
        VecSmul (st[24], acc, st[24], len[21]);
        VecDot  (st[24], st[24], alhood, len[21]);
        *data = (float)len[21];
        *units = *data * LOG(acc);
    }
    else
    {
        VecSadd (st[22], EPS, st[26], len[21]);
        VecLog  (st[26], st[24], len[21]);
        VecDiv  (st[24], st[26], st[24], len[21]);
        VecDot  (st[24], st[22], units, len[21]);
        VecRecip(st[26], st[24], len[21]);
        VecDot  (st[24], st[22], data, len[21]);
        VecSub  (st[21], st[25], st[24], len[21]);
        VecMul  (st[24], st[22], st[24], len[21]);
        VecDot  (st[24], st[24], alhood, len[21]);
    }
    *alhood /= TWO;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemPoi
 * 
 * Purpose:  Generate residuals and likelihood from mock Poisson data, for
 *                pr(Data|Mock) = exp( units - alhood ) ,
 *  
 *               <22>  := sqrt(<21>+1) / <25>
 *                      = sqrt(gradgrad(-log(lhood))
 *                      = sqrt(data) / mockdata = accuracies
 *               <24>  := ( <21> - <25> ) / sqrt(<21>+1)
 *                      = - grad (-log Lhood) / sqrt(gradgrad(-log Lhood))
 *                      = Scaled residuals
 *              alhood := SUM ( <25> - <21> + <21> log(<21>/<25>) )
 *  
 *           Also generate dimensionalities
 *               data  := Number of nonzero accuracies
 *               units := Sum log ( nonzero accuracies )
 *  
 *            Rescale <23> to keep the map constant.
 * 
 * Notes: (1) This routine contains the entire definition of likelihood and
 *            related quantities, with zero likelihood minimum at F=D.
 *        (2) Proper operation of this routine requires positive mock data.
 *        (3) Adding 1 to D when setting Poisson accuracies can be justified
 *            by averaging gradgradL over a reasonable range of values of F.
 *        (4) The Poisson contribution ( D log D - D - log D! ) to
 *            UNITS is here approximated by  -(1/2) log( 1 + 2*pi*D ) .
 *            Importantly, error = 0 for D=0,  error -> 0 as D -> infinity.
 *            Max error = 0.007 (less than 1%) at D=1.
 * 
 * History:
 *      JS              7 Dec 1988     First MEMSYS3 version
 *      MKC/JS         19 Dec 1990     First release
 *      MKC            13 Feb 1991     Division by zero in <29> fixed
 *      JS             01 May 1999     Detailed revision
 *-----------------------------------------------------------------------------
 */
static int  MemPoi(
const MEMSYS* psMemsys, 
float*   alhood,  /*   O  L = - log likelihood */
float*   data,    /*   O  Number of data */
float*   units)   /*   O  Sum log accuracy */
{
    float** st  = psMemsys->psOutput->afSt;
    int*    len = psMemsys->psOutput->aiLen;
    int     (*Umemx)(float,float*,float*) = psMemsys->psProcs->pfnUmemx;
    float   RPI2 = ONE / (EIGHT * ATAN(ONE));
    int     CALLvalue = 0;

#if FLOWCHART
    printf("   MemPoi\n");
#endif
    if( Umemx )
    {
        VecMul   (st[23], st[29], st[23], len[21]);
        CALL( VecNonlin(st[25], st[24], st[29], len[21], Umemx) )
        VecMov   (st[24], st[25], len[21]);
    }
    else
    {
        VecMul  (st[23], st[22], st[23], len[21]);
    }
    VecDiv  (st[21], st[25], st[22], len[21]);
    VecSadd (st[22], EPS, st[22], len[21]);
    VecLog  (st[22], st[22], len[21]);
    VecMul  (st[21], st[22], st[22], len[21]);
    VecSub  (st[22], st[21], st[22], len[21]);
    VecAdd  (st[25], st[22], st[22], len[21]);
    VecSum  (st[22], alhood, len[21]);
    VecSadd (st[21], RPI2, st[22], len[21]);
    VecLog  (st[22], st[24], len[21]);
    VecSum  (st[24], units, len[21]);
    VecSadd (st[21], ONE, st[22], len[21]);
    VecSqrt (st[22], st[22], len[21]);
    VecSub  (st[21], st[25], st[24], len[21]);
    VecDiv  (st[24], st[22], st[24], len[21]);
    VecDiv  (st[22], st[25], st[22], len[21]);
    if( Umemx )
    {
        VecMul  (st[29], st[22], st[29], len[21]);
        VecSadd (st[29], EPS, st[26], len[21]);
        VecDiv  (st[23], st[26], st[23], len[21]);
    }
    else
    {
        VecDiv  (st[23], st[22], st[23], len[21]);
    }
    *data = (float)len[21];
    *units /= -TWO;
    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function: MemVrand
 * 
 * Purpose:  <m> := random vector
 *  
 *   If    ( ncorr < 0 ) then
 *                     components are uncorrelated random signs
 *   elseif( ncorr >= 0 ) then
 *                     components are from N(0,1),
 *                     with possible correlation with previous samples.
 *  
 *   During several successive calls to this routine with the same value of
 *   ncorr (psMemsys, strictly positive), let the m'th and n'th vector results
 *   be
 *       (m)                             (n)
 *      x    (i=0,1,2,....)     and     x    (i=0,1,2,....)
 *       i                               i
 *   All values are individually from N(0,1), but have covariance structure
 *                          -
 *         (m)  (n)        |  max( 1 - abs(m-n)/(ncorr+1) , 0 )    if  i=j
 *      < x    x   >   =   |
 *         i    j          |    0   if i not equal to j
 *                          -
 *   Examples: ncorr = 2  gives (for each i=j) covariance
 *  
 *                         | 3  2  1  0  0 |
 *         (m)  (n)        | 2  3  2  1  0 |
 *      < x    x   >   =   | 1  2  3  2  1 | / 3
 *                         | 0  1  2  3  2 |
 *                         | 0  0  1  2  3 |
 *  
 *             between the first five samples.
 *             ncorr = 0 would have given the uncorrelated identity.
 * 
 * Notes: (1) For correlated Gaussian, CPU time proportional to ncorr
 *        (2) In this implementation, the vector samples depend on addressing
 *            details such as the block-size.
 * 
 * History:
 *      MKC/JS         19 Dec 1990     First release
 *      JS             01 May 1999     Detailed revision (ncorr redefined)
 *-----------------------------------------------------------------------------
 */
static int  MemVrand(
const MEMSYS* psMemsys, 
int       m,     /* I    Output area number                   */
int       n,     /* I    Work area (not m if ncorr > 0)       */
int       ncorr) /* I    < 0  for random signs, or
                         = 0  for normal, uncorrelated
                         > 0  for normal, with correlation */
{
    float**   st    = psMemsys->psOutput->afSt;
    int*      len   = psMemsys->psOutput->aiLen;
    unsigned* krand = psMemsys->psOutput->iRand;
    unsigned  xrand[112];
    float     temp;
    int       i;

#if FLOWCHART
    printf("   MemVrand(%2d)\n", m);
#endif
/* Random signs */
    if( ncorr < 0 )
    {
        VecRand(krand, st[m], len[m]);
        VecSign(st[m], st[m], len[m]);
/* Normal deviates */
    }
    else if( ncorr == 0 )
    {
/* Independent samples */
        VecRand(krand, st[m], len[m]);
    }
    else
    {
/* Correlated samples, properly normalised. Areas m and n MUST BE DIFFERENT. */
        temp = SQRT(ONE / (ncorr + ONE));
        VecRand(krand, st[n], len[n]);
        VecSmul(st[n], temp, st[m], len[n]);
/* After the first pass, store the random generator */
        for( i = 0; i < 112; ++i )
            xrand[i]= krand[i];
        for( i = 0; i < ncorr; ++i ) 
        {
            VecRand (krand, st[n], len[n]);
            VecSmula(st[n], temp, st[m], st[m], len[n]);
        }
/* After all ncorr passes are done, recover the random generator so
 * that the last ncorr-1 samples will be repeated at the start of the 
 * next call.  This is the source of the correlation.
 */
        for( i = 0; i < 112; ++i )
            krand[i]= xrand[i];
    }
    return 0;
}
