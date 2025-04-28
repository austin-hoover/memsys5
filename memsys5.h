/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Project:   MemSys5  Quantified Maximum Entropy System.
 * 
 * Filename:  memsys5.h
 * 
 * Purpose:   Header for memsys5.c version 1.41
 *            MemSys structure type definitions for
 *            MEMPTR  MEMIN  MEMOUT  MEMPROCS  MEMSYS,
 *            MemSys function prototypes, and return codes.
 *
 * Notes:     The MemSys5 structures use Hungarian notation, in which names
 *            have initial capitalisation, preceded by a prefix specifying the
 *            variable type :-
 *                    p = pointer to         a = array of
 *                    s = structure          f = floating point
 *                    i = int                c = character
 *                  pfn = function pointer  fp = file pointer
 *
 *            The MemSys5 functions all return type "int".
 *            A return value of 0 is "success".
 *            A positive return value is informative.
 *            A negative return value marks an error state.
 * 
 * History:   21 Feb 1992, 01 May 1999
 * 
 *            Copyright (c) 1992-1999 Maximum Entropy Data Consultants Ltd.
 *-----------------------------------------------------------------------------
 */
#ifndef MEMSYS5H
#define MEMSYS5H

/**************/
/* Structures */
/**************/

struct sMEMPTR;
struct sMEMIN;
struct sMEMOUT;
struct sMEMPROCS;

typedef struct
{
    struct sMEMPTR*     psVpoint;
    struct sMEMIN*      psInput;
    struct sMEMOUT*     psOutput;
    struct sMEMPROCS*   psProcs;
    FILE*               fpLog;
} MEMSYS;                 /* Global structure with all MemSys structures */

typedef struct sMEMPTR
{
    int      iNcell;      /* I    # result cells                             */
    float*   afResult;    /*   O  MaxEnt result                     [iNcell] */
    float*   afDefault;   /*(I)   Model (if psInput->fDef = 0)      [iNcell] */
    float*   afMask;      /*(I)   Quantify mask (for MemMask5)      [iNcell] */

    int      iNdata;      /* I    # data                                     */
    float*   afData;      /* I    Data                              [iNdata] */
    float*   afAcc;       /*(I)   Accuracies (if psInput->fAcc = 0) [iNdata] */
    float*   afLagrange;  /*   O  Lagrange multipliers for result   [iNdata] */
    float*   afMock;      /*   O  Mock data (for memrun=4)          [iNdata] */
    float*   afResiduals; /*   O  Acc * (Data - Mock)               [iNdata] */
    float*   afNonlinear; /*  (O) (required if psProcs->Umemx)      [iNdata] */

    void*    psUser;      /* I    Transform information                      */
} MEMPTR;                 /* MemSys vector pointers                          */

typedef struct sMEMIN
{
    int      iBayes;      /* I    Stopping switch    */
    int      iEntropy;    /* I    Entropy switch     */
    int      iGauss;      /* I    Likelihood switch  */
    int      iNrand;      /* I    # random vecs      */
    int      iIseed;      /* I    Regulariser seed   */
    float    fAim;        /* I    Constraint target  */
    float    fRate;       /* I    Dist & beta cntrl  */
    float    fDef;        /* I    Default (if >0)    */
    float    fAcc;        /* I    Accuracies (if >0) */
    float    fUtol;       /* I    Tolerance, O(0.1)  */
} MEMIN;                  /* MemSys input parameters */

typedef struct sMEMOUT
{
    float    fS;          /*   O  Entropy            */
    float    fTest;       /*   O  Gradient misfit    */
    float    fChisq;      /*   O  Chisquared         */
    float    fScale;      /*   O  Scaling sigma(D)   */
    float    fPlow;       /*   O  < evidence         */
    float    fPhigh;      /*   O  > evidence         */
    float    fPdev;       /*   O  Std dev evidence   */
    float    fGlow;       /*   O  < good             */
    float    fGhigh;      /*   O  > good             */
    float    fGdev;       /*   O  Std dev good       */
    float    fOmega;      /*   O  Termination        */
    float    fAlpha;      /* I O  Regularising parm  */
    int      iNtrans;     /* I O  # of transforms    */

    float    fHhigh;      /* I O  Internal MemSys bubble mismatch  */
    float    fXtable[8];  /* I O  Internal MemSys log(alpha)       */
    float    fYtable[8];  /* I O  Internal MemSys log(omega)       */
    float    fVtable[8];  /* I O  Internal MemSys variance(omega)  */
    int      iNtable;     /* I O  Internal MemSys length of table  */
    int      iIstat;      /* I O  Internal MemSys return code      */
    unsigned iRand[112];  /* I O  Internal MemSys randomiser state */
    float*   afSt[41];    /*  (O) Internal MemSys addresses        */
    int      aiLen[41];   /*  (O) Internal Memsys lengths          */
} MEMOUT;                 /* MemSys output parameters              */

typedef struct sMEMPROCS
{
    int (*pfnOpus)                            /* Object-to-Data transform  */
                  (       void*   psUser,     /* I    Transform structure  */
                    const float*  afObject,   /* I    Input MaxEnt area    */
                          float*  afData );   /*   O  Output data  area    */

    int (*pfnTropus)                          /* Transpose opus            */
                  (       void*   psUser,     /* I    Transform structure  */
                    const float*  afData,     /* I    Input  data   area   */
                          float*  afObject ); /*   O  Output MaxEnt area   */

    int (*pfnUmemx)                           /* Data nonlinearity         */
                  (       float   fX,         /* I    Linear transform (x) */
                          float*  fPhi,       /*   O  Nonlinear mock F(x)  */
                          float*  fDphi );    /*   O  Differential dF/dx   */

    int (*pfnUinfo)                           /* Manage diagnostics        */
                  (       FILE*   fpLog,      /* I    File for diagnostics */
                    const char*   cString);   /* I    Diagnostic string    */

} MEMPROCS;               /* MemSys procedure calls                        */

/**************/
/* Prototypes */
/**************/

extern int  Mem5                         /* MaxEnt iterate    */
              ( const MEMSYS* psMemsys,  /* Global structure  */
                int           iMemrun ); /* Start or continue */

extern int  MemMovie5                    /* Movie sample      */
              ( const MEMSYS* psMemsys,  /* Global structure  */
                int           iNcorr );  /* Sample overlap    */

extern int  MemMask5                     /* Quantify MaxEnt mask   */
              ( const MEMSYS* psMemsys,  /* Global structure       */
                float*        fMean,     /* Mean of required sum   */
                float*        fErrlo,    /* Lower limit on std dev */
                float*        fErrhi );  /* Upper limit on std dev */

extern int  MemTrq                       /* Transpose consistency check */
              ( const MEMSYS* psMemsys,  /* Global structure            */
                int           iSeed,     /* Random seed                 */
                float*        fErr );    /* Dimensionless inconsistency */

/****************/
/* Return codes */
/****************/

#define E_MEM_INPUT          -131
#define E_MALLOC             -130

#endif
