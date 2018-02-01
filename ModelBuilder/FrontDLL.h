/**********************************************************************/
/* Frontline Systems Solver DLL (Dynamic Link Library) Version 4.5    */
/* Copyright 1995-2002 Frontline Systems, Inc.  All Rights Reserved.  */
/* Include this header file in your C or C++ program for MS Windows.  */
/**********************************************************************/

/*  Version 4.5 includes support for multiple, dynamically loaded,
	"dual use" large-scale Solver Engines that can be used with either
	the Premium Solver Platform V5.0 or the Solver DLL Platform 4.5.

    Version 3.5 includes support for Nonsmooth Solver engines, and
    API changes for multi-problem, recursive and multi-threaded use. */

/***************** VERSION COMPATIBILITY INFORMATION ******************/

/*  API FUNCTIONS:
  	The following API functions were added in Version 4.5:  
	loadengine()

	The following API functions were added in Version 3.5:  
	varstat(), constat(), lpread(), getuse(), reportuse()

    The following API functions (originally from V1.0) were deleted:
	setepgap(), getepgap(), setepopt(), getepopt() */

/*  SYMBOLIC NAMES:
	In Version 4.5, PROB_xxx symbols have numeric values different
	from Version 3.5.

	In Version 5.0, new "Solver Result" codes compatible with the
	Premium Solver Platform V5.0 will be returned for the "pstat" 
	argument of solution().  To maintain the same PSTAT_xxx return
	codes used in previous Solver DLL versions, call the API function
	setintparam( lp, PARAM_PSTAT, 1) */

/*  ARGUMENTS:
    In Version 3.5, all API functions except loadlp() and loadnlp() 
	take an HPROBLEM argument.  Hence, the following functions from 
	Version 3.0 have argument list changes: getproblimits(), 
	infointparam(), setintparam(), getintparam(), infodblparam(), 
	setdblparam(), getdblparam(), setdefaults(), setlpcallbackfunc(), 
	getlpcallbackfunc(), setmipcallbackfunc(), getmipcallbackfunc().
	Also, getlpcallbackfunc() and getmipcallbackfunc() return a value */

/**********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
   #define _FAR_
   #define _HUGE_
   #define _CC  __stdcall
#else /* WIN16 */
   #define _FAR_  _far
   #define _HUGE_ _huge
   #define _CC  __export __far __pascal
#endif

typedef long   INTARG;  /* same in Win16 and Win32 */
typedef double REALARG; /* same on all platforms...*/
typedef unsigned char _FAR_ * LPBYTEARG;
typedef INTARG  _FAR_ * LPINTARG;
typedef REALARG _FAR_ * LPREALARG;

/* Define huge pointers for int and REALARG types.  HPREALARG
   is used for matval[], and HPINTARG is used for matind[] to allow
   for more than 32768 nonzero coefficients, in 16-bit Windows. */
typedef INTARG  _HUGE_ * HPINTARG;
typedef REALARG _HUGE_ * HPREALARG;

#ifndef LPSTR
   #define LPSTR char _FAR_ *
#endif
#ifndef LPVOID
   #define LPVOID void _FAR_ *
#endif
#ifndef NULL
   #define NULL 0
#endif

/* Define the SAFEARRAY type for arrays passed from/to Visual Basic */
#ifndef __oaidl_h__
typedef struct  tagSAFEARRAYBOUND
    {
    unsigned long cElements;
    long lLbound;
    }	SAFEARRAYBOUND;
typedef struct  tagSAFEARRAY
    {
    unsigned short cDims;
    unsigned short fFeatures;
    unsigned long cbElements;
    unsigned long cLocks;
    void _FAR_ * pvData;
    SAFEARRAYBOUND rgsabound[ 1 ];
    }	SAFEARRAY;
typedef SAFEARRAY _FAR_ * LPSAFEARRAY;
#endif
typedef LPSAFEARRAY _FAR_ * PARAM_SAFEARRAY; /* VB passes this */

/* Define the value which represents "plus infinity" for bounds */
#define INFBOUND 1E30

/* Define problem type codes - linear, quadratic, nonlinear, nonsmooth, integer */
typedef enum {
   PROB_LP     =  8,
   PROB_MIP    =  10,
   PROB_QP     =  16,
   PROB_QPMIP  =  18,
   PROB_NLP    =  32,
   PROB_NLPMIP =  34,
   PROB_NSP    =  64,
   PROB_NSPMIP =  66
   } PROBTYPE;

/* Define problem structure codes - dense/sparse with and w/o occurrence flags */
typedef enum {
   STRUC_DENSE      = 0,
   STRUC_DENSEFLAG  = 1,
   STRUC_SPARSE     = 2,
   STRUC_SPARSEFLAG = 3
   } PROBSTRUC;

/* Define return codes for the "pstat" argument of solution() */
typedef enum {
   PSTAT_LIC_PROBLEM     = -1, /* Solver Platform or Engine license invalid*/
   PSTAT_CONTINUE        = 0,  /* Used in callback to continue solving     */
   PSTAT_OPTIMAL         = 1,  /* An optimal solution has been found       */
   PSTAT_INFEASIBLE      = 2,  /* No feasible solution could be found      */
   PSTAT_UNBOUNDED       = 3,  /* The objective value is unbounded         */
   PSTAT_RESTART         = 4,  /* Restart of solution process requested    */
   PSTAT_IT_LIM_FEAS     = 5,  /* Iteration limit exceeded, feasible sol   */
   PSTAT_IT_LIM_INFEAS   = 6,  /* Iteration limit exceeded, no feasible sol*/
   PSTAT_TIME_LIM_FEAS   = 7,  /* Time limit exceeded, feasible solution   */
   PSTAT_TIME_LIM_INFEAS = 8,  /* Time limit exceeded, no feasible solution*/
   PSTAT_ABORT_FEAS      = 12, /* User aborted solve, feasible solution    */
   PSTAT_ABORT_INFEAS    = 13, /* User aborted solve, no feasible solution */

   PSTAT_FRACT_CHANGE    = 65, /* Objective function changing too slowly   */
   PSTAT_NO_REMEDIES     = 66, /* All remedies failed to find better point */
   PSTAT_FLOAT_ERROR     = 67, /* Error when evaluating problem functions  */
   PSTAT_MEM_LIM         = 68, /* Could not allocate enough main memory    */
   PSTAT_ENTRY_ERROR     = 69, /* Attempt to use DLL from 2 or more threads*/
   PSTAT_PROBABLE        = 70, /* MSL (Bayesian) global opt test satisfied */
   PSTAT_BOUNDS_MISSING  = 71, /* Bounds missing for EV/MSL Require Bounds */
   PSTAT_BOUNDS_CONFLICT = 72, /* <= = >= bounds conflict with bin/alldif  */
   PSTAT_BOUNDS_INCONSIST= 73, /* lb[i] > ub[i] for some variable bound i  */
   PSTAT_DERIVATIVE_ERR  = 74, /* jacobian func could not compute gradients*/

   PSTAT_MIP_OPTIMAL         = 101,/* Optimal integer solution found       */
   PSTAT_MIP_OPTIMAL_TOL     = 102,/* Integer solution found within epgap  */
   PSTAT_MIP_INFEASIBLE      = 103,/* No feasible integer solution         */
   PSTAT_MIP_SOL_LIM         = 104,/* Limit of integer solutions exceeded  */
   PSTAT_MIP_NODE_LIM_FEAS   = 105,/* Node limit exceeded, integer solution*/
   PSTAT_MIP_NODE_LIM_INFEAS = 106,/* Node limit exceeded, no int solution */
   PSTAT_MIP_TIME_LIM_FEAS   = 107,/* Time limit exceeded, integer solution*/
   PSTAT_MIP_TIME_LIM_INFEAS = 108,/* Time limit exceeded, no int solution */

   PSTAT_LINEAR_INVALID = 200, /* Linear model is not a valid assumption   */
   PSTAT_BAD_DATASET    = 201, /* Inconsistent dataset - contact Frontline */
   PSTAT_API_ERROR      = 202, /* Internal API error - contact Frontline   */
   PSTAT_API_DEFAULT    = 203  /* Internal API error - contact Frontline   */
   } PSTAT;

#define PSTAT_USER_ABORT  PSTAT_ABORT_FEAS  /* may be returned by callback */

/* Define symbolic names for parameters controlled by set/get/infointparam
   and set/get/infodblparam.  When PARAM_ARGCK, PARAM_ARRAY, and PARAM_PSTAT 
   are used with setintparam() and getintparam(), the HPROBLEM argument is
   ignored and the effect is global, on all problem instances and threads */
typedef enum {
   PARAM_ARGCK      =  990, /* on argument errors: 1-MsgBox, 0-retval only*/
   PARAM_PSTAT      =  993, /* use PSTAT_xxx codes in solution() pstat arg*/
   PARAM_ARRAY      =  995, /* arrays are: 0-C-style arrays, 1-SAFEARRAYs */
   PARAM_USERP      =  997, /* 0-eval/test, 1-auto report, 2-user reports */
   PARAM_IISBND     =  999, /* IIS finder: 0-includes bnds, 1-omits bounds*/
   PARAM_ITLIM      = 1020, /* limit on LP iterations - per (sub)problem  */
   PARAM_NDLIM      = 2017, /* limit on Branch & Bound nodes explored     */
   PARAM_MIPLIM     = 2015, /* limit on Branch & Bound IP solutions found */
   PARAM_SAMPSZ     = 6001, /* Evol Solver population (sample) size       */
   PARAM_NOIMP      = 6002, /* Evol Solver max time w/no improvement      */
   PARAM_SCAIND     = 1033, /* scaling: -1-none, 0-normal, 1-aggressive   */
   PARAM_CRAIND     = 1007, /* crashing: 0-none, 1-crash initial basis    */
   PARAM_SUBALG     = 2026, /* B&B LP subprob algorithm: 1-Primal, 2-Dual */
   PARAM_RELAX      = 2501, /* solve relaxation: 1-on, 0-use int vars     */
   PARAM_PREPRO     = 2502, /* B&B Preprocessing and Probing: 1-on, 0-off */
   PARAM_OPTFIX     = 2503, /* B&B Optimality Fixing:  1-on, 0-off        */
   PARAM_REORDR     = 2504, /* B&B Branch Variable Reordering: 1-on, 0-off*/
   PARAM_IMPBND     = 2505, /* B&B Bounds Improvement: 1-on, 0-off        */
   PARAM_REQBND     = 6003, /* EV Solver vars: 1-require bounds, 0-don't  */
   PARAM_NSFEAS     = 6004, /* EV Solver feasibility: 1-strict, 0-penalty */
   PARAM_NSLOCL     = 6005, /* EVLS:0-Random, 1-Direct, 2-Nonlin, 3-Linear*/
   PARAM_RSEED      = 6006, /* EV Solver/MSL random number generator seed */

   PARAM_TILIM      = 1038, /* global solution time limit  */
   PARAM_EPOPT      = 1014, /* LP optimality tolerance     */
   PARAM_EPPIV      = 1091, /* LP pivot tolerance          */
   PARAM_EPSOL      = 1092, /* LSLP solution tolerance     */
   PARAM_EPRHS      = 1016, /* LP feasibility tolerance    */
   PARAM_EPGAP      = 2009, /* integer tolerance/MIP gap   */
   PARAM_CUTLO      = 2006, /* known incumbent for max MIP */
   PARAM_CUTHI      = 2007, /* known incumbent for min MIP */

   PARAM_EPNEWT     = 5001, /* NLP constraint tolerance    */
   PARAM_EPCONV     = 5002, /* NLP slow change stopping tol*/
   PARAM_EPSTEP     = 5003, /* NLP differencing step size  */
   PARAM_LINVAR     = 5010, /* NLP recognize linear vars   */
   PARAM_DERIV      = 5011, /* NLP derivative computation  */
   PARAM_ESTIM      = 5012, /* NLP estimates for basic vars*/
   PARAM_DIREC      = 5013, /* NLP search direction option */
   PARAM_MULTI      = 5501, /* NLP global multistart search*/
   PARAM_MTOPO      = 5502, /* NLP global topographicsearch*/
   PARAM_MUTATE     = 6010, /* Evol mutation probability   */
   PARAM_EXTINC     = 6011, /* Evol extinction percentage  */

   /* Parameters unique to "dual use" large-scale Solver Engines  */

   PARAM_LSLP_CRASH   =   1007, /* Same as PARAM_CRAIND, see above*/

   PARAM_LSSQP_LINCON =  11001, /* Treat constraints as linear    */
   PARAM_LSSQP_LINOBJ =  11002, /* Treat objective as linear      */

   PARAM_LSGRG_RLXBND =  12001, /* Relax bounds req. on variables */

   PARAM_OPTQ_CHKDUP  =  13001, /* Check for duplicated solutions */
   PARAM_OPTQ_TAGUCH  =  13002, /* Use Taguchi design experiments */
   PARAM_OPTQ_OBJPRE  =  13003, /* Precision of objective function*/
   PARAM_OPTQ_VARPRE  =  13004, /* Precision of decision variables*/
   PARAM_OPTQ_BNDFRQ  =  13005, /* Boundary frequency for search  */

   PARAM_LGO_GOITER   =  14001, /* Global phase maximum iterations*/
   PARAM_LGO_GOITWO   =  14002, /* Global phase max iters w/o impr*/
   PARAM_LGO_GOCUT    =  14003, /* Global search objective cutoff */
   PARAM_LGO_LOCUT    =  14004, /* Local search objective cutoff  */
   PARAM_LGO_GOCONV   =  14005, /* Global search convergence toler*/
   PARAM_LGO_GOMETH   =  14006, /* Global search method: */
   /*0-local search, 1-global B&B, 2-adaptive random, 3-multistart*/

   } PARAM;

/* Define "wherefrom" codes for callback functions  */
typedef enum {
   CALLBACK_PRIMAL  = 1,      /* callback on LP/QP pivot   */
   CALLBACK_NLP     = 11,     /* callback on NLP/NSP iter  */
   CALLBACK_MIP     = 101     /* callback on MIP branch    */
   } CALLBACK_CODE;

/* Define info request codes for callback functions */
typedef enum {
   CBINFO_PRIMAL_OBJ     = 1, /* callback: objective       */
   CBINFO_PRIMAL_INFMEAS = 3, /* callback: infeasibility   */
   CBINFO_PRIMAL_FEAS    = 5, /* callback: whether feasible*/
   CBINFO_ITCOUNT        = 7, /* callback: iteration count */
   CBINFO_BEST_INTEGER   = 8, /* callback: MIP incumbent   */
   CBINFO_NODE_COUNT     = 10,/* callback: nodes explored  */
   CBINFO_INCUMB_COUNT   = 11,/* callback: MIP solutions   */
   CBINFO_MIP_ITERATIONS = 12 /* callback: MIP iterations  */
   } CBINFO;


typedef struct PROBLEM _FAR_ * HPROBLEM;

/* Define the procedure pointer type for callback functions */
typedef INTARG (_CC * _CCPROC) (HPROBLEM lpinfo, INTARG wherefrom);

/* Define the callbacks to compute function values, and the objective
   value, objective gradient and Jacobian matrix of partial derivatives */
typedef INTARG (_CC * _FUNCEVAL) (HPROBLEM lp, INTARG numcols,
   INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG var,
   INTARG varone, INTARG vartwo);

typedef INTARG (_CC * _JACOBIAN) (HPROBLEM lp, INTARG numcols,
   INTARG numrows, INTARG nzspace, LPREALARG objval, LPREALARG obj,
   LPINTARG matbeg, LPINTARG matcnt, HPINTARG matind, HPREALARG matval,
   LPREALARG var, LPBYTEARG objtype, LPBYTEARG matvaltype);

/* These callback typedefs are useful for Visual Basic.
   They pass OLE SAFEARRAYs rather than C-style arrays */
typedef INTARG (_CC * _SAFUNCEVAL) (HPROBLEM lp, INTARG numcols,
   INTARG numrows, LPREALARG objval, LPSAFEARRAY _FAR_ * lhs, 
   LPSAFEARRAY _FAR_ * var, INTARG varone, INTARG vartwo);

typedef INTARG (_CC * _SAJACOBIAN) (HPROBLEM lp, INTARG numcols,
   INTARG numrows, INTARG nzspace, LPREALARG objval, LPSAFEARRAY _FAR_ * obj,
   LPSAFEARRAY _FAR_ * matbeg, LPSAFEARRAY _FAR_ * matcnt, 
   LPSAFEARRAY _FAR_ * matind, LPSAFEARRAY _FAR_ * matval,
   LPSAFEARRAY _FAR_ * var, LPSAFEARRAY _FAR_ * objtype,
   LPSAFEARRAY _FAR_ * matvaltype);


/* The problem structure:  loadlp() and loadnlp() load their arguments
   into a PROBLEM structure and return an HPROBLEM pointer to this data
   structure.  Other API calls take this HPROBLEM pointer as an argument.
   unloadprob() frees memory and clears addresses in the structure members.
*/

#ifdef __cplusplus
   class FRONTAPI;
#endif

typedef struct PROBLEM
   {
   struct PROBLEM _FAR_ * lp; /* points to itself, for compatibility */
   LPSTR probname;
   INTARG probtype, probstruc, numcols, numrows, numints, objsen;
   INTARG linearvars, linearfcns, smoothvars, smoothfcns;
   LPREALARG objx, rhsx;
   LPBYTEARG sense;
   LPINTARG matbeg, matcnt;
   HPINTARG matind;
   HPREALARG matval;
   LPREALARG var, lb, ub, rngval;
   INTARG colspace, rowspace, nzspace, nzobj;
   LPINTARG qmatbeg, qmatcnt;
   HPINTARG qmatind;
   HPREALARG qmatval;
   INTARG qnzspace;
   LPBYTEARG ctype, objtype, matvaltype, fcntype, vartype;
   INTARG diffgrps, diffvars;
   LPINTARG diffindex, different;
   REALARG objective, sum_of_infeas, best_integer;
   INTARG iterations, branches, total_iters, integer_soln;
   INTARG stat;
   REALARG objval;
   LPREALARG x, piout, slack, dj;
   struct {
      INTARG begin, end;
      LPREALARG std, lower, upper;
      } obj;
   struct {
      INTARG begin, end;
      LPREALARG std, lower, upper;
      } rhs;
   struct {
      INTARG numrows, numcols, stat, write;
      LPINTARG rowind, rowbdstat, colind, colbdstat;
      } iis;
   struct {
      PARAM_SAFEARRAY objx, matbeg, matcnt, matind, matval;
	  PARAM_SAFEARRAY objtype, matvaltype;
      } sa;
   _FUNCEVAL funceval;
   _JACOBIAN jacobian;
#ifdef __cplusplus
   FRONTAPI _FAR_ *lpAPI;
#else
   LPVOID lpAPI;
#endif
   }
   PROBLEM;


/* Prototypes for the DLL's callable entry points... */


/* Define the entry points for the Linear & Quadratic Solver */
HPROBLEM _CC loadlp (LPSTR probname,
   INTARG numcols, INTARG numrows, INTARG objsen,
   LPREALARG obj, LPREALARG rhs, LPBYTEARG sense,
   LPINTARG matbeg, LPINTARG matcnt, HPINTARG matind, HPREALARG matval,
   LPREALARG lb, LPREALARG ub, LPREALARG rngval,
   INTARG colspace, INTARG rowspace, INTARG nzspace);
   
INTARG _CC loadquad (HPROBLEM lp, LPINTARG qmatbeg, LPINTARG qmatcnt,
   HPINTARG qmatind, HPREALARG qmatval, INTARG qnzspace, LPREALARG var);

INTARG _CC loadctype (HPROBLEM lp, LPBYTEARG ctype); /* MIP problems */

/* Define the entry points for the Nonlinear & NonSmooth Solver */
HPROBLEM _CC loadnlp (LPSTR probname,
   INTARG numcols, INTARG numrows, INTARG objsen,
   LPREALARG obj, LPREALARG rhs, LPBYTEARG sense,
   LPINTARG matbeg, LPINTARG matcnt, HPINTARG matind, HPREALARG matval,
   LPREALARG var, LPREALARG lb, LPREALARG ub, LPREALARG rngval,
   INTARG nzspace, _FUNCEVAL funceval, _JACOBIAN jacobian);

INTARG _CC loadnltype (HPROBLEM lp, LPBYTEARG objtype, LPBYTEARG matvaltype);

INTARG _CC testnltype (HPROBLEM lp, INTARG numtests, LPREALARG testvals,
					   LPINTARG pstat, LPBYTEARG objtype, LPBYTEARG matvaltype);

/* Define the optimize and solution routines */
INTARG _CC optimize (HPROBLEM lp);

INTARG _CC mipoptimize (HPROBLEM lp);

INTARG _CC solution (HPROBLEM lp, LPINTARG pstat, LPREALARG pobj,
   LPREALARG x, LPREALARG piout, LPREALARG slack, LPREALARG dj);

INTARG _CC objsa (HPROBLEM lp, INTARG begidx, INTARG endidx,
   LPREALARG lower, LPREALARG upper);

INTARG _CC rhssa (HPROBLEM lp, INTARG begidx, INTARG endidx,
   LPREALARG lower, LPREALARG upper);

/* varstat() and constat() are new for NonSmooth Solver engines */
INTARG _CC varstat (HPROBLEM lp, INTARG begidx, INTARG endidx,
   LPREALARG mid, LPREALARG disp, LPREALARG lower, LPREALARG upper);

INTARG _CC constat (HPROBLEM lp, INTARG begidx, INTARG endidx,
   LPREALARG mid, LPREALARG disp, LPREALARG lower, LPREALARG upper);

/* Define the IIS (Irreducibly Infeasible Set) finding routines */
INTARG _CC findiis (HPROBLEM lp, LPINTARG iisnumrows_p, LPINTARG iisnumcols_p);

INTARG _CC getiis (HPROBLEM lp, LPINTARG iisstat_p,
   LPINTARG rowind, LPINTARG rowbdstat, LPINTARG iisnumrows_p,
   LPINTARG colind, LPINTARG colbdstat, LPINTARG iisnumcols_p);

/* unloadprob MUST be called to free memory */
INTARG _CC unloadprob (HPROBLEM _FAR_ *lp_p);

/* routines to set and get parameter values, including size limits.
   Note- getproblimits(), getuse(), reportuse(), infointparam(), and 
   infodblparam() always ignore their HPROBLEM argument; setintparam()
   and getintparam() ignore their HPROBLEM argument when they are
   called with PARAM_ARGCK, PARAM_ARRAY and PARAM_USERP parameters.
   In all other situations, you must first define a problem with 
   loadlp() or loadnlp(), and pass its HPROBLEM pointer to the 
   other routines. */

INTARG _CC loadengine (HPROBLEM lp, INTARG type, LPSTR enginename);

INTARG _CC getproblimits (HPROBLEM lp, INTARG type,
   LPINTARG numcols_p, LPINTARG numrows_p, LPINTARG numints_p);

INTARG _CC getuse (HPROBLEM lp, LPINTARG loadprob_p, LPINTARG optimize_p,
   LPINTARG verify_p, LPINTARG repload_p, LPINTARG repopt_p, LPINTARG repdate_p);

INTARG _CC reportuse (HPROBLEM lp, LPSTR probname, LPSTR filename,
   LPSTR profilename, LPSTR password);

INTARG _CC infointparam (HPROBLEM lp, INTARG whichparam,
   LPINTARG defvalue_p, LPINTARG minvalue_p, LPINTARG maxvalue_p);

INTARG _CC infodblparam (HPROBLEM lp, INTARG whichparam,
   LPREALARG defvalue_p, LPREALARG minvalue_p, LPREALARG maxvalue_p);

INTARG _CC setintparam (HPROBLEM lp, INTARG whichparam, INTARG newvalue);

INTARG _CC getintparam (HPROBLEM lp, INTARG whichparam, LPINTARG value_p);

INTARG _CC setdblparam (HPROBLEM lp, INTARG whichparam, REALARG newvalue);

INTARG _CC getdblparam (HPROBLEM lp, INTARG whichparam, LPREALARG value_p);

INTARG _CC setdefaults (HPROBLEM lp);

/* routines to define and use the LP and MIP callback functions */

INTARG _CC setlpcallbackfunc (HPROBLEM lp, _CCPROC callback);

INTARG _CC getlpcallbackfunc (HPROBLEM lp, _CCPROC *callback_p);

INTARG _CC setnlpcallbackfunc (HPROBLEM lp, _CCPROC callback);

INTARG _CC getnlpcallbackfunc (HPROBLEM lp, _CCPROC *callback_p);

INTARG _CC setmipcallbackfunc (HPROBLEM lp, _CCPROC callback);

INTARG _CC getmipcallbackfunc (HPROBLEM lp, _CCPROC *callback_p);

INTARG _CC getcallbackinfo (HPROBLEM lpinfo, INTARG wherefrom,
   INTARG infonumber, void _FAR_ *result_p);

/* routines to read and write a file summarizing an LP/QP problem in
   algebraic form */

INTARG _CC lpread (HPROBLEM lp, LPSTR filename, LPINTARG objsen_p,
   LPINTARG numcols_p, LPINTARG numrows_p, LPINTARG numints_p,
   LPINTARG matcnt, LPINTARG qmatcnt);

INTARG _CC lpwrite (HPROBLEM lp, LPSTR filename);

INTARG _CC lprewrite (HPROBLEM lp, LPSTR filename);

INTARG _CC iiswrite (HPROBLEM lp, LPSTR filename);

#ifdef __cplusplus
}
#endif
