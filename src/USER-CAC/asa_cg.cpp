/* =========================================================================
   ============================== ASA_CG ===================================
   =========================================================================
       ________________________________________________________________
      | A Conjugate Gradient (cg_descent) based Active Set Algorithm   |
      |                                                                |
      |                  Version 1.0 (May 18, 2008)                    |
      |                  Version 1.1 (June 29, 2008)                   |
      |                  Version 1.2 (September 10, 2008)              |
      |                  Version 1.3 (September 25, 2009)              |
      |                  Version 2.0 (March 30, 2011)                  |
      |                  Version 2.1 (April 4, 2011)                   |
      |                  Version 2.2 (April 14, 2011)                  |
      |                  Version 3.0 (September 19, 2013)              |
      |                                                                |
      |        William W. Hager    and      Hongchao Zhang             |
      |        hager@math.ufl.edu         hzhang@math.ufl.edu          |
      |  Department of Mathematics      Department of Mathematics      |
      |     University of Florida       Louisiana State University     |
      |  Gainesville, Florida 32611      Baton Rouge, Louisiana        |
      |     352-392-0281 x 244                                         |
      |                                                                |
      |      Copyright by William W. Hager and Hongchao Zhang          |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |________________________________________________________________|
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|*/

#include "asa_user.h"
#include "asa_cg.h"
#include "cg_blas.h"


using namespace LAMMPS_NS;
/* begin external variables */
double one [1], zero [1] ;
BLAS_INT blas_one [1] ;
/* end external variables */
int asa_cg /*  return:
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f| in cg)
                       2 (cg iterations in all passes or
                          in one pass exceeded their limit)
                       3 (slope always negative in line search in cg)
                       4 (number secant iterations exceed nsecant in cg)
                       5 (search direction not a descent direction in cg)
                       6 (excessive updating of eps in cg)
                       7 (Wolfe conditions never satisfied in cg)
                       8 (debugger is on and the function value increases in cg)
                       9 (no cost or gradient improvement in
                          2n + Parm->nslow iterations)
                      10 (out of memory in cg)
                      11 (function nan or +-INF and could not be repaired in cg)
                      12 (invalid choice for memory parameter for cg)
                      13 (out of memory)
                      14 (cbb iterations in all passes or
                          in one pass exceeded their limit)
                      15 (line search failed in cbb iteration)
                      16 (search direction in cbb is not descent direction)
                      17 (function value became nan in cbb)*/
(
    double            *x, /* input: starting guess, output: the solution */
    double           *lo, /* lower bounds */
    double           *hi, /* upper bounds */
    ASA_INT                n, /* problem dimension */
    asa_stat       *Stat, /* structure with statistics (can be NULL) */
    asacg_parm    *CParm, /* user parameters, NULL = use default parameters */
    asa_parm      *AParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* |Proj (x_k - g_k) - x_k|_inf <= grad_tol */

    double (*valgrad) (asa_objective *), /* function and gradient
                                            NULL = use value & grad routines */
    double         *Work, /* NULL => allocate real work space
                             if DynamicMemory = FALSE
                                 LBFGS = 1  => need 2*mem*(n+1) + 5*n + m
                                 memory > 0 => need (mem+7)*n +(3*mem+9)*mem+5+m
                                               mem = MIN(memory, n)
                                 memory = 0 => need 5*n + m
                            if DymamicMemory = TRUE, need  5n + m */
    ASA_INT          *iWork,   /* NULL => allocate integer work space
                             otherwise provide space to n integers */
    LAMMPS_NS::PairCAC* objpoint,
    LAMMPS_NS::AtomVecCAC* avec_objpoint
)
{
    int gp, ident, status, mem ;
    ASA_INT i, j, nfree, cbb_totit, cg_totit, *ifree ;
    double alpha, gj, pert_lo, pert_hi, t, tl, th, gnorm, ginorm, pgnorm, xnorm,
           xj, xg, xp, *work, *d, *g, *xtemp, *gtemp, *pg ;
    asacg_parm *cgParm, cgParmStruc ;
    asa_parm *asaParm, asaParmStruc ;
    asa_com Com ;
    

    /* assign values to the external variables */
    one [0] = (double) 1 ;
    zero [0] = (double) 0 ;
    blas_one [0] = (BLAS_INT) 1 ;

/* initialize the parameters */

    if ( CParm == NULL )
    {
        cgParm = &cgParmStruc ;
        asa_cg_default (cgParm) ;
    }
    else cgParm = CParm ;
    if ( cgParm->PrintParms ) asa_printcgParms (cgParm) ;

    if ( AParm == NULL )
    {
        asaParm = &asaParmStruc ;
        asa_default (asaParm) ;
    }
    else asaParm = AParm ;

    if ( asaParm->PrintParms ) asa_printParms (asaParm) ;

    /* abort after maxit iterations of cbb in one pass */
    if ( asaParm->maxit_fac == ASA_INF ) Com.pgmaxit = ASA_INT_INF ;
    else Com.pgmaxit = (ASA_INT) (((double) n)*asaParm->maxit_fac) ;

    /* abort after totit iterations of cbb in all passes */
    if ( asaParm->totit_fac == ASA_INF ) cbb_totit = ASA_INT_INF ;
    else cbb_totit = (ASA_INT) (((double) n)*asaParm->totit_fac) ;

    /* abort after maxfunc function evaluation in one pass of cbb */
    if ( asaParm->maxfunc_fac == ASA_INF ) Com.pgmaxfunc = ASA_INT_INF ;
    else Com.pgmaxfunc = (ASA_INT) (((double) n)*asaParm->maxfunc_fac) ;

    cg_totit = cgParm->maxit ;

    Com.eps = cgParm->eps ;
    Com.neps = 0 ;
    Com.PertRule = cgParm->PertRule ;

    pert_lo = asaParm->pert_lo ;
    pert_hi = asaParm->pert_hi ;
    if(objpoint!=NULL){
    Com.user = objpoint->Objective ;
	(objpoint->Objective)->n = n ;
    }
    if(avec_objpoint!=NULL){
    Com.user = avec_objpoint->Objective ;
	(avec_objpoint->Objective)->n = n ;    
    }
    Com.tau1 = asaParm->tau1 ;
    Com.tau2 = asaParm->tau2 ;

    Com.cgParm = cgParm ;
    Com.asaParm = asaParm ;
    Com.x = x ;
    Com.n = n ;             /* problem dimension */
    Com.n5 = n % 5 ;
    Com.nf = (ASA_INT) 0 ;      /* number of function evaluations */
    Com.ng = (ASA_INT) 0 ;      /* number of gradient evaluations */
    Com.cbbiter = (ASA_INT) 0 ; /* number of cbb iterations evaluations */
    Com.cgiter = (ASA_INT) 0 ;  /* number of cg iterations */
    Com.AWolfe = cgParm->AWolfe ; /* do not touch user's AWolfe */
    Com.AArmijo = asaParm->AArmijo ; /* do not touch user's AArmijo */

    Com.valgrad = NULL ;
    Com.DimReduce = FALSE ;
    Com.objpoint=objpoint;
    Com.avec_objpoint=avec_objpoint;
    /* allocate integer work array */
    if ( iWork == NULL )
    {
        ifree = Com.ifree = (ASA_INT *) malloc (n*sizeof (ASA_INT)) ;
    }
    else
    {
        ifree = Com.ifree = iWork ;
    }
    if ( ifree == NULL )
    {
        if ( asaParm->PrintFinal || asaParm->PrintLevel >= 1 )
        {
            printf ("Insufficient memory for specified problem dimension %e\n",
                     (double) n) ;
        }
        status = 13 ;
        return (status) ;
    }

    mem = cgParm->memory ;/* cg_descent corresponds to mem = 0 */

    if ( (mem != 0) && (mem < 3) )
    {
        if ( asaParm->PrintFinal || asaParm->PrintLevel >= 1 )
        {
            printf ("memory = %i not permitted in cg_descent\n", mem) ;
            printf ("memory must be at least 3\n") ;
        }
        status = 12 ;
        return (status) ;
    }

    /* allocate work array */
    mem = MIN2 (mem, n) ;
    if ( Work == NULL )
    {
        if ( mem == 0 || asaParm->DynamicMemory)
        {
            work = (double *) malloc ((5*n+asaParm->m)*sizeof (double)) ;
        }
        else if ( cgParm->LBFGS || (mem >= n) )
        /* use L-BFGS to solve subproblem*/
        {
            i = 2*mem*(n+1) + 5*n + asaParm->m ;
            work = (double *) malloc (i*sizeof (double)) ;
        }
        else /* limited memory CG_DESCENT to solve subproblem */
        {
            i = (mem+7)*n + (3*mem+9)*mem + 5 + asaParm->m ;
            work = (double *) malloc (i*sizeof (double)) ;
        }
    }
    else work = Work ;
    if ( work == NULL )
    {
        if ( asaParm->PrintFinal || asaParm->PrintLevel >= 1 )
        {
            printf ("Insufficient memory for specified problem dimension %e\n",
                     (double) n) ;
        }
        status = 13 ;
        return (status) ;
    }

    d = Com.d = work ;
    g = Com.g = d+n ;
    xtemp = Com.xtemp = g+n ;
    gtemp = Com.gtemp = xtemp+n ;
    pg = Com.pg = gtemp+n ;
    Com.lastfvalues = pg+n ; /* size asaParm->m */
    Com.lo = lo ;
    Com.hi = hi ;
    Com.cbbiter = 0 ;
    Com.cbbfunc = 0 ;
    Com.cbbgrad = 0 ;
    Com.cgiter = 0 ;
    Com.cgfunc = 0 ;
    Com.cggrad = 0 ;

    ident = FALSE ;
    xnorm = ZERO ;
    for (j = 0; j < n; j++)
    {
        t = x [j] ;
        if      ( t > hi [j] ) t = hi [j] ;
        else if ( t < lo [j] ) t = lo [j] ;
        x [j] = t ;
        if ( xnorm < fabs (t) ) xnorm = fabs (t) ;
    }

    Com.f = asa_fg (g, x, &Com) ;

    pgnorm = ZERO ;
    gnorm = ZERO ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        gj = g [j] ;
        xg = xj - gj ;
        if      ( xg > hi [j] ) xp = hi [j] - xj ;
        else if ( xg < lo [j] ) xp = lo [j] - xj ;
        else                    xp = -gj ;
        pg [j] = xp ;
        pgnorm = MAX2 (pgnorm, fabs (xp)) ;
        gnorm = MAX2 (gnorm, fabs (gj)) ;
    }
    if ( asaParm->StopRule ) Com.tol = MAX2 (pgnorm*asaParm->StopFac, grad_tol) ;
    else                     Com.tol = grad_tol ;

    Com.pgnorm = Com.pgnorm_start = pgnorm ;
    if ( asa_tol (pgnorm, &Com) )
    {
        status = 0 ;
        goto Exit ;
    }

    if ( xnorm != ZERO ) Com.alpha = alpha = xnorm/gnorm ;
    else                 Com.alpha = alpha = ONE/gnorm ;

    /* compute gradient norm for inactive variables */
    ginorm = ZERO ;
    nfree = 0 ;
    gp = FALSE ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        tl = lo [j] ;
        th = hi [j] ;
        gj = g [j] ;
        xg = xj - alpha*gj ;
        if      ( (xg >= th) && (th-xj > pert_hi) ) gp = TRUE ;
        else if ( (xg <= tl) && (xj-tl > pert_lo) ) gp = TRUE ;
        if ( (xj-tl > pert_lo) && (th - xj > pert_hi) )
        {
            ginorm = MAX2 (ginorm, fabs (gj)) ;
            ifree [nfree] = j ;
            nfree++ ;
        }
    }
    Com.ginorm = ginorm ;
    Com.nfree = nfree ;

    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("\ninitial f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n",
                 Com.f, pgnorm, ginorm) ;
        printf ("            nfree = %i xnorm = %14.6e gp = %i\n",
                (int) nfree, xnorm, gp) ;
    }

    if ( (ginorm < Com.tau1*pgnorm) || gp || asaParm->GradProjOnly )
    {
        Com.cbbfunc = 1 ;
        Com.cbbgrad = 1 ;
        goto Grad_proj ;
    }
    else
    {
        Com.cgfunc = 1 ;
        Com.cggrad = 1 ;
        goto CG_descent ;
    }

    Grad_proj:
    if ( asaParm->PrintLevel >= 1 ) printf ("\nGradProj:\n") ;
    Com.DimReduce = FALSE ;
    status = asa_grad_proj(&Com) ;
    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("exit Grad_proj\n") ;
    }
    if ( Com.cbbiter >= cbb_totit ) status = 14 ;
    if ( status >= 0 ) goto Exit ;

    /* extract free variable */
    nfree = 0 ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
        {
            ifree [nfree] = j ;
            nfree++ ;
        }
    }
    Com.nfree = nfree ;

    CG_descent:
    if ( nfree != n )
    {
       asa_shrink_all (&Com) ;
       asa_copy (xtemp+nfree, x+nfree, n-nfree) ;
       Com.DimReduce = TRUE ;
    }
    else Com.DimReduce = FALSE ;

    if ( asaParm->PrintLevel >= 1 ) printf ("\nCG:\n") ;
    status = asa_descent (&Com) ;

    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("exit the CG subroutine\n") ;
    }
    if ( Com.DimReduce ) asa_expand_all (&Com) ;
    if ( Com.cgiter >= cg_totit ) status = 2 ;

    if ( status >= -2 ) goto Exit ;

    /* ginorm < tau2* pgnorm without hitting boundary */
    if ( status == -5 )
    {
        Com.alpha = asa_init_bbstep (&Com) ;
        goto Grad_proj ;

    }
    /* ginorm >= tau2* pgnorm and many components of x hit boundary  */
    else if ( status == -4 )
    {
        ginorm = ZERO ;
        nfree = 0 ;
        for (j = 0 ; j < n; j++)
        {
            xj = x [j] ;
            if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
            {
                t = fabs (g [j]) ;
                ginorm = MAX2 (ginorm, t) ;
                ifree [nfree] = j ;
                nfree++ ;
            }
        }
        Com.nfree = nfree ;
        Com.ginorm = ginorm ;

        if ( ginorm >= Com.tau1*Com.pgnorm ) goto CG_descent ;
        else
        {
           if ( asaParm->PrintLevel >= 1 ) printf ("ginorm < tau1* pgnorm\n") ;
           Com.alpha = asa_init_bbstep (&Com) ;
           goto Grad_proj ;
        }
    }
    /* ginorm >= tau2*pgnorm and only one component of x hits boundary */
    else if ( status == -3 )
    {
        if ( pgnorm < asaParm->pgdecay*MAX2 (Com.pgnorm_start, ONE) )
        {
            ident = asa_identify (x, g, Com.pgnorm, &Com) ;
        }
        if ( ident )
        {
            ident = FALSE ;
            ginorm = ZERO ;
            nfree = 0 ;
            for (j = 0 ; j < n; j++)
            {
                xj = x [j] ;
                if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                {
                    t = fabs (g [j]) ;
                    ginorm = MAX2 (ginorm, t) ;
                    ifree [nfree] = j ;
                    nfree++ ;
                }
            }
            Com.nfree = nfree ;
            Com.ginorm = ginorm ;
            if ( ginorm >= Com.tau1*Com.pgnorm ) goto CG_descent ;
            else
            {
               if ( asaParm->PrintLevel >= 1 )
                   printf ("ginorm < tau1* pgnorm\n" ) ;
               Com.alpha = asa_init_bbstep (&Com) ;
               goto Grad_proj ;
            }
        }
        else
        {
            Com.alpha = asa_init_bbstep (&Com) ;
            goto Grad_proj ;
        }
    }

    Exit:
    if ( (asaParm->PrintFinal) || (asaParm->PrintLevel >= 1) )
    {
        const char mess1 [] = "Possible causes of this error message:" ;
        const char mess2 [] = "   - your tolerance may be too strict: "
                              "grad_tol = " ;
        const char mess4 [] = "   - your gradient routine has an error" ;
        const char mess5 [] = "   - the parameter epsilon in "
                              "asa_descent_c.parm is too small" ;
        printf ("\nFinal convergence status = %d\n", status);
        if ( status == -2 )
        {
            printf ("Function value became nan at cg iteration %10.0e\n",
                     (double) Com.cgiter) ;
        }
        else if ( status == -1 )
        {
            printf ("Function value of starting point is nan at "
                     "cg iteration %10.0f\n", (double) Com.cgiter) ;
        }
        else if ( status == 0 )
        {
            printf ("Convergence tolerance for gradient satisfied\n") ;
        }
        else if ( status == 1 )
        {
            printf ("Terminating in cg since change in function value "
                    "<= feps*|f|\n") ;
        }
        else if ( status == 2 )
        {
            printf ("Number of iterations exceed specified limits "
                    "for cg routine\n") ;
            printf ("Iterations: %10.0f maxit: %10.0f totit: %10.0f\n",
                    (double) Com.cgiter, (double) Com.cgmaxit,
                    (double) cg_totit) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
        }
        else if ( status == 3 )
        {
            printf ("Slope always negative in cg line search\n") ;
            printf ("%s\n", mess1) ;
            printf ("   - your cost function has an error\n") ;
            printf ("%s\n", mess4) ;
        }
        else if ( status == 4 )
        {
            printf ("Line search fails in cg, too many secant steps\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
        }
        else if ( status == 5 )
        {
            printf ("Search direction not a descent direction in cg\n") ;
        }
        else if ( status == 6 )
        /* line search fails in cg, excessive eps updating */
        {
            printf ("Line search fails in cg iteration,\n") ;
            printf ("due to excessive updating of eps\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n\n", mess4) ;
        }

        else if ( status == 7 ) /* line search fails */
        {
            printf ("Line search fails in cg iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n%s\n\n", mess4, mess5) ;
        }
        else if ( status == 8 )
        {
            printf ("Debugger is on, function value does not improve in cg\n") ;
            printf ("new value: %25.16e old value: %25.16e\n",
                Com.f_debug, Com.f0) ;
        }
        else if ( status == 9 )
        {
            printf ("%i iterations without strict improvement in cost "
                    "or gradient in cg\n\n", Com.nslow) ;
        }
        else if ( status == 10 )
        {
            printf ("Insufficient memory for specified problem dimension %e"
                    " in cg\n", (double) Com.nfree) ;
        }
        else if ( status == 11 )
        {
            printf ("Function nan and could not be repaired in cg\n\n") ;
        }
        else if ( status == 12 )
        {
            printf ("memory = %i is an invalid choice for cg parameter\n",
                     cgParm->memory) ;
            printf ("memory should be either 0 or greater than 2\n\n") ;
        }
        else if ( status == 13 )
        {
            printf ("Insufficient memory\n") ;
        }
        else if ( status == 14 )
        {
            printf ("Number of iterations or function evaluation exceed\n"
                          "specified limits for cbb routine\n") ;
            printf ("Iterations: %e maxit: %e totit: %e\n",
                     (double) Com.cbbiter, (double) Com.pgmaxit,
                     (double) cbb_totit) ;
            printf ("Total function evaluations: %e maxfunc: %e\n",
                     (double) Com.nf, (double) Com.pgmaxfunc);
        }
        if ( status == 15 ) /* line search fails in cbb iteration */
        {
            printf ("Line search fails in cbb iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n", mess4) ;
        }

        if ( status == 16 )
        {
            printf ("Search direction not descent direction in "
                    "asa_grad_proj\n") ;
            printf ("directional derivative: %e\n", Com.gtd) ;
        }
        if ( status == 17 )
        {
             printf ("At cbb iteration %e function value became nan\n",
                      (double) Com.cbbiter) ;
        }

        printf ("projected gradient max norm: %13.6e\n", Com.pgnorm) ;
        printf ("function value:              %13.6e\n", Com.f) ;
        printf ("\nTotal cg  iterations:           %10.0f\n",
                (double) Com.cgiter) ;
        printf ("Total cg  function evaluations: %10.0f\n",
                (double) Com.cgfunc) ;
        printf ("Total cg  gradient evaluations: %10.0f\n",
                (double) Com.cggrad) ;
        printf ("Total cbb iterations:           %10.0f\n",
                (double) Com.cbbiter) ;
        printf ("Total cbb function evaluations: %10.0f\n",
                (double) Com.cbbfunc) ;
        printf ("Total cbb gradient evaluations: %10.0f\n",
                    (double) Com.cbbgrad) ;
        printf ("------------------------------------------\n") ;
        printf ("Total function evaluations:     %10.0f\n",
                (double) Com.nf) ;
        printf ("Total gradient evaluations:     %10.0f\n",
                (double) Com.ng) ;
        printf ("==========================================\n\n") ;
    }
    if ( iWork == NULL ) free (ifree) ;
    if ( Work == NULL ) free (work) ;
    if ( Stat != NULL )
    {
        Stat->f = Com.f ;
        Stat->pgnorm = Com.pgnorm ;
        Stat->cgiter = Com.cgiter ;
        Stat->cgfunc = Com.cgfunc ;
        Stat->cggrad = Com.cggrad ;
        Stat->cbbiter = Com.cbbiter ;
        Stat->cbbfunc = Com.cbbfunc ;
        Stat->cbbgrad = Com.cbbgrad ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_descent =========================================================
   =========================================================================
   cg_descent conjugate gradient algorithm with modifications to handle the
   bound constraints. Based on cg_descent Version 6.0.
   ========================================================================= */
 int asa_descent /*  return:
                      -5 (ginorm < tau2*pgnorm without hitting boundary)
                      -4 (ginorm >=tau2*pgnorm, many x components hit boundary)
                      -3 (ginorm >=tau2*pgnorm, one x component hits boundary)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f|)
                       2 (total iterations exceeded maxit)
                       3 (slope always negative in line search)
                       4 (number secant iterations exceed nsecant)
                       5 (search direction not a descent direction)
                       6 (excessive updating of eps)
                       7 (Wolfe conditions never satisfied)
                       8 (debugger is on and the function value increases)
                       9 (no cost or gradient improvement in
                          2n + Parm->nslow iterations)
                      10 (out of memory)
                      11 (function nan or +-INF and could not be repaired)
                      12 (invalid choice for memory parameter) */
(
    asa_com *Com
)
{
    ASA_INT     i, iter, IterRestart, maxit, n, n5, nfree, nf, ng, nrestart,
            nrestartsub ;
    int     nslow, slowlimit, IterQuad, status, PrintLevel, QuadF ;
    double  delta2, Qk, Ck, fbest, gbest, pgnorm, ginorm,
            f, ftemp, gnorm, xnorm, gnorm2, dnorm2, denom, bdist,
            t, dphi, dphi0, alpha, talpha,
            xj, gj, xg, xp, sts, sty, sk, xi,
            ykyk, ykgk, dkyk, beta, QuadTrust,
            *x, *d, *g, *xtemp, *gtemp, *lo, *hi, *pg, *work ;

    /* new variables added in Version 6.0 */
    int     l1, l2, j, k, mem, memsq, memk, memk_begin, mlast, mlast_sub,
            mp, mp_begin, mpp, nsub, spp, spp1, SkFstart, SkFlast,
            Subspace, UseMemory, Restart, LBFGS, InvariantSpace, IterSub,
            IterSubStart, IterSubRestart, FirstFull, SubSkip, SubCheck,
            StartSkip, StartCheck, DenseCol1, memk_is_mem, d0isg, qrestart ;
    double  gHg, scale, gsubnorm2,  ratio, stgkeep,
            alphaold, zeta, yty, ytg, t1, t2, t3, t4,
           *Rk, *Re, *Sk, *SkF, *stemp, *Yk, *SkYk,
           *dsub, *gsub, *gsubtemp, *gkeep, *tau, *vsub, *wsub ;

    asacg_parm *Parm ;
    asa_parm *asaParm ;

/* initialization */

    x = Com->x ;
    lo = Com->lo ;
    hi = Com->hi ;
    n = Com->n ;
    d = Com->d ;
    g = Com->g ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    pg = Com->pg ;
    nfree = Com->nfree ;
    nf = Com->nf ;
    ng = Com->ng ;
    pgnorm = Com->pgnorm ;
    ginorm = Com->ginorm ;
    Parm = Com->cgParm ;
    asaParm = Com->asaParm ;
    PrintLevel = Parm->PrintLevel ;

    if ( PrintLevel >= 1 )
    {
        printf ("Dimension in CG, nfree = %i\n", (int) nfree) ;
    }

    mem = MIN2 (nfree, Parm->memory) ;
    if ( asaParm->DynamicMemory )      /*dynamically assign memory */
    {
        if ( Parm->LBFGS || (mem >= nfree) ) /* use L-BFGS */
        {
            work = (double *) malloc (2*mem*(nfree+1)*sizeof (double)) ;
        }
        else /* limited memory CG_DESCENT */
        {
            i = (mem+2)*nfree + (3*mem+9)*mem + 5 ;
            work = (double *) malloc (i*sizeof (double)) ;
        }
        if ( work == NULL )
        {
            status = 10 ;
            goto Exit ;
        }
    }
    else
    {
        work = Com->lastfvalues + asaParm->m ;
    }

    qrestart = MIN2 (nfree, Parm->qrestart) ;
    Com->Wolfe = FALSE ; /* initially Wolfe line search not performed */
    QuadF = FALSE ;     /* initially function assumed to be nonquadratic */
    LBFGS = FALSE ;
    UseMemory = FALSE ;/* do not use memory */
    Subspace = FALSE ; /* full space, check subspace condition if UseMemory */
    FirstFull = FALSE ;/* not first full iteration after leaving subspace */
    memk = 0 ;         /* number of vectors in current memory */

    /* the conjugate gradient algorithm is restarted every nrestart iteration */
    nrestart = (ASA_INT) (((double) nfree)*Parm->restart_fac) ;

    /* abort when number of iterations reaches maxit in one pass through cg */
    Com->cgmaxit = maxit = Parm->maxit ;

    if ( mem > 0 )
    {
        if ( (mem == nfree) || Parm->LBFGS )
        {
            LBFGS = TRUE ;      /* use L-BFGS */
            mlast = -1 ;
            Sk = work ;
            Yk = Sk + mem*nfree ;
            SkYk = Yk + mem*nfree ;
            tau = SkYk + mem ;
        }
        else
        {
            UseMemory = TRUE ; /* previous search direction will be saved */
            SubSkip = 0 ;      /* number of iterations to skip checking memory*/
            SubCheck = mem*Parm->SubCheck ; /* number of iterations to check */
            StartCheck = 0 ;   /* start checking memory at iteration 0 */
            InvariantSpace = FALSE ; /* iterations not in invariant space */
            FirstFull = TRUE ;       /* first iteration in full space */
            nsub = 0 ;               /* initial subspace dimension */
            memsq = mem*mem ;
            SkF = work ;    /* directions in memory (x_k+1 - x_k) */
            stemp = SkF + mem*nfree ;/* stores x_k+1 - x_k */
            gkeep = stemp + nfree ;  /* store grad when first direction != -g */
            Sk = gkeep + nfree ;   /* Sk = Rk at start of LBFGS in subspace */
            Rk = Sk + memsq ;  /* upper triangular factor in SkF = Zk*Rk */
            /* zero out Rk to ensure lower triangle is 0 */
            asa_init (Rk, ZERO, memsq) ;
            Re = Rk + memsq ;  /* end column of Rk, used for new direction */
            Yk = Re + mem+1 ;
            SkYk = Yk + memsq+mem+2 ; /* dot products sk'yk in the subspace */
            tau = SkYk + mem ;       /* stores alpha in Nocedal and Wright */
            dsub = tau + mem ;       /* direction projection in subspace */
            gsub = dsub + mem ;      /* gradient projection in subspace */
            gsubtemp = gsub + mem+1 ;/* new gsub before update */
            wsub = gsubtemp + mem ;  /* mem+1 work array for triangular solve */
            vsub = wsub + mem+1 ;    /* mem work array for triangular solve */
        }
    }

    fbest = ASA_INF ;
    gbest = ASA_INF ;
    nslow = 0 ;
    slowlimit = 2*nfree + Parm->nslow ;
    n5 = nfree % 5 ;

    Ck = ZERO ;
    Qk = ZERO ;

    f = Com->f ;
    Com->f0 = f + f ;
    Com->SmallCost = fabs (f)*Parm->SmallCost ;

    /* compute inf-norm of x and distance to bounds */
    xnorm = ZERO ;
    bdist = ASA_INF ;
    for (i = 0; i < nfree; i++)
    {
        xi = x [i] ;
        if ( xnorm < fabs (xi) ) xnorm = fabs (xi) ;
        t = xi - lo [i] ;
        if ( bdist > t ) bdist = t ;
        t = hi [i] - xi ;
        if ( bdist > t ) bdist = t ;
    }
    Com->bdist = bdist ;

    /* set d = -g, compute gnorm  = infinity norm of g and
                           gnorm2 = square of 2-norm of g */
    gnorm = asa_update_inf2 (g, g, d, &gnorm2, nfree) ;
    dnorm2 = gnorm2 ;

    if ( PrintLevel >= 2 )
    {
        printf ("iter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e "
                "memk: %i\n", (int) 0, f, pgnorm, ginorm, memk) ;
    }

    dphi0 = -gnorm2 ;
    delta2 = 2*Parm->delta - ONE ;
    alpha = Parm->step ;
    if ( alpha == ZERO )
    {
        if ( xnorm == ZERO )
        {
            if ( f != ZERO ) alpha = 2.*fabs (f)/gnorm2 ;
            else             alpha = ONE ;
        }
        else    alpha = Parm->psi0*xnorm/gnorm ;
    }

    Com->df0 = -2.0*fabs(f)/alpha ;

    Restart = FALSE ;  /* do not restart the algorithm */
    IterRestart = 0 ;  /* counts number of iterations since last restart */
    IterSub = 0 ;      /* counts number of iterations in subspace */
    IterQuad = 0 ;     /* counts number of iterations that function change
                          is close to that of a quadratic */

/*  start the conjugate gradient iteration
    alpha starts as old step, ends as final step for current iteration
    f is function value for alpha = 0
    Com->QuadOK = TRUE means that a quadratic step was taken */

    for (iter = 1; iter <= maxit; iter++)
    {
        /* save old alpha to simplify formula computing subspace direction */
        alphaold = alpha ;
        Com->QuadOK = FALSE ;
        alpha = Parm->psi2*alpha ;
        if ( Com->bdist > 0 )
        {
            Com->minflag = FALSE ;
            Com->minstep = Com->bdist/sqrt(dnorm2) ;
        }
        else asa_maxstep (x, d, Com) ;
        if ( f != ZERO ) t = fabs ((f-Com->f0)/f) ;
        else             t = ONE ;
        Com->UseCubic = TRUE ;
        if ( (t < Parm->CubicCutOff) || !Parm->UseCubic ) Com->UseCubic = FALSE;
        if ( Parm->QuadStep )
        {
            /* test if quadratic interpolation step should be tried */
            if ( ((t > Parm->QuadCutOff)&&(fabs(f) >= Com->SmallCost)) || QuadF)
            {
                if ( QuadF )
                {
                    if ( asaParm->HardConstraint )
                    {
                        talpha = Parm->psi1*alpha ;
                        t = asaParm->parm1*talpha ;
                        if ( Com->minstep < t )
                        {
                            if ( !Com->minflag ) asa_maxstep (x, d, Com) ;
                        }
                        if ( Com->minstep >= t )
                        {
                            t = MIN2 (talpha, asaParm->parm2*Com->minstep) ;
                            if ( t < talpha )
                            {
                                if ( !Com->minflag )
                                {
                                    asa_maxstep (x, d, Com) ;
                                    talpha = MIN2 (talpha,
                                                  asaParm->parm2*Com->minstep) ;
                                }
                                talpha = t ;
                            }
                            Com->alpha = talpha ;
                            status = asacg_evaluate ("g", "y", Com) ;
                            if ( status ) goto Exit ;
                            if ( Com->df > dphi0 )
                            {
                                alpha = -dphi0/((Com->df-dphi0)/Com->alpha) ;
                                Com->QuadOK = TRUE ;
                            }
                            else if ( LBFGS )
                            {
                                if ( memk >= nfree )
                                {
                                    alpha = ONE ;
                                    Com->QuadOK = TRUE ;
                                }
                                else alpha = 2. ;
                            }
                            else if ( Subspace )
                            {
                                if ( memk >= nsub )
                                {
                                    alpha = ONE ;
                                    Com->QuadOK = TRUE ;
                                }
                                else alpha = 2. ;
                            }
                        }
                        else if ( LBFGS )
                        {
                            if ( memk >= nfree )
                            {
                                alpha = ONE ;
                                Com->QuadOK = TRUE ;
                            }
                            else alpha = 2. ;
                        }
                        else if ( Subspace )
                        {
                            if ( memk >= nsub )
                            {
                                alpha = ONE ;
                                Com->QuadOK = TRUE ;
                            }
                            else alpha = 2. ;
                        }
                    }
                    else
                    {
                        Com->alpha = Parm->psi1*alpha ;
                        status = asacg_evaluate ("g", "y", Com) ;
                        if ( status ) goto Exit ;
                        if ( Com->df > dphi0 )
                        {
                            alpha = -dphi0/((Com->df-dphi0)/Com->alpha) ;
                            Com->QuadOK = TRUE ;
                        }
                        else if ( LBFGS )
                        {
                            if ( memk >= nfree )
                            {
                                alpha = ONE ;
                                Com->QuadOK = TRUE ;
                            }
                            else  alpha = 2. ;
                        }
                        else if ( Subspace )
                        {
                            if ( memk >= nsub )
                            {
                                alpha = ONE ;
                                Com->QuadOK = TRUE ;
                            }
                            else  alpha = 2. ;
                        }
                    }

                }
                else
                {
                    if ( asaParm->HardConstraint )
                    {
                        t = MAX2 (Parm->psi_lo, Com->df0/(dphi0*Parm->psi2)) ;
                        talpha = MIN2 (t, Parm->psi_hi)*alpha ;
                        t = asaParm->parm1*talpha ;
                        if ( Com->minstep < t )
                        {
                            if ( !Com->minflag ) asa_maxstep (x, d, Com) ;
                        }
                        if ( Com->minstep >= t )
                        {
                            t = MIN2 (talpha, asaParm->parm2*Com->minstep) ;
                            if ( t < talpha )
                            {
                                if ( !Com->minflag )
                                {
                                    asa_maxstep (x, d, Com) ;
                                    talpha = MIN2 (talpha,
                                                  asaParm->parm2*Com->minstep) ;
                                }
                                talpha = t ;
                            }
                            Com->alpha = talpha ;
                            status = asacg_evaluate ("f", "y", Com) ;
                            if ( status ) goto Exit ;
                            ftemp = Com->f ;
                            denom = 2.*(((ftemp-f)/Com->alpha)-dphi0) ;
                            if ( denom > ZERO )
                            {
                                t = -dphi0*Com->alpha/denom ;
                                /* safeguard */
                                if ( ftemp >= f )
                                  alpha = MAX2 (t, Com->alpha*Parm->QuadSafe) ;
                                else  alpha = t ;
                                Com->QuadOK = TRUE ;
                            }
                        }
                    }
                    else
                    {
                        t = MAX2 (Parm->psi_lo, Com->df0/(dphi0*Parm->psi2)) ;
                        Com->alpha = MIN2 (t, Parm->psi_hi)*alpha ;
                        status = asacg_evaluate ("f", "y", Com) ;
                        if ( status ) goto Exit ;
                        ftemp = Com->f ;
                        denom = 2.*(((ftemp-f)/Com->alpha)-dphi0) ;
                        if ( denom > ZERO )
                        {
                            t = -dphi0*Com->alpha/denom ;
                            /* safeguard */
                            if ( ftemp >= f )
                              alpha = MAX2 (t, Com->alpha*Parm->QuadSafe) ;
                            else  alpha = t ;
                            Com->QuadOK = TRUE ;
                        }
                    }
                }
            }
        }
        Com->f0 = f ;                          /* f0 saved as prior value */
        Com->df0 = dphi0 ;

        if ( PrintLevel >= 3 )
        {
            if ( Com->minflag )
            {
                printf ("minstep =%14.6e, maxstep =%14.6e\n",
                         Com->minstep, Com->maxstep) ;
            }
            else printf ("bdist  =%14.6e\n", Com->bdist) ;
            printf ("QuadOK: %2i initial a: %14.6e f0: %14.6e dphi0: %14.6e\n",
                    Com->QuadOK, alpha, Com->f0, dphi0) ;
            if ( (alpha > Com->minstep) && Com->QuadOK )
            {
                printf("Quadratic step > minstep to boundary\n") ;
            }
        }

        /* parameters in Wolfe, approximate Wolfe conditions, and in update */
        Qk = Parm->Qdecay*Qk + ONE ;
        Ck = Ck + (fabs (f) - Ck)/Qk ;        /* average cost magnitude */

        if ( Parm->PertRule ) Com->fpert = f + Com->eps*fabs (f) ;
        else                  Com->fpert = f + Com->eps ;

        Com->wolfe_hi = Parm->delta*dphi0 ;
        Com->wolfe_lo = Parm->sigma*dphi0 ;
        Com->awolfe_hi = delta2*dphi0 ;
        Com->alpha = alpha ;/* either double prior step or quadratic fit step */

        /* perform line search */
        status = asa_line (Com) ;

        /*try approximate Wolfe line search if ordinary Wolfe fails */
        if ( (status > 0) && !Com->AWolfe )
        {
            if ( PrintLevel >= 2 )
            {
                 printf ("\nWOLFE LINE SEARCH FAILS\n") ;
            }
            if ( status != 3 )
            {
                Com->AWolfe = TRUE ;
                status = asa_line (Com) ;
            }
        }

        alpha = Com->alpha ;
        f = Com->f ;
        dphi = Com->df ;

        if ( (status > 0) || (status == -1) || (status == -2) ) goto Exit ;

        /*Test for convergence to within machine epsilon
          [set feps to zero to remove this test] */

        if ( (-alpha*dphi0 <= Parm->feps*fabs (f)) && (status == 0) )
        {
            status = 1 ;
            goto Exit ;
        }

        /* test how close the cost function changes are to that of a quadratic
           QuadTrust = 0 means the function change matches that of a quadratic*/
        t = alpha*(dphi+dphi0) ;
        if ( fabs (t) <= Parm->qeps*MIN2 (Ck, ONE) ) QuadTrust = ZERO ;
        else QuadTrust = fabs((2.0*(f-Com->f0)/t)-ONE) ;
        if ( QuadTrust <= Parm->qrule) IterQuad++ ;
        else                           IterQuad = 0 ;

        if ( IterQuad == qrestart ) QuadF = TRUE ;
        IterRestart++ ;

        if ( !Com->AWolfe )
        {
            if ( fabs (f-Com->f0) < Parm->AWolfeFac*Ck )
            {
                Com->AWolfe = TRUE ;
                if ( Com->Wolfe ) Restart = TRUE ;
            }
        }

        if ( (mem > 0) && !LBFGS )
        {
            if ( UseMemory )
            {
                if ( (iter - StartCheck > SubCheck) && !Subspace )
                {
                    StartSkip = iter ;
                    UseMemory = FALSE ;
                    if ( SubSkip == 0 ) SubSkip = mem*Parm->SubSkip ;
                    else                SubSkip *= 2 ;
                    if ( PrintLevel >= 2 )
                    {
                        printf ("skip subspace %i iterations\n", SubSkip) ;
                    }
                }
            }
            else
            {
                if ( iter - StartSkip > SubSkip )
                {
                    StartCheck = iter ;
                    UseMemory = TRUE ;
                    memk = 0 ;
                }
            }
        }

        if ( !UseMemory )
        {
            if ( !LBFGS )
            {
                if ( (IterRestart >= nrestart) || ((IterQuad == qrestart)
                     && (IterQuad != IterRestart)) ) Restart = TRUE ;
            }
        }
        else
        {
            if ( Subspace ) /* the iteration is in the subspace */
            {
                IterSubRestart++ ;

                /* compute project of g into subspace */
                gsubnorm2 = ZERO ;
                mp = SkFstart ;
                j = nsub - mp ;

                /* multiply basis vectors by new gradient */
                asa_matvec (wsub, SkF, gtemp, nsub, nfree, 0) ;

                /* rearrange wsub and store in gsubtemp
                   (elements associated with old vectors should
                    precede elements associated with newer vectors */
                asa_copy0 (gsubtemp, wsub+mp, j) ;
                asa_copy0 (gsubtemp+j, wsub, mp) ;

                /* solve Rk'y = gsubtemp */
                asa_trisolve (gsubtemp, Rk, mem, nsub, 0) ;
                gsubnorm2 = asa_dot0 (gsubtemp, gsubtemp, nsub) ;
                gnorm2 = asa_dot (gtemp, gtemp, nfree);
                ratio = sqrt(gsubnorm2/gnorm2) ;
                if ( ratio < ONE - Parm->eta1  ) /* Exit Subspace */
                {
                   if ( PrintLevel >= 2 )
                   {
                       printf ("iter: %i exit subspace\n", (int) iter) ;
                   }
                   FirstFull = TRUE ; /* first iteration in full space */
                   Subspace = FALSE ; /* leave the subspace */
                   InvariantSpace = FALSE ;
                   /* check the subspace condition for SubCheck iterations
                      starting from the current iteration (StartCheck) */
                   StartCheck = iter ;
                }
                else
                {
                   /* Check if a restart should be done in subspace */
                   if ( IterSubRestart == nrestartsub ) Restart = TRUE ;
                }
            }
            else  /* in full space */
            {
                if ( (IterRestart == 1) || FirstFull ) memk = 0 ;
                if ( (memk == 1) && InvariantSpace )
                {
                     memk = 0 ;
                     InvariantSpace = FALSE ;
                }
                if (memk < mem )
                {
                    memk_is_mem = FALSE ;
                    SkFstart = 0 ;
                    /* SkF stores basis vector of the form alpha*d
                       We factor SkF = Zk*Rk where Zk has orthonormal columns
                       and Rk is upper triangular. Zk is not stored; wherever
                       it is needed, we use SkF * inv (Rk) */
                    if (memk == 0)
                    {
                        mlast = 0 ;  /* starting pointer in the memory */
                        memk = 1 ;   /* dimension of current subspace */

                        t = sqrt(dnorm2) ;
                        zeta = alpha*t ;
                        Rk [0] = zeta ;
                        asa_scale (SkF, d, alpha, nfree) ;
                        Yk [0] = (dphi - dphi0)/t ;
                        gsub [0] = dphi/t ;
                        SkYk [0] = alpha*(dphi-dphi0) ;
                        FirstFull = FALSE ;
                        if ( IterRestart > 1 )
                        {
                           /* Need to save g for later correction of first
                              column of Yk. Since g does not lie in the
                              subspace and the first column is dense */
                           asa_copy (gkeep, g, nfree) ;
                           /* Also store dot product of g with the first
                              direction vector -- this saves a later dot
                              product when we fix the first column of Yk */
                           stgkeep = dphi0*alpha ;
                           d0isg = FALSE ;
                        }
                        else d0isg = TRUE ;
                    }
                    else
                    {
                        mlast = memk ; /* starting pointer in the memory */
                        memk++ ;       /* total number of Rk in the memory */
                        mpp = mlast*nfree ;
                        spp = mlast*mem ;
                        asa_scale (SkF+mpp, d, alpha, nfree) ;

                        /* check if the alphas are far from 1 */
                        if ((fabs(alpha-5.05)>4.95)||(fabs(alphaold-5.05)>4.95))
                        {
                            /* multiply basis vectors by new direction vector */
                            asa_matvec (Rk+spp, SkF, SkF+mpp, mlast, nfree, 0) ;

                            /* solve Rk'y = wsub to obtain the components of the
                               new direction vector relative to the orthonormal
                               basis Z in S = ZR, store in next column of Rk */
                            asa_trisolve (Rk+spp, Rk, mem, mlast, 0) ;
                        }
                        else /* alphas are close to 1 */
                        {
                            t1 = -alpha ;
                            t2 = beta*alpha/alphaold ;
                            for (j = 0; j < mlast; j++)
                            {
                                Rk [spp+j] = t1*gsub [j] + t2*Rk [spp-mem+j] ;
                            }
                        }
                        t = alpha*alpha*dnorm2 ;
                        t1 = asa_dot0 (Rk+spp, Rk+spp, mlast) ;
                        if (t <= t1)   zeta = t*Parm->eta2 ;
                        else           zeta = sqrt(t-t1);

                        Rk [spp+mlast] = zeta ;
                        t = - zeta/alpha ; /* t = asa_dot0 (Zk+mlast*n, g, n)*/
                        Yk [spp-mem+mlast] = t ;
                        gsub [mlast] = t ;

                        /* multiply basis vectors by new gradient */
                        asa_matvec (wsub, SkF, gtemp, mlast, nfree, 0) ;
                        /* exploit dphi for last multiply */
                        wsub [mlast] = alpha*dphi ;
                        /* solve for new gsub */
                        asa_trisolve (wsub, Rk, mem, memk, 0) ;
                        /* subtract old gsub from new gsub = column of Yk */
                        asacg_Yk (Yk+spp, gsub, wsub, NULL, memk) ;

                        SkYk [mlast] = alpha*(dphi-dphi0) ;
                    }
                }
                else  /* memk = mem */
                {
                    memk_is_mem = TRUE ;
                    mlast = mem-1 ;
                    asa_scale (stemp, d, alpha, nfree) ;
                    /* compute projection of s_k = alpha_k d_k into subspace
                       check if the alphas are far from 1 */
                    if ((fabs(alpha-5.05)>4.95)||(fabs(alphaold-5.05)>4.95))
                    {
                        mp = SkFstart ;
                        j = mem - mp ;

                        /* multiply basis vectors by sk */
                        asa_matvec (wsub, SkF, stemp, mem, nfree, 0) ;
                        /* rearrange wsub and store in Re = end col Rk */
                        asa_copy0 (Re, wsub+mp, j) ;
                        asa_copy0 (Re+j, wsub, mp) ;

                        /* solve Rk'y = Re */
                        asa_trisolve (Re, Rk, mem, mem, 0) ;
                    }
                    else /* alphas close to 1 */
                    {
                        t1 = -alpha ;
                        t2 = beta*alpha/alphaold ;
                        for (j = 0; j < mem; j++)
                        {
                            Re [j] = t1*gsub [j] + t2*Re [j-mem] ;
                        }
                    }

                    /* t = 2-norm squared of s_k */
                    t = alpha*alpha*dnorm2 ;
                    /* t1 = 2-norm squared of projection */
                    t1 = asa_dot0 (Re, Re, mem) ;
                    if (t <= t1)   zeta = t*Parm->eta2 ;
                    else           zeta = sqrt(t-t1);

                    /* dist from new search direction to prior subspace*/
                    Re [mem] = zeta ;

                    /* projection of prior g on new orthogonal
                       subspace vector */
                    t = -zeta/alpha ; /* t = asa_dot(Zk+mpp, g, n)*/
                    gsub [mem] = t ;
                    Yk [memsq] = t ;  /* also store it in Yk */

                    spp = memsq + 1 ;
                    mp = SkFstart ;
                    j = mem - mp ;

                    /* multiply basis vectors by gtemp */
                    asa_matvec (vsub, SkF, gtemp, mem, nfree, 0) ;

                    /* rearrange and store in wsub */
                    asa_copy0 (wsub, vsub+mp, j) ;
                    asa_copy0 (wsub+j, vsub, mp) ;

                    /* solve Rk'y = wsub */
                    asa_trisolve (wsub, Rk, mem, mem, 0) ;
                    wsub [mem] = (alpha*dphi - asa_dot0 (wsub, Re, mem))/zeta;

                    /* add new column to Yk, store new gsub */
                    asacg_Yk (Yk+spp, gsub, wsub, NULL, mem+1) ;

                    /* store sk (stemp) at SkF+SkFstart */
                    asa_copy (SkF+SkFstart*nfree, stemp, nfree) ;
                    SkFstart++ ;
                    if ( SkFstart == mem ) SkFstart = 0 ;

                    mp = SkFstart ;
                    for (k = 0; k < mem; k++)
                    {
                        spp = (k+1)*mem + k ;
                        t1 = Rk [spp] ;
                        t2 = Rk [spp+1] ;
                        t = sqrt(t1*t1 + t2*t2) ;
                        t1 = t1/t ;
                        t2 = t2/t ;

                        /* update Rk */
                        Rk [k*mem+k] = t ;
                        for (j = (k+2); j <= mem; j++)
                        {
                            spp1 = spp ;
                            spp = j*mem + k ;
                            t3 = Rk [spp] ;
                            t4 = Rk [spp+1] ;
                            Rk [spp1] = t1*t3 + t2*t4 ;
                            Rk [spp+1] = t1*t4 - t2*t3 ;
                        }
                        /* update Yk */
                        if ( k < 2 ) /* mem should be greater than 2 */
                        {
                            /* first 2 rows are dense */
                            spp = k ;
                            for (j = 1; j < mem; j++)
                            {
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            spp1 = spp ;
                            spp = mem*mem + 1 + k ;
                            t3 = Yk [spp] ;
                            t4 = Yk [spp+1] ;
                            Yk [spp1] = t1*t3 + t2*t4 ;
                            Yk [spp+1] = t1*t4 -t2*t3 ;
                        }
                        else if ( (k == 2) && (2 < mem-1))
                        {
                            spp = k ;

                            /* col 1 dense since the oldest direction
                               vector has been dropped */
                            j = 1 ;
                            spp1 = spp ;
                            spp = j*mem + k ;
                            /* single nonzero percolates down the column */
                            t3 = Yk [spp] ;  /* t4 = 0. */
                            Yk [spp1] = t1*t3 ;
                            Yk [spp+1] = -t2*t3 ;
                            /* process rows in Hessenberg part of matrix */
                            for (j = 2; j < mem; j++)
                            {
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            spp1 = spp ;
                            spp = mem*mem + 1 + k ;
                            t3 = Yk [spp] ;
                            t4 = Yk [spp+1] ;
                            Yk [spp1] = t1*t3 + t2*t4 ;
                            Yk [spp+1] = t1*t4 -t2*t3 ;
                        }
                        else if ( k < (mem-1) )
                        {
                            spp = k ;

                            /* process first column */
                            j = 1 ;
                            spp1 = spp ;
                            spp = j*mem + k ;
                            t3 = Yk [spp] ;  /* t4 = 0. */
                            Yk [spp1] = t1*t3 ;
                            Yk [spp+1] = -t2*t3 ;

                            /* process rows in Hessenberg part of matrix */
                            j = k-1 ;
                            spp = (j-1)*mem+k ;
                            spp1 = spp ;
                            spp = j*mem + k ;
                            t3 = Yk [spp] ;
                            Yk [spp1] = t1*t3 ; /* t4 = 0. */
                            /* Yk [spp+1] = -t2*t3 ;*/
                            /* Theoretically this element is zero */
                            for (j = k; j < mem; j++)
                            {
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            spp1 = spp ;
                            spp = mem*mem + 1 + k ;
                            t3 = Yk [spp] ;
                            t4 = Yk [spp+1] ;
                            Yk [spp1] = t1*t3 + t2*t4 ;
                            Yk [spp+1] = t1*t4 -t2*t3 ;
                        }
                        else /* k = mem-1 */
                        {
                            spp = k ;

                            /* process first column */
                            j = 1 ;
                            spp1 = spp ;
                            spp = j*mem + k ;
                            t3 = Yk [spp] ; /* t4 = 0. */
                            Yk [spp1] = t1*t3 ;

                            /* process rows in Hessenberg part of matrix */
                            j = k-1 ;
                            spp = (j-1)*mem+k ;
                            spp1 = spp ;
                            spp = j*mem + k ;
                            t3 = Yk [spp] ; /* t4 = 0. */
                            Yk [spp1] = t1*t3 ;

                            j = k ;
                            spp1 = spp ;
                            spp = j*mem + k ; /* j=mem-1 */
                            t3 = Yk [spp] ;
                            t4 = Yk [spp+1] ;
                            Yk [spp1] = t1*t3 + t2*t4 ;

                            spp1 = spp ;
                            spp = mem*mem + 1 + k ; /* j=mem */
                            t3 = Yk [spp] ;
                            t4 = Yk [spp+1] ;
                            Yk [spp1] = t1*t3 + t2*t4 ;
                        }
                        /* update g in subspace */
                        if ( k < (mem-1) )
                        {
                            t3 = gsub [k] ;
                            t4 = gsub [k+1] ;
                            gsub [k] = t1*t3 + t2*t4 ;
                            gsub [k+1] = t1*t4 -t2*t3 ;
                        }
                        else /* k = mem-1 */
                        {
                            t3 = gsub [k] ;
                            t4 = gsub [k+1] ;
                            gsub [k] = t1*t3 + t2*t4 ;
                        }
                    }

                    /* update SkYk */
                    for (k = 0; k < mlast; k++) SkYk [k] = SkYk [k+1] ;
                    SkYk [mlast] = alpha*(dphi-dphi0) ;
                }

                /* calculate t = ||gsub|| / ||gtemp||  */
                gsubnorm2 = asa_dot0 (gsub, gsub, memk) ;
                gnorm2 = asa_dot (gtemp, gtemp, nfree) ;
                ratio = sqrt ( gsubnorm2/gnorm2 ) ;
                if ( ratio > ONE-Parm->eta2) InvariantSpace = TRUE ;

                /* check to see whether to enter subspace */
                if ( ((memk > 1) && InvariantSpace) ||
                     ((memk == mem) && (ratio > ONE-Parm->eta0)) )
                {
                    if ( PrintLevel >= 2 )
                    {
                        if ( InvariantSpace )
                        {
                            printf ("iter: %i invariant space, "
                                    "enter subspace\n", (int) iter) ;
                        }
                        else
                        {
                            printf ("iter: %i enter subspace\n", (int) iter) ;
                        }
                    }
                    /* if the first column is dense, we need to correct it
                       now since we do not know the entries until the basis
                       is determined */
                    if ( !d0isg && !memk_is_mem )
                    {
                        wsub [0] = stgkeep ;
                        /* mlast = memk -1 */
                        asa_matvec (wsub+1, SkF+nfree, gkeep, mlast, nfree, 0) ;
                        /* solve Rk'y = wsub */
                        asa_trisolve (wsub, Rk, mem, memk, 0) ;
                        /* corrected first column of Yk */
                        Yk [1] -= wsub [1] ;
                        asa_scale0 (Yk+2, wsub+2, -ONE, memk-2) ;
                    }
                    if ( d0isg && !memk_is_mem ) DenseCol1 = FALSE ;
                    else                         DenseCol1 = TRUE ;

                    Subspace = TRUE ;
                    /* reset subspace skipping to 0, need to test invariance */
                    SubSkip = 0 ;
                    IterSubRestart = 0 ;
                    IterSubStart = IterSub ;
                    nsub = memk ; /* dimension of subspace */
                    nrestartsub = (int) (((double) nsub)*Parm->restart_fac) ;
                    mp_begin = mlast ;
                    memk_begin = nsub ;
                    SkFlast = (SkFstart+nsub-1) % mem ;
                    asa_copy0 (gsubtemp, gsub, nsub) ;
                    /* Rk contains the sk for subspace, initialize Sk = Rk */
                    asa_copy (Sk, Rk, (int) mem*nsub) ;
                }
                else
                {
                   if ( (IterRestart == nrestart) ||
                       ((IterQuad == qrestart) && (IterQuad != IterRestart)) )
                   {
                       Restart = TRUE ;
                   }
                }
            } /* done checking the full space */
        } /* done using the memory */

            Com->bdist -= alpha*sqrt(dnorm2) ;
            /* compute the ginorm and pgnorm of gtemp */
            pgnorm = asa_update_pg (xtemp, gtemp, &ginorm, Com->bdist,
                                    lo, hi, nfree, n, Com) ;

            if ( asa_tol (pgnorm, Com) )
            {
                status = 0 ;
                for (j = 0; j < nfree; j++)
                {
                    xj = xtemp [j] ;
                    x [j] = xj ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                for (; j < n; j++)
                {
                    xj = x [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                goto Exit1 ;
            }
            if ( ginorm < pgnorm*Com->tau2 ) status = -5 ;
            if ( status < -2 )
            {
                sts = ZERO ;
                sty = ZERO ;
                for (j = 0; j < nfree; j++)
                {
                    t = xtemp[j] ;
                    sk = t - x [j] ;
                    x [j] = t ;
                    sts += sk*sk ;

                    t = gtemp [j] ;
                    sty += sk*(t-g [j]) ;
                    g [j] = t ;
                }
                Com->sts = sts ;
                Com->sty = sty ;
                goto Exit ;
            }

        /* compute search direction */
        if ( LBFGS )
        {
            if ( IterRestart == nrestart ) /* restart the l-bfgs method */
            {
                IterRestart = 0 ;
                IterQuad = 0 ;
                mlast = -1 ;
                memk = 0 ;

                /* copy xtemp to x */
                asa_copy (x, xtemp, nfree) ;

                /* set g = gtemp, d = -g, compute 2-norm of g */
                gnorm2 = asa_update_2 (g, gtemp, d, nfree) ;

                dnorm2 = gnorm2 ;
                dphi0 = -gnorm2 ;
            }
            else
            {
                mlast = (mlast+1) % mem ;
                spp = mlast*nfree ;
                asa_step (Sk+spp, xtemp, x, -ONE, nfree) ;
                asa_step (Yk+spp, gtemp, g, -ONE, nfree) ;
                SkYk [mlast] = alpha*(dphi-dphi0) ;
                if (memk < mem) memk++ ;

                /* copy xtemp to x */
                asa_copy (x, xtemp, nfree) ;

                /* copy gtemp to g and compute 2-norm of g */
                gnorm2 = asa_update_2 (g, gtemp, NULL, nfree) ;

                /* calculate Hg = H g, saved in gtemp */
                mp = mlast ;  /* memk is the number of vectors in the memory */
                for (j = 0; j < memk; j++)
                {
                    mpp = mp*nfree ;
                    t = asa_dot (Sk+mpp, gtemp, nfree)/SkYk[mp] ;
                    tau [mp] = t ;
                    asa_daxpy (gtemp, Yk+mpp, -t, nfree) ;
                    mp -=  1;
                    if ( mp < 0 ) mp = mem-1 ;
                }
                /* scale = (alpha*dnorm2)/(dphi-dphi0) ; */
                scale = SkYk[mlast]/asa_dot (Yk+mlast*nfree, Yk+mlast*nfree,
                                             nfree) ;

                asa_scale (gtemp, gtemp, scale, nfree) ;

                for (j = 0; j < memk; j++)
                {
                    mp +=  1;
                    if ( mp == mem ) mp = 0 ;
                    mpp = mp*nfree ;
                    t = asa_dot (Yk+mpp, gtemp, nfree)/SkYk[mp] ;
                    asa_daxpy (gtemp, Sk+mpp, tau [mp]-t, nfree) ;
                }

                /* set d = -gtemp, compute 2-norm of gtemp */
                dnorm2 = asa_update_2 (NULL, gtemp, d, nfree) ;
                dphi0 = -asa_dot (g, gtemp, nfree) ;
            }
        } /* end of LBFGS */

        else if ( Subspace ) /* compute search direction in subspace */
        {
            IterSub++ ;

            /* set x = xtemp, g = gtemp */
            asa_update_xy (x, xtemp, g, gtemp, nfree) ;

            if ( Restart ) /*restart in subspace*/
            {
                Restart = FALSE ;
                IterRestart = 0 ;
                IterSubRestart = 0 ;
                IterQuad = 0 ;
                mp_begin = -1 ;
                memk_begin = 0 ;
                memk = 0 ;

                if ( PrintLevel >= 2 ) printf ("RESTART Sub-CG\n") ;

                /* search direction d = -Zk gsub, gsub = Zk' g, dsub = -gsub
                                 => d =  Zk dsub = SkF (Rk)^{-1} dsub */
                asa_scale0 (dsub, gsubtemp, -ONE, nsub) ;
                asa_copy0 (gsub, gsubtemp, nsub) ;
                asa_copy0 (vsub, dsub, nsub) ;
                asa_trisolve (vsub, Rk, mem, nsub, 1) ;
                /* rearrange and store in wsub */
                mp = SkFlast ;
                j = nsub - (mp+1) ;
                asa_copy0 (wsub, vsub+j, mp+1) ;
                asa_copy0 (wsub+(mp+1), vsub, j) ;
                asa_matvec (d, SkF, wsub, nsub, nfree, 1) ;

                dphi0 = -gsubnorm2 ; /* gsubnorm2 was calculated before */
                dnorm2 = gsubnorm2 ;
            }
            else  /* continue in subspace without restart */
            {
                mlast_sub = (mp_begin + IterSubRestart) % mem ;

                if (IterSubRestart > 0 ) /* not first iteration in subspace  */
                {
                    /* add new column to Yk memory,
                       calculate yty, Sk, Yk and SkYk */
                    spp = mlast_sub*mem ;
                    asa_scale0 (Sk+spp, dsub, alpha, nsub) ;
                    /* yty = (gsubtemp-gsub)'(gsubtemp-gsub),
                       set gsub = gsubtemp */
                    asacg_Yk (Yk+spp, gsub, gsubtemp, &yty, nsub) ;
                    SkYk [mlast_sub] = alpha*(dphi - dphi0) ;
                    scale = SkYk [mlast_sub]/yty ;
                }
                else
                {
                    yty = asa_dot0 (Yk+mlast_sub*mem, Yk+mlast_sub*mem, nsub) ;
                    scale = SkYk [mlast_sub]/yty ;
                }

                /* calculate gsubtemp = H gsub */
                mp = mlast_sub ;
                /* memk = size of the L-BFGS memory in subspace */
                memk = MIN2 (memk_begin + IterSubRestart, mem) ;
                l1 = MIN2 (IterSubRestart, memk) ;
                /* l2 = number of triangular columns in Yk with a zero */
                l2 = memk - l1 ;
                /* l1 = number of dense column in Yk (excluding first) */
                l1++ ;
                l1 = MIN2 (l1, memk) ;

                /* process dense columns */
                for (j = 0; j < l1; j++)
                {
                    mpp = mp*mem ;
                    t = asa_dot0 (Sk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                    tau [mp] = t ;
                    /* update gsubtemp -= t*Yk+mpp */
                    asa_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                    mp-- ;
                    if ( mp < 0 ) mp = mem-1 ;
                }

                /* process columns from triangular (Hessenberg) matrix */
                for (j = 1; j < l2; j++)
                {
                    mpp = mp*mem ;
                    t = asa_dot0 (Sk+mpp, gsubtemp, mp+1)/SkYk[mp] ;
                    tau [mp] = t ;
                    /* update gsubtemp -= t*Yk+mpp */
                    if ( mp == 0 && DenseCol1 )
                    {
                        asa_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                    }
                    else
                    {
                        asa_daxpy0 (gsubtemp, Yk+mpp, -t, MIN2(mp+2,nsub)) ;
                    }
                    mp-- ;
                    if ( mp < 0 ) mp = mem-1 ;
                }
                asa_scale0 (gsubtemp, gsubtemp, scale, nsub) ;

                /* process columns from triangular (Hessenberg) matrix */
                for (j = 1; j < l2; j++)
                {
                    mp++ ;
                    if ( mp == mem ) mp = 0 ;
                    mpp = mp*mem ;
                    if ( mp == 0 && DenseCol1 )
                    {
                        t = asa_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                    }
                    else
                    {
                        t = asa_dot0 (Yk+mpp, gsubtemp,
                                      MIN2(mp+2,nsub))/SkYk[mp] ;
                    }
                    /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                    asa_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, mp+1) ;
                }

                /* process dense columns */
                for (j = 0; j < l1; j++)
                {
                    mp++ ;
                    if ( mp == mem ) mp = 0 ;
                    mpp = mp*mem ;
                    t = asa_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk [mp] ;
                    /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                    asa_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, nsub) ;
                } /* done computing H gsubtemp */

                /* compute d = Zk dsub = SkF (Rk)^{-1} dsub */
                asa_scale0 (dsub, gsubtemp, -ONE, nsub) ;
                asa_copy0 (vsub, dsub, nsub) ;
                asa_trisolve (vsub, Rk, mem, nsub, 1) ;
                /* rearrange and store in wsub */
                mp = SkFlast ;
                j = nsub - (mp+1) ;
                asa_copy0 (wsub, vsub+j, mp+1) ;
                asa_copy0 (wsub+(mp+1), vsub, j) ;

                asa_matvec (d, SkF, wsub, nsub, nfree, 1) ;
                dphi0 = -asa_dot0  (gsubtemp, gsub, nsub) ;

                dnorm2 = asa_dot0 (dsub, dsub,nsub);

            }
        } /* end of subspace search direction */
        else  /* compute the search direction in the full space */
        {
            if ( Restart ) /*restart in fullspace*/
            {
                Restart = FALSE ;
                IterRestart = 0 ;
                IterQuad = 0 ;
                if ( PrintLevel >= 2 ) printf ("RESTART CG\n") ;

                /* set x = xtemp */
                asa_copy (x, xtemp, nfree) ;

                if ( UseMemory )
                {
                   /* set g = gtemp, d = -gtemp,
                      gnorm2 was already computed above */
                   asa_update_dg0 (g, gtemp, d, NULL, nfree) ;
                }
                else
                {
                    /* set g = gtemp, d = -gtemp, compute 2-norm of gtemp*/
                   asa_update_dg0 (g, gtemp, d, &gnorm2, nfree) ;
                }

                dphi0 = -gnorm2 ;
                dnorm2 = gnorm2 ;
                beta = ZERO ;
            }
            else if ( !FirstFull ) /* normal fullspace step*/
            {
                /* set x = xtemp */
                asa_copy (x, xtemp, nfree) ;

                /* ykyk = ||gtemp-g||_2^2, and ykgk = (gtemp-g) dot gnew */
                asa_update_ykyk (g, gtemp, &ykyk, &ykgk, nfree) ;


                dkyk = dphi - dphi0 ;
                if ( Parm->AdaptiveBeta ) t = 2. - ONE/(0.1*QuadTrust + ONE) ;
                else                      t = Parm->theta ;
                beta = (ykgk - t*dphi*ykyk/dkyk)/dkyk ;

                /* faster: initialize dnorm2 = gnorm2 at start, then
                           dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
                           gnorm2 = ||g_{k+1}||^2
                           dnorm2 = ||d_{k+1}||^2
                           dpi = g_{k+1}' d_k */

                /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
                beta = MAX2 (beta, Parm->BetaLower*dphi0/dnorm2) ;

                /* update search direction d = -gtemp + beta*dold */
                if ( UseMemory )
                {
                    /* update search direction d = -gtemp + beta*dold and
                       g = gtemp, and compute 2-norm of d,
                       2-norm of g computed above */
                    dnorm2 = asa_update_dg (g, gtemp, d, beta, NULL, nfree) ;
                }
                else
                {
                    /* update search direction d = -g + beta*dold, and
                       g = gtemp, and compute 2-norms of d and g */
                    dnorm2 = asa_update_dg (g, gtemp, d, beta, &gnorm2, nfree) ;
                }

                dphi0 = -gnorm2 + beta*dphi ;
                if ( Parm->debug ) /* Check that dphi0 = d'g */
                {
                    t = ZERO ;
                    for (i = 0; i < nfree; i++)  t = t + d [i]*g [i] ;
                    if ( fabs(t-dphi0) > Parm->debugtol*fabs(dphi0) )
                    {
                        printf("Warning, dphi0 != d'g!\n");
                        printf("dphi0:%13.6e, d'g:%13.6e\n",dphi0, t) ;
                    }
                }
            }
            else /* FirstFull = TRUE, precondition after leaving subspace */
            {

                /* ykyk = ||gtemp-g||_2^2, and ykgk = (gtemp-g) dot gnew */
                asa_update_ykyk (g, gtemp, &ykyk, &ykgk, nfree) ;

                /* set x = xtemp, g = gtemp */
                asa_update_xy (x, xtemp, g, gtemp, nfree) ;

                mlast_sub = (mp_begin + IterSubRestart) % mem ;
                /* save Sk */
                spp = mlast_sub*mem ;
                asa_scale0 (Sk+spp, dsub, alpha, nsub) ;
                /* calculate yty, save Yk, set gsub = gsubtemp */
                asacg_Yk (Yk+spp, gsub, gsubtemp, &yty, nsub) ;
                ytg = asa_dot0  (Yk+spp, gsub, nsub) ;
                t = alpha*(dphi - dphi0) ;
                SkYk [mlast_sub] = t ;

                /* scale = t/ykyk ; */
                scale = t/yty ;

                /* calculate gsubtemp = H gsub */
                mp = mlast_sub ;
                /* memk = size of the L-BFGS memory in subspace */
                memk = MIN2 (memk_begin + IterSubRestart, mem) ;
                l1 = MIN2 (IterSubRestart, memk) ;
                /* l2 = number of triangular columns in Yk with a zero */
                l2 = memk - l1 ;
                /* l1 = number of dense column in Yk (excluding first) */
                l1++ ;
                l1 = MIN2 (l1, memk) ;

                /* process dense columns */
                for (j = 0; j < l1; j++)
                {
                    mpp = mp*mem ;
                    t = asa_dot0 (Sk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                    tau [mp] = t ;
                    /* update gsubtemp -= t*Yk+mpp */
                    asa_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                    mp-- ;
                    if ( mp < 0 ) mp = mem-1 ;
                }

                /* process columns from triangular (Hessenberg) matrix */
                for (j = 1; j < l2; j++)
                {
                    mpp = mp*mem ;
                    t = asa_dot0 (Sk+mpp, gsubtemp, mp+1)/SkYk[mp] ;
                    tau [mp] = t ;
                    /* update gsubtemp -= t*Yk+mpp */
                    if ( mp == 0 && DenseCol1 )
                    {
                        asa_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                    }
                    else
                    {
                        asa_daxpy0 (gsubtemp, Yk+mpp, -t, MIN2(mp+2,nsub)) ;
                    }
                    mp-- ;
                    if ( mp < 0 ) mp = mem-1 ;
                }
                asa_scale0 (gsubtemp, gsubtemp, scale, nsub) ;

                /* process columns from triangular (Hessenberg) matrix */
                for (j = 1; j < l2; j++)
                {
                    mp++ ;
                    if ( mp == mem ) mp = 0 ;
                    mpp = mp*mem ;
                    if ( mp == 0 && DenseCol1 )
                    {
                        t = asa_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                    }
                    else
                    {
                        t = asa_dot0 (Yk+mpp, gsubtemp,
                                      MIN2(mp+2,nsub))/SkYk[mp];
                    }
                    /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                    asa_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, mp+1) ;
                }

                /* process dense columns */
                for (j = 0; j < l1; j++)
                {
                    mp++ ;
                    if ( mp == mem ) mp = 0 ;
                    mpp = mp*mem ;
                    t = asa_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk [mp] ;
                    /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                    asa_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, nsub) ;
                } /* done computing H gsubtemp */

                /* compute beta */
                dkyk = dphi - dphi0 ;
                if ( Parm->AdaptiveBeta ) t = 2. - ONE/(0.1*QuadTrust + ONE) ;
                else                      t = Parm->theta ;
                t1 = MAX2(ykyk-yty, ZERO) ; /* Theoretically t1 = ykyk-yty */
                scale = (alpha*dkyk)/ykyk ; /* = sigma */
                beta = scale*((ykgk - ytg) - t*dphi*t1/dkyk)/dkyk ;
             /* beta = MAX2 (beta, Parm->BetaLower*dphi0/dnorm2) ; */
                beta = MAX2 (beta, Parm->BetaLower*(dphi0*alpha)/dkyk) ;

                /* compute search direction
                   d = -Zk (H - sigma)ghat - sigma g + beta d

                   Note: d currently contains last 2 terms so only need
                         to add the Zk term. Above gsubtemp = H ghat */

                /* form vsub = sigma ghat - H ghat = sigma ghat - gsubtemp */
                asa_scale0 (vsub, gsubtemp, -ONE, nsub) ;
                asa_daxpy0 (vsub, gsub, scale, nsub) ;
                asa_trisolve (vsub, Rk, mem, nsub, 1) ;

                /* rearrange vsub and store in wsub */
                mp = SkFlast ;
                j = nsub - (mp+1) ;
                asa_copy0 (wsub, vsub+j, mp+1) ;
                asa_copy0 (wsub+(mp+1), vsub, j) ;


                /* save old direct d in gtemp */
                asa_copy (gtemp, d, nfree) ;

                /* d = Zk (sigma - H)ghat */
                asa_matvec (d, SkF, wsub, nsub, nfree, 1) ;

                /* incorporate the new g and old d terms in new d */
                asa_daxpy (d, g, -scale, nfree) ;
                asa_daxpy (d, gtemp, beta, nfree) ;

                gHg = asa_dot0  (gsubtemp, gsub, nsub) ;
                t1 = MAX2(gnorm2 -gsubnorm2, ZERO) ;
                dphi0 = -gHg - scale*t1 + beta*dphi ;
                /* dphi0 = asa_dot (d, g, nfree) could be inaccurate */
                dnorm2 = asa_dot (d, d, nfree) ;

            }  /* end of preconditioned step */
        }  /* search direction has been computed */

        /* test for slow convergence */
        if ( (f < fbest) || (gnorm2 < gbest) )
        {
            nslow = 0 ;
            if ( f < fbest ) fbest = f ;
            if ( gnorm2 < gbest ) gbest = gnorm2 ;
        }
        else nslow++ ;
        if ( nslow > slowlimit )
        {
            Com->nslow = nslow ;
            status = 9 ;
            goto Exit ;
        }

        if ( PrintLevel >= 2 )
        {
            printf ("\niter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e "
                    "memk: %i Subspace: %i\n",
                    (int) iter, f, pgnorm, ginorm, memk, Subspace) ;
        }

        if ( Parm->debug )
        {
            if ( f > Com->f0 + Parm->debugtol*Ck )
            {
                status = 8 ;
                goto Exit ;
            }
        }

        if ( dphi0 > ZERO )
        {
           status = 5 ;
           goto Exit ;
        }


    }
    status = 2 ;

Exit:
    if ( status < -2 )
    {
        for (j = nfree; j < n; j++) g [j] = gtemp [j] ;
    }
    else
    {
        pgnorm = ZERO ;
        for (j = 0; j < n; j++)
        {
            xj = xtemp [j] ;
            x [j] = xj ;
            gj = gtemp [j] ;
            g [j] = gj ;
            xg = xj - gj ;
            if      ( xg > hi [j] ) xp = hi [j] - xj ;
            else if ( xg < lo [j] ) xp = lo [j] - xj ;
            else                    xp = -gj ;
            pgnorm = MAX2 (pgnorm, fabs (xp)) ;
            pg [j] = xp ;
        }
    }

Exit1:
    Com->pgnorm = pgnorm ;
    Com->ginorm = ginorm ;
    Com->f = f ;
    Com->f_debug = f ;
    Com->cgfunc += Com->nf - nf ;
    Com->cggrad += Com->ng - ng ;
    Com->cgiter += iter ;
    if ( PrintLevel >= 3 )
    {
        printf ("\niter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n",
                (int) iter, f, pgnorm, ginorm) ;
    }
    if ( PrintLevel >= 2 )
    {
        printf ("\nCG Termination status: %i\n", status) ;
        if ( status == -5 )
        {
            printf ("ginorm < tau2*pgnorm without hitting boundary\n") ;
        }
        if ( status == -4 )
        {
            printf ("ginorm >= tau2*pgnorm, many x components hit boundary\n") ;
        }
        else if ( status == -3 )
        {
            printf ("ginorm >= tau2*pgnorm, one x component hits boundary\n") ;
        }
        printf ("proj gradient max norm: %13.6e\n", pgnorm) ;
        printf ("function value:         %13.6e\n", f) ;
        printf ("cg iterations:          %13.6e\n", (double) iter) ;
        printf ("function evaluations:   %13.6e\n", (double) Com->nf - nf) ;
        printf ("gradient evaluations:   %13.6e\n", (double) Com->ng - ng) ;
    }
    if ( asaParm->DynamicMemory ) free (work) ;
    return (status) ;
}

/* =========================================================================
   === asa_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   ========================================================================= */
 int asa_Wolfe
(
    double       alpha , /* stepsize */
    double           f , /* function value associated with stepsize alpha */
    double        dphi , /* derivative value associated with stepsize alpha */
    asa_com        *Com
)
{
    if ( dphi >= Com->wolfe_lo )
    {

        /* test original Wolfe conditions */
        if ( f - Com->f0 <= alpha*Com->wolfe_hi )
        {
            if ( Com->cgParm->PrintLevel >= 4 )
            {
                printf ("Wolfe conditions hold\n") ;
/*                printf ("wolfe f: %14.6e f0: %14.6e dphi: %14.6e\n",
                         f, Com->f0, dphi) ;*/
            }
            return (1) ;
        }
        /* test approximate Wolfe conditions */
        else if ( Com->AWolfe )
        {
            if ( (f <= Com->fpert) && (dphi <= Com->awolfe_hi) )
            {
                if ( Com->cgParm->PrintLevel >= 4 )
                {
                    printf ("Approximate Wolfe conditions hold\n") ;
/*                    printf ("f: %14.6e fpert: %14.6e dphi: %14.6e awolf_hi: "
                            "%14.6e\n", f, Com->fpert, dphi, Com->awolfe_hi) ;*/
                }
                return (1) ;
            }
        }
    }
    return (0) ;
}

/* =========================================================================
   === asa_tol =============================================================
   =========================================================================
   Check for convergence
   ========================================================================= */
int asa_tol
(
    double      pgnorm, /* projected gradient sup-norm */
    asa_com       *Com
)
{
    /*StopRule = T => |grad|_infty <=max (tol, |grad|_infty*StopFac)
                 F => |grad|_infty <= tol*(1+|f|)) */
    if ( Com->asaParm->StopRule )
    {
        if ( pgnorm <= Com->tol ) return (1) ;
    }
    else if ( pgnorm <= Com->tol*(ONE + fabs (Com->f)) ) return (1) ;
    return (0) ;
}

/* =========================================================================
   === asa_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   Return:
      -2 (function nan)
       0 (Wolfe or approximate Wolfe conditions satisfied)
       3 (slope always negative in line search)
       4 (number line search iterations exceed nline)
       6 (excessive updating of eps)
       7 (Wolfe conditions never satisfied)

   ========================================================================= */
 int asa_line
(
    asa_com       *Com  /* cg com structure */
)
{
    int AWolfe, iter, ngrow, PrintLevel, qb, qb0, status, toggle ;
    double alpha, a, a1, a2, b, bmin, B, da, db, d0, d1, d2, dB, df, f, fa, fb,
           fB, a0, b0, da0, db0, fa0, fb0, width, rho, minstep, phiold ;
    const char *s1, *s2, *fmt1, *fmt2 ;
    asacg_parm *Parm ;

    AWolfe = Com->AWolfe ;
    Parm = Com->cgParm ;
    PrintLevel = Parm->PrintLevel ;
    minstep = Com->minstep ;

    if ( PrintLevel >= 2 )
    {
        if ( AWolfe )
        {
            printf ("Approximate Wolfe line search\n") ;
            printf ("=============================\n") ;
        }
        else
        {
            printf ("Wolfe line search\n") ;
            printf ("=================\n") ;
        }
    }

    /* evaluate function or gradient at Com->alpha (starting guess) */
    if ( Com->alpha > minstep )
    {
        if ( !Com->minflag )
        {
            asa_maxstep (Com->x, Com->d, Com) ;
            minstep = Com->minstep ;
            if ( Com->alpha > minstep )
            {
               Com->alpha = minstep ;
               Com->QuadOK = FALSE ;
            }
        }
        else
        {
            Com->alpha = minstep ;
            Com->QuadOK = FALSE ;
        }
    }
    b = Com->alpha ;

    if ( Com->QuadOK )
    {
        status = asacg_evaluate ("fg", "y", Com) ;
      b = Com->alpha ;
        fb = Com->f ;
        if ( !AWolfe ) fb -= b*Com->wolfe_hi ;
        qb = TRUE ; /* function value at b known */
    }
    else
    {
        status = asacg_evaluate ("g", "y", Com) ;
      b = Com->alpha ;
        qb = FALSE ;
    }
    if ( status ) return (status) ; /* function is undefined */

    if ( AWolfe )
    {
        db = Com->df ;
        d0 = da = Com->df0 ;
    }
    else
    {
        db = Com->df - Com->wolfe_hi ;
        d0 = da = Com->df0 - Com->wolfe_hi ;
    }
    a = ZERO ;
    a1 = ZERO ;
    d1 = d0 ;
    fa = Com->f0 ;
    if ( PrintLevel >= 2 )
    {
        fmt1 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb: %13.6e "
               "da: %13.6e db: %13.6e\n" ;
        fmt2 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb:  x.xxxxxxxxxx "
               "da: %13.6e db: %13.6e\n" ;
        if ( Com->QuadOK ) s2 = "OK" ;
        else               s2 = "" ;
        if ( qb ) printf (fmt1, "start    ", s2, a, b, fa, fb, da, db);
        else      printf (fmt2, "start    ", s2, a, b, fa, da, db) ;
    }

    /* if a quadratic interpolation step performed, check Wolfe conditions */
    if ( (Com->QuadOK) && (Com->f <= Com->f0) )
    {
        if ( asa_Wolfe (b, Com->f, Com->df, Com) ) return (0) ;
    }

    /* if a Wolfe line search and the Wolfe conditions have not been satisfied*/
    if ( !AWolfe ) Com->Wolfe = TRUE ;

    /*Find initial interval [a,b] such that
      da <= 0, db >= 0, fa <= fpert = [(f0 + eps*fabs (f0)) or (f0 + eps)] */
    rho = Com->rho ;
    ngrow = 1 ;
    while ( db < ZERO )
    {
        if ( !qb )
        {
            status = asacg_evaluate ("f", "n", Com) ;
            if ( status ) return (status) ;
            if ( AWolfe ) fb = Com->f ;
            else          fb = Com->f - b*Com->wolfe_hi ;
            qb = TRUE ;
        }
        if ( fb > Com->fpert ) /* contract interval [a, b] */
        {
            status = asacg_contract (&a, &fa, &da, &b, &fb, &db, Com) ;
            if ( status == 0 ) return (0) ;   /* Wolfe conditions hold */
            if ( status == -2 ) goto Line ; /* db >= 0 */
            if ( Com->neps > Parm->neps ) return (6) ;
        }

        /* expansion phase */
        ngrow++ ;
        if ( ngrow > Parm->ntries ) return (3) ;
        /* update interval (a replaced by b) */
        a = b ;
        fa = fb ;
        da = db ;
        /* store old values of a and corresponding derivative */
        d2 = d1 ;
        d1 = da ;
        a2 = a1 ;
        a1 = a ;

        bmin = rho*b ;
        if ( (ngrow == 2) || (ngrow == 3) || (ngrow == 6) )
        {
            if ( d1 > d2 )
            {
                if ( ngrow == 2 )
                {
                    b = a1 - (a1-a2)*(d1/(d1-d2)) ;
                }
                else
                {
                    if ( (d1-d2)/(a1-a2) >= (d2-d0)/a2 )
                    {
                        /* convex derivative, secant overestimates minimizer */
                        b = a1 - (a1-a2)*(d1/(d1-d2)) ;
                    }
                    else
                    {
                        /* concave derivative, secant underestimates minimizer*/
                        b = a1 - Parm->SecantAmp*(a1-a2)*(d1/(d1-d2)) ;
                    }
                }
                /* safeguard growth */
                b = MIN2 (b, Parm->ExpandSafe*a1) ;
            }
            else rho *= Parm->RhoGrow ;
        }
        else rho *= Parm->RhoGrow ;
        b = MAX2 (bmin, b) ;
        Com->alphaold = Com->alpha ;
        if ( b > minstep )
        {
            if ( !Com->minflag )
            {
                asa_maxstep (Com->x, Com->d, Com) ;
                minstep = Com->minstep ;
                if ( b > minstep ) b = minstep ;
            }
            else   b = minstep ;
        }
        Com->alpha = b ;
        if ( b != Com->alphaold )
        {
            status = asacg_evaluate ("g", "p", Com) ;
  b = Com->alpha ;
            if ( status ) return (status) ;
            qb = FALSE ;
            if ( AWolfe ) db = Com->df ;
            else          db = Com->df - Com->wolfe_hi ;
            if ( PrintLevel >= 3 )
            {
                if ( Com->QuadOK ) s2 = "OK" ;
                else               s2 = "" ;
                printf (fmt2, "expand   ", s2, a, b, fa, da, db) ;
            }
        }
        else /* a new constraint is active */
        {
            do /* while statement */
            {
                Com->alphaold = Com->alpha ;
                phiold = Com->f ;
                if ( Com->alpha < Com->maxstep )
                {
                    Com->alpha = Com->rho*Com->alphaold ;
                    status = asacg_evaluate ("f", "e", Com) ;
                    if ( status ) return (status) ;
                }
            } while ( Com->f < phiold ) ;
            if ( Com->alphaold == minstep )
            {
                Com->alpha = minstep ;
                asa_step (Com->xtemp, Com->x, Com->d, minstep, Com->nfree) ;
                status = -3 ;
            }
            else
            {
                Com->alpha = Com->alphaold ;
                status = asacg_evaluate ("g", "e", Com) ;
                if ( status ) return (status) ;
                status = -4 ;
            }
            Com->f = phiold ;
            return (status) ;
        }

    }

    /* we now have fa <= fpert, da <= 0, db >= 0 */
Line:
    toggle = 0 ;
    width = b - a ;
    qb0 = FALSE ;
    for (iter = 0; iter < Parm->nline; iter++)
    {
        /* determine the next iterate */
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
        {
            Com->QuadOK = TRUE ;
            if ( Com->UseCubic && qb )
            {
                s1 = "cubic    " ;
                alpha = asacg_cubic (a, fa, da, b, fb, db) ;
                if ( alpha < ZERO ) /* use secant method */
                {
                    s1 = "secant   " ;
                    if      ( -da < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if ( da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                 alpha = -1. ;
                }
            }
            else
            {
                s1 = "secant   " ;
                if      ( -da < db ) alpha = a - (a-b)*(da/(da-db)) ;
                else if ( da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                else                 alpha = -1. ;
            }
            width = Parm->gamma*(b - a) ;
        }
        else if ( toggle == 1 ) /* iteration based on smallest value*/
        {
            Com->QuadOK = TRUE ;
            if ( Com->UseCubic )
            {
                s1 = "cubic    " ;
                if ( Com->alpha == a ) /* a is most recent iterate */
                {
                    alpha = asacg_cubic (a0, fa0, da0, a, fa, da) ;
                }
                else if ( qb0 )        /* b is most recent iterate */
                {
                    alpha = asacg_cubic (b, fb, db, b0, fb0, db0) ;
                }
                else alpha = -1. ;

                /* if alpha no good, use cubic between a and b */
                if ( (alpha <= a) || (alpha >= b) )
                {
                    if ( qb ) alpha = asacg_cubic (a, fa, da, b, fb, db) ;
                    else alpha = -1. ;
                }

                /* if alpha still no good, use secant method */
                if ( alpha < ZERO )
                {
                    s1 = "secant   " ;
                    if      ( -da < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if ( da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                 alpha = -1. ;
                }
            }
            else /* ( use secant ) */
            {
                s1 = "secant   " ;
                if ( (Com->alpha == a) && (da > da0) ) /* use a0 if possible */
                {
                    alpha = a - (a-a0)*(da/(da-da0)) ;
                }
                else if ( db < db0 )                   /* use b0 if possible */
                {
                    alpha = b - (b-b0)*(db/(db-db0)) ;
                }
                else /* secant based on a and b */
                {
                    if      ( -da < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if ( da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                 alpha = -1. ;
                }

                if ( (alpha <= a) || (alpha >= b) )
                {
                    if      ( -da < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if ( da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                 alpha = -1. ;
                }
            }
        }
        else
        {
            alpha = .5*(a+b) ; /* use bisection if b-a decays slowly */
            s1 = "bisection" ;
            Com->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) )
        {
            alpha = .5*(a+b) ;
            s1 = "bisection" ;
            if ( (alpha == a) || (alpha == b) ) return (7) ;
            Com->QuadOK = FALSE ; /* bisection was used */
        }

        if ( toggle == 0 ) /* save values for next iteration */
        {
            a0 = a ;
            b0 = b ;
            da0 = da ;
            db0 = db ;
            fa0 = fa ;
            if ( qb )
            {
                fb0 = fb ;
                qb0 = TRUE ;
            }
        }

        toggle++ ;
        if ( toggle > 2 ) toggle = 0 ;

        Com->alpha = alpha ;
        status = asacg_evaluate ("fg", "n", Com) ;
        if ( status ) return (status) ;
        Com->alpha = alpha ;
        f = Com->f ;
        df = Com->df ;
        if ( Com->QuadOK )
        {
            if ( asa_Wolfe (alpha, f, df, Com) )
            {
                if ( PrintLevel >= 3 )
                {
                    printf ("             a: %13.6e f: %13.6e df: %13.6e %1s\n",
                             alpha, f, df, s1) ;
                }
                return (0) ;
            }
        }
        if ( !AWolfe )
        {
            f -= alpha*Com->wolfe_hi ;
            df -= Com->wolfe_hi ;
        }
        if ( df >= ZERO )
        {
            b = alpha ;
            fb = f ;
            db = df ;
            qb = TRUE ;
        }
        else if ( f <= Com->fpert )
        {
            a = alpha ;
            da = df ;
            fa = f ;
        }
        else
        {
            B = b ;
            if ( qb ) fB = fb ;
            dB = db ;
            b = alpha ;
            fb = f ;
            db = df ;
            /* contract interval [a, alpha] */
            status = asacg_contract (&a, &fa, &da, &b, &fb, &db, Com) ;
            if ( status == 0 ) return (0) ;
            if ( status == -1 ) /* eps reduced, use [a, b] = [alpha, b] */
            {
                if ( Com->neps > Parm->neps ) return (6) ;
                a = b ;
                fa = fb ;
                da = db ;
                b = B ;
                if ( qb ) fb = fB ;
                db = dB ;
            }
            else qb = TRUE ;
        }
        if ( PrintLevel >= 3 )
        {
            if ( Com->QuadOK ) s2 = "OK" ;
            else               s2 = "" ;
            if ( !qb ) printf (fmt2, s1, s2, a, b, fa, da, db) ;
            else       printf (fmt1, s1, s2, a, b, fa, fb, da, db) ;
        }
    }
    return (4) ;
}

/* =========================================================================
   ==== asacg_contract =====================================================
   =========================================================================
   The input for this routine is an interval [a, b] with the property that
   fa <= fpert, da >= 0, db >= 0, and fb >= fpert. The returned status is

  11  function or derivative not defined
   0  if the Wolfe conditions are satisfied
  -1  if a new value for eps is generated with the property that for the
      corresponding fpert, we have fb <= fpert
  -2  if a subinterval, also denoted [a, b], is generated with the property
      that fa <= fpert, da >= 0, and db <= 0

   NOTE: The input arguments are unchanged when status = -1
   ========================================================================= */
 int asacg_contract
(
    double    *A, /* left side of bracketing interval */
    double   *fA, /* function value at a */
    double   *dA, /* derivative at a */
    double    *B, /* right side of bracketing interval */
    double   *fB, /* function value at b */
    double   *dB, /* derivative at b */
    asa_com  *Com  /* cg com structure */
)
{
    int AWolfe, iter, PrintLevel, toggle, status ;
    double a, alpha, b, old, da, db, df, d1, dold, f, fa, fb, f1, fold,
           t, width ;
    const char *s ;
    asacg_parm *Parm ;

    AWolfe = Com->AWolfe ;
    Parm = Com->cgParm ;
    PrintLevel = Parm->PrintLevel ;
    a = *A ;
    fa = *fA ;
    da = *dA ;
    b = *B ;
    fb = *fB ;
    db = *dB ;
    f1 = fb ;
    d1 = db ;
    toggle = 0 ;
    width = ZERO ;
    for (iter = 0; iter < Parm->nshrink; iter++)
    {
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
        {
            /* cubic based on bracketing interval */
            alpha = asacg_cubic (a, fa, da, b, fb, db) ;
            toggle = 0 ;
            width = Parm->gamma*(b-a) ;
            if ( iter ) Com->QuadOK = TRUE ; /* at least 2 cubic iterations */
        }
        else if ( toggle == 1 )
        {
            Com->QuadOK = TRUE ;
            /* cubic based on most recent iterate and smallest value */
            if ( old < a ) /* a is most recent iterate */
            {
                alpha = asacg_cubic (a, fa, da, old, fold, dold) ;
            }
            else           /* b is most recent iterate */
            {
                alpha = asacg_cubic (a, fa, da, b, fb, db) ;
            }
        }
        else
        {
            alpha = .5*(a+b) ; /* use bisection if b-a decays slowly */
            Com->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) )
        {
            alpha = .5*(a+b) ;
            Com->QuadOK = FALSE ; /* bisection was used */
        }

        toggle++ ;
        if ( toggle > 2 ) toggle = 0 ;

        Com->alpha = alpha ;
        status = asacg_evaluate ("fg", "n", Com) ;
        if ( status ) return (status) ;
        f = Com->f ;
        df = Com->df ;

        if ( Com->QuadOK )
        {
            if ( asa_Wolfe (alpha, f, df, Com) ) return (0) ;
        }
        if ( !AWolfe )
        {
            f -= alpha*Com->wolfe_hi ;
            df -= Com->wolfe_hi ;
        }
        if ( df >= ZERO )
        {
            *B = alpha ;
            *fB = f ;
            *dB = df ;
            *A = a ;
            *fA = fa ;
            *dA = da ;
            return (-2) ;
        }
        if ( f <= Com->fpert ) /* update a using alpha */
        {
            old = a ;
            a = alpha ;
            fold = fa ;
            fa = f ;
            dold = da ;
            da = df ;
        }
        else                     /* update b using alpha */
        {
            old = b ;
            b = alpha ;
            fb = f ;
            db = df ;
        }
        if ( PrintLevel >= 3 )
        {
            if ( Com->QuadOK ) s = "OK" ;
            else               s = "" ;
            printf ("contract  %2s a: %13.6e b: %13.6e fa: %13.6e fb: "
                    "%13.6e da: %13.6e db: %13.6e\n", s, a, b, fa, fb, da, db) ;
        }
    }

    /* see if the cost is small enough to change the PertRule */
    if ( fabs (fb) <= Com->SmallCost ) Com->PertRule = FALSE ;

    /* increase eps if slope is negative after Parm->nshrink iterations */
    t = Com->f0 ;
    if ( Com->PertRule )
    {
        if ( t != ZERO )
        {
            Com->eps = Parm->egrow*(f1-t)/fabs (t) ;
            Com->fpert = t + fabs (t)*Com->eps ;
        }
        else Com->fpert = 2.*f1 ;
    }
    else
    {
        Com->eps = Parm->egrow*(f1-t) ;
        Com->fpert = t + Com->eps ;
    }
    if ( PrintLevel >= 2 )
    {
        printf ("--increase eps: %e fpert: %e\n", Com->eps, Com->fpert) ;
    }
    Com->neps++ ;
    return (-1) ;
}
/* =========================================================================
   ==== asacg_cubic ========================================================
   =========================================================================
   Compute the minimizer of a Hermite cubic. If the computed minimizer
   outside [a, b], return -1 (it is assumed that a >= 0).
   ========================================================================= */
 double asacg_cubic
(
    double  a,
    double fa, /* function value at a */
    double da, /* derivative at a */
    double  b,
    double fb, /* function value at b */
    double db  /* derivative at b */
)
{
    double c, d1, d2, delta, t, v, w ;
    delta = b - a ;
    if ( delta == ZERO ) return (a) ;
    v = da + db - 3.*(fb-fa)/delta ;
    t = v*v - da*db ;
    if ( t < ZERO ) /* complex roots, use secant method */
    {
         if ( fabs (da) < fabs (db) ) c = a - (a-b)*(da/(da-db)) ;
         else if ( da != db )         c = b - (a-b)*(db/(da-db)) ;
         else                         c = -1 ;
         return (c) ;
    }

    if ( delta > ZERO ) w = sqrt(t) ;
    else                w =-sqrt(t) ;
    d1 = da + v - w ;
    d2 = db + v + w ;
    if ( (d1 == ZERO) && (d2 == ZERO) ) return (-1.) ;
    if ( fabs (d1) >= fabs (d2) ) c = a + delta*da/d1 ;
    else                          c = b - delta*db/d2 ;
    return (c) ;
}



/* =========================================================================
   ==== asacg_Yk ==============================================================
   =========================================================================
   Compute y = gnew - gold, set gold = gnew, compute y'y
   ========================================================================= */
 void asacg_Yk
(
    double    *y, /*output vector */
    double *gold, /* initial vector */
    double *gnew, /* search direction */
    double  *yty, /* y'y */
    ASA_INT        n  /* length of the vectors */
)
{
    ASA_INT n5, i ;
    double s, t ;
    n5 = n % 5 ;
    if ( (y != NULL) && (yty == NULL) )
    {
        for (i = 0; i < n5; i++)
        {
            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
        }
        for (; i < n; )
        {
            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;
        }
    }
    else if ( (y == NULL) && (yty != NULL) )
    {
        s = ZERO ;
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;
        }
        *yty = s ;
    }
    else
    {
        s = ZERO ;
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;
        }
        *yty = s ;
    }

    return ;
}

/* =========================================================================
   ==== asacg_evaluate =======================================================
   Evaluate the function and/or gradient.  Also, possibly check if either is nan
   and if so, then reduce the stepsize. Only used at the start of an iteration.
   Return:
      11 (function nan)
       0 (successful evaluation)
   =========================================================================*/

 int asacg_evaluate
(
    const char    *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    const char     *nan, /* y means check function/derivative values for nan */
    asa_com   *Com
)
{
    ASA_INT n ;
    int i ;
    double alpha, *d, *g, *gtemp, *x, *xtemp ;
    asacg_parm *Parm ;
	asa_objective *user ;
    Parm = Com->cgParm ;
    n = Com->n ;
    x = Com->x ;
    g = Com->g ;
    d = Com->d ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    alpha = Com->alpha ;
    user = Com->user ;
    user->x = xtemp ;
    user->g = gtemp ;

    /* check to see if values are nan */
    if ( !strcmp (nan, "y") || !strcmp (nan, "p") || !strcmp (nan, "e"))
    {
        if ( !strcmp (what, "f") ) /* compute function */
        {
            if ( !strcmp (nan, "e") ) /* alpha is expansion step in cg_line */
            {
                asa_project (xtemp, x, d, alpha, Com) ;
            }
            else
            {
                asa_step (xtemp, x, d, alpha, Com->nfree) ;
            }
            /* provisional function value */
            if ( Com->DimReduce )
            {
               /* Expand xtemp to the full space*/
               asa_expandx (xtemp, Com) ;
               /* Evaluate function */
               user->ifree = Com->ifree ;
               user->nfree = Com->nfree ;
               if(Com->objpoint!=NULL)
               Com->f = Com->objpoint->myvalue(user) ;
               if(Com->avec_objpoint!=NULL)
               Com->f = Com->avec_objpoint->myvalue(user) ;
               /* Shrink xtemp to the reduced space */
               asa_shrinkx (xtemp, Com) ;
            }
            else
            {
               /* Evaluate function */
               user->ifree = NULL ;
               user->nfree = Com->n ;
               if(Com->objpoint!=NULL)
               Com->f = Com->objpoint->myvalue(user) ;
               if(Com->avec_objpoint!=NULL)
               Com->f = Com->avec_objpoint->myvalue(user) ;
            }
            Com->nf++ ;
            /* reduce stepsize if function value is nan */
            if ( (Com->f != Com->f) || (Com->f >= ASA_INF) || (Com->f <= -ASA_INF) )
            {
                for (i = 0; i < Parm->ntries; i++)
                {
                    /* contract from good alpha */
                    if ( !strcmp (nan, "p") || !strcmp (nan, "e") )
                    {
                        alpha = Com->alphaold + .8*(alpha - Com->alphaold) ;
                    }
                    else                      /* multiply by nan_decay */
                    {
                        alpha *= Parm->nan_decay ;
                    }
                    /* initial alpha is expansion step in cg_line */
                    if ( !strcmp (nan, "e") )
                    {
                        asa_project (xtemp, x, d, alpha, Com) ;
                    }
                    else
                    {
                        asa_step (xtemp, x, d, alpha, Com->nfree) ;
                    }
                    /* provisional function value */
                    if ( Com->DimReduce )
                    {
                       /* Expand xtemp to the full space*/
                       asa_expandx (xtemp, Com) ;
                       /* Evaluate function */
                       user->ifree = Com->ifree ;
                       user->nfree = Com->nfree ;
                       if(Com->objpoint!=NULL)
                       Com->f =Com->objpoint->myvalue(user) ;
                       if(Com->avec_objpoint!=NULL)
                       Com->f =Com->avec_objpoint->myvalue(user) ;
                       /* Shrink xtemp to the reduced space */
                       asa_shrinkx (xtemp, Com) ;
                    }
                    else
                    {
                       /* Evaluate function */
                       user->ifree = NULL ;
                       user->nfree = Com->n ;
                       if(Com->objpoint!=NULL)
                       Com->f = Com->objpoint->myvalue(user) ;
                       if(Com->avec_objpoint!=NULL)
                       Com->f = Com->avec_objpoint->myvalue(user) ;
                    }
                    Com->nf++ ;
                    if ( (Com->f == Com->f) && (Com->f < ASA_INF) &&
                         (Com->f > -ASA_INF) ) break ;
                }
                if ( i == Parm->ntries ) return (11) ;
            }
            Com->alpha = alpha ;
        }
        else if ( !strcmp (what, "g") ) /* compute gradient */
        {
            if ( !strcmp (nan, "e") ) /* alpha is expansion step in cg_line */
            {
                asa_project (xtemp, x, d, alpha, Com) ;
            }
            else
            {
                asa_step (xtemp, x, d, alpha, Com->nfree) ;
            }
            if ( Com->DimReduce )
            {
                /* Expand xtemp to the full space*/
                asa_expandx (xtemp, Com) ;
                /* Evaluate gradient */
                user->ifree = Com->ifree ;
                user->nfree = Com->nfree ;
               if(Com->objpoint!=NULL) 
               Com->objpoint->mygrad(user) ;
               if(Com->avec_objpoint!=NULL) 
               Com->avec_objpoint->mygrad(user) ;
                /* Shrink x and g to the reduced space */
                asa_shrinkxg (xtemp, gtemp, Com) ;
            }
            else
            {
               /* Evaluate gradient */
               user->ifree = NULL ;
               user->nfree = Com->n ;
               if(Com->objpoint!=NULL) 
               Com->objpoint->mygrad(user) ;
               if(Com->avec_objpoint!=NULL) 
               Com->avec_objpoint->mygrad(user) ;
            }
            Com->ng++ ;
            Com->df = asa_dot (gtemp, d, Com->nfree) ;
            /* reduce stepsize if derivative is nan */
            if ( (Com->df != Com->df) || (Com->df >= ASA_INF) || (Com->df <= -ASA_INF) )
            {
                for (i = 0; i < Parm->ntries; i++)
                {
                    /* contract from good alpha */
                    if ( !strcmp (nan, "p") || !strcmp (nan, "e") )
                    {
                        alpha = Com->alphaold + .8*(alpha - Com->alphaold) ;
                    }
                    else  /* multiply by nan_decay */
                    {
                        alpha *= Parm->nan_decay ;
                    }
                    /* initial alpha is expansion step in cg_line */
                    if ( !strcmp (nan, "e") )
                    {
                        asa_project (xtemp, x, d, alpha, Com) ;
                    }
                    else
                    {
                        asa_step (xtemp, x, d, alpha, Com->nfree) ;
                    }
                    if ( Com->DimReduce )
                    {
                        /* Expand xtemp to the full space*/
                        asa_expandx (xtemp, Com) ;
                        /* Evaluate gradient */
                        user->ifree = Com->ifree ;
                        user->nfree = Com->nfree ;
                        if(Com->objpoint!=NULL)
                        Com->objpoint->mygrad(user) ;
                        if(Com->avec_objpoint!=NULL)
                        Com->avec_objpoint->mygrad(user) ;
                        /* Shrink x and g to the reduced space */
                        asa_shrinkxg (xtemp, gtemp, Com) ;
                    }
                    else
                    {
                        /* Evaluate gradient */
                        user->ifree = NULL ;
                        user->nfree = Com->n ;
                        if(Com->objpoint!=NULL)
                        Com->objpoint->mygrad(user) ;
                        if(Com->avec_objpoint!=NULL)
                        Com->avec_objpoint->mygrad(user) ;
                    }
                    Com->ng++ ;
                    Com->df = asa_dot (gtemp, d, Com->nfree) ;
                    if ( (Com->df == Com->df) && (Com->df < ASA_INF) &&
                         (Com->df > -ASA_INF) ) break ;
                }
                if ( i == Parm->ntries ) return (11) ;
                Com->rho = Parm->nan_rho ;
            }
            else Com->rho = Parm->rho ;
            Com->alpha = alpha ;
        }
        else                            /* compute function and gradient */
        {
            asa_step (xtemp, x, d, alpha, Com->nfree) ;
            if ( Com->DimReduce )
            {
                /* Expand xtemp to the full space*/
                asa_expandx (xtemp, Com) ;
                /* Evaluate function and gradient */
                user->ifree = Com->ifree ;
                user->nfree = Com->nfree ;
                if ( Com->valgrad != NULL )
                {
                   //Com->f = Com->valgrad (user) ;
                }
                else
                {
                   if(Com->objpoint!=NULL){
                   Com->objpoint->mygrad(user) ;
                   Com->f = Com->objpoint->myvalue(user) ;
                   }
                   if(Com->avec_objpoint!=NULL){
                   Com->avec_objpoint->mygrad(user) ;
                   Com->f = Com->avec_objpoint->myvalue(user) ;
                   }
                }
                /* Shrink xtemp and gtemp to the reduced space */
                asa_shrinkxg (xtemp, gtemp, Com) ;
            }
            else
            {
                /* Evaluate function and gradient */
                user->ifree = NULL ;
                user->nfree = Com->n ;
                if ( Com->valgrad != NULL )
                {
                  // Com->f = Com->valgrad (user) ;
                }
                else
                {
                   if(Com->objpoint!=NULL){
                   Com->objpoint->mygrad(user) ;
                   Com->f = Com->objpoint->myvalue(user) ;
                   }
                   if(Com->avec_objpoint!=NULL){
                   Com->avec_objpoint->mygrad(user) ;
                   Com->f = Com->avec_objpoint->myvalue(user) ;
                   }
                }
            }
            Com->df = asa_dot (gtemp, d, Com->nfree) ;
            Com->nf++ ;
            Com->ng++ ;
            /* reduce stepsize if function value or derivative is nan */
            if ( (Com->df !=  Com->df) || (Com->f != Com->f) ||
                 (Com->df >=  ASA_INF)     || (Com->f >= ASA_INF)    ||
                 (Com->df <= -ASA_INF)     || (Com->f <= -ASA_INF))
            {
                for (i = 0; i < Parm->ntries; i++)
                {
                    if ( !strcmp (nan, "p") ) /* contract from good alpha */
                    {
                        alpha = Com->alphaold + .8*(alpha - Com->alphaold) ;
                    }
                    else                      /* multiply by nan_decay */
                    {
                        alpha *= Parm->nan_decay ;
                    }
                    asa_step (xtemp, x, d, alpha, Com->nfree) ;
                    if ( Com->DimReduce )
                    {
                        /* Expand xtemp to the full space*/
                        asa_expandx (xtemp, Com) ;
                        /* Evaluate function and gradient */
                        user->ifree = Com->ifree ;
                        user->nfree = Com->nfree ;
                        if ( Com->valgrad != NULL )
                        {
                          // Com->f = Com->valgrad (user) ;
                        }
                        else
                        {
                           if(Com->objpoint!=NULL){
                           Com->objpoint->mygrad(user) ;
                           Com->f = Com->objpoint->myvalue(user) ;
                           }
                           if(Com->avec_objpoint!=NULL){
                            Com->avec_objpoint->mygrad(user) ;
                            Com->f = Com->avec_objpoint->myvalue(user) ;
                           }
                        }
                        /* Shrink xtemp and gtemp to the reduced space */
                        asa_shrinkxg (xtemp, gtemp, Com) ;
                    }
                    else
                    {
                        /* Evaluate function and gradient */
                        user->ifree = NULL ;
                        user->nfree = Com->n ;
                        if ( Com->valgrad != NULL )
                        {
                           //Com->f = Com->valgrad (user) ;
                        }
                        else
                        {
                           if(Com->objpoint!=NULL){
                           Com->objpoint->mygrad(user) ;
                           Com->f = Com->objpoint->myvalue(user) ;
                           }
                           if(Com->avec_objpoint!=NULL){
                           Com->avec_objpoint->mygrad(user) ;
                           Com->f = Com->avec_objpoint->myvalue(user) ;
                           }
                        }
                    }
                    Com->df = asa_dot (gtemp, d, Com->nfree) ;
                    Com->nf++ ;
                    Com->ng++ ;
                    if ( (Com->df == Com->df) && (Com->f == Com->f) &&
                         (Com->df <  ASA_INF)     && (Com->f <  ASA_INF)    &&
                         (Com->df > -ASA_INF)     && (Com->f > -ASA_INF) ) break ;
                }
                if ( i == Parm->ntries ) return (11) ;
                Com->rho = Parm->nan_rho ;
            }
            else Com->rho = Parm->rho ;
            Com->alpha = alpha ;
        }
    }
    else                                /* evaluate without nan checking */
    {
        if ( !strcmp (what, "fg") )     /* compute function and gradient */
        {
            if ( alpha == ZERO )        /* evaluate at x */
            {
                user->x = x ;
                user->g = g ;
                if ( Com->DimReduce )
                {
                    /* Expand x to the full space*/
                    asa_expandx (x, Com) ;
                    /* Evaluate function and gradient */
                    user->ifree = Com->ifree ;
                    user->nfree = Com->nfree ;
                    if ( Com->valgrad != NULL )
                    {
                        //Com->f = Com->valgrad (user) ;
                    }
                    else
                    {
                        if(Com->objpoint!=NULL){
                        Com->objpoint->mygrad(user) ;
                        Com->f = Com->objpoint->myvalue(user) ;
                        }
                        if(Com->avec_objpoint!=NULL){
                        Com->avec_objpoint->mygrad(user) ;
                        Com->f = Com->avec_objpoint->myvalue(user) ;
                        }
                    }
                    /* Shrink x and g to the reduced space */
                    asa_shrinkxg (x, g, Com) ;
                }
                else
                {
                    /* Evaluate function and gradient */
                    user->ifree = NULL ;
                    user->nfree = Com->n ;
                    if ( Com->valgrad != NULL )
                    {
                      // Com->f = Com->valgrad (user) ;
                    }
                    else
                    {
                       if(Com->objpoint!=NULL){
                       Com->objpoint->mygrad(user) ;
                       Com->f = Com->objpoint->myvalue(user) ;
                       }
                       if(Com->avec_objpoint!=NULL){
                       Com->avec_objpoint->mygrad(user) ;
                       Com->f = Com->avec_objpoint->myvalue(user) ;
                       }
                    }
                }
            }
            else
            {
                asa_step (xtemp, x, d, alpha, Com->nfree) ;
                if ( Com->DimReduce )
                {
                    /* Expand xtemp to the full space*/
                    asa_expandx (xtemp, Com) ;
                    /* Evaluate function and gradient */
                    user->ifree = Com->ifree ;
                    user->nfree = Com->nfree ;
                    if ( Com->valgrad != NULL )
                    {
                      //  Com->f = Com->valgrad (user) ;
                    }
                    else
                    {
                        if(Com->objpoint!=NULL){
                        Com->objpoint->mygrad(user) ;
                        Com->f = Com->objpoint->myvalue(user) ;
                        }
                        if(Com->avec_objpoint!=NULL){
                        Com->avec_objpoint->mygrad(user) ;
                        Com->f = Com->avec_objpoint->myvalue(user) ;
                        }
                    }
                    /* Shrink xtemp and gtemp to the reduced space */
                    asa_shrinkxg (xtemp, gtemp, Com) ;
                }
                else
                {
                    /* Evaluate function and gradient */
                    user->ifree = NULL ;
                    user->nfree = Com->n ;
                    if ( Com->valgrad != NULL )
                    {
                     //  Com->f = Com->valgrad (user) ;
                    }
                    else
                    {
                       if(Com->objpoint!=NULL){
                       Com->objpoint->mygrad(user) ;
                       Com->f = Com->objpoint->myvalue(user) ;
                       }
                       if(Com->avec_objpoint!=NULL){
                       Com->avec_objpoint->mygrad(user) ;
                       Com->f = Com->avec_objpoint->myvalue(user) ;
                       }
                    }
                }
                Com->df = asa_dot (gtemp, d, Com->nfree) ;
            }
            Com->nf++ ;
            Com->ng++ ;
            if ( (Com->df != Com->df) || (Com->f != Com->f) ||
                 (Com->df == ASA_INF)     || (Com->f == ASA_INF)    ||
                 (Com->df ==-ASA_INF)     || (Com->f ==-ASA_INF) ) return (11) ;
        }
        else if ( !strcmp (what, "f") ) /* compute function */
        {
            asa_step (xtemp, x, d, alpha, Com->nfree) ;
            if ( Com->DimReduce )
            {
               /* Expand xtemp to the full space*/
               asa_expandx (xtemp, Com) ;
               /* Evaluate function */
               user->ifree = Com->ifree ;
               user->nfree = Com->nfree ;
               if(Com->objpoint!=NULL)
               Com->f = Com->objpoint->myvalue(user) ;
               if(Com->avec_objpoint!=NULL)
               Com->f = Com->avec_objpoint->myvalue(user) ;
               /* Shrink xtemp to the reduced space */
               asa_shrinkx (xtemp, Com) ;
            }
            else
            {
               /* Evaluate function */
               user->ifree = NULL ;
               user->nfree = Com->n ;
               if(Com->objpoint!=NULL)
               Com->f = Com->objpoint->myvalue(user) ;
               if(Com->avec_objpoint!=NULL)
               Com->f = Com->avec_objpoint->myvalue(user) ;
            }
            Com->nf++ ;
            if ( (Com->f != Com->f) || (Com->f == ASA_INF) || (Com->f ==-ASA_INF) )
                return (11) ;
        }
        else
        {
            asa_step (xtemp, x, d, alpha, Com->nfree) ;
            if ( Com->DimReduce )
            {
                /* Expand xtemp to the full space*/
                asa_expandx (xtemp, Com) ;
                /* Evaluate gradient */
                user->ifree = Com->ifree ;
                user->nfree = Com->nfree ;
                if(Com->objpoint!=NULL)
                Com->objpoint->mygrad(user) ;
                if(Com->avec_objpoint!=NULL)
                Com->avec_objpoint->mygrad(user) ;
                /* Shrink x and g to the reduced space */
                asa_shrinkxg (xtemp, gtemp, Com) ;
            }
            else
            {
               /* Evaluate gradient */
               user->ifree = NULL ;
               user->nfree = Com->n ;
               if(Com->objpoint!=NULL)
               Com->objpoint->mygrad(user) ;
               if(Com->avec_objpoint!=NULL)
               Com->avec_objpoint->mygrad(user) ;
            }
            Com->df = asa_dot (gtemp, d, Com->nfree) ;
            Com->ng++ ;
            if ( (Com->df != Com->df) || (Com->df == ASA_INF) || (Com->df ==-ASA_INF) )
                return (11) ;
        }
    }
    return (0) ;
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Start of routines that could use the BLAS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* =========================================================================
   ==== asa_matvec =========================================================
   =========================================================================
   Compute y = A*x or A'*x where A is a dense rectangular matrix
   ========================================================================= */
 void asa_matvec
(
    double *y, /* product vector */
    double *A, /* dense matrix */
    double *x, /* input vector */
    int     n, /* number of columns of A */
    ASA_INT     m, /* number of rows of A */
    int     w  /* T => y = A*x, F => y = A'*x */
)
{
/* if the blas have not been installed, then hand code the produce */
#ifdef NOBLAS
    ASA_INT j, l ;
    l = 0 ;
    if ( w )
    {
        asa_scale0 (y, A, x [0], (int) m) ;
        for (j = 1; j < n; j++)
        {
            l += m ;
            asa_daxpy0 (y, A+l, x [j], (int) m) ;
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            y [j] = asa_dot0 (A+l, x, (int) m) ;
            l += m ;
        }
    }
#endif

/* if the blas have been installed, then possibly call gdemv */
#ifndef NOBLAS
    ASA_INT j, l ;
    BLAS_INT M, N ;
    if ( w || (!w && (m*n < MATVEC_START)) )
    {
        l = 0 ;
        if ( w )
        {
            asa_scale (y, A, x [0], m) ;
            for (j = 1; j < n; j++)
            {
                l += m ;
                asa_daxpy (y, A+l, x [j], m) ;
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                y [j] = asa_dot0 (A+l, x, (int) m) ;
                l += m ;
            }
        }
    }
    else
    {
        M = (BLAS_INT) m ;
        N = (BLAS_INT) n ;
        /* only use transpose mult with blas
        CG_DGEMV ("n", &M, &N, one, A, &M, x, blas_one, zero, y, blas_one) ;*/
        CG_DGEMV ("t", &M, &N, one, A, &M, x, blas_one, zero, y, blas_one) ;
    }
#endif

    return ;
}

/* =========================================================================
   ==== asa_trisolve =====================================================
   =========================================================================
   Solve Rx = y or R'x = y where R is a dense upper triangular matrix
   ========================================================================= */
 void asa_trisolve
(
    double *x, /* right side on input, solution on output */
    double *R, /* dense matrix */
    int     m, /* leading dimension of R */
    int     n, /* dimension of triangular system */
    int     w  /* T => Rx = y, F => R'x = y */
)
{
    int i, l ;
    if ( w )
    {
        l = m*n ;
        for (i = n; i > 0; )
        {
            i-- ;
            l -= (m-i) ;
            x [i] /= R [l] ;
            l -= i ;
            asa_daxpy0 (x, R+l, -x [i], i) ;
        }
    }
    else
    {
        l = 0 ;
        for (i = 0; i < n; i++)
        {
            x [i] = (x [i] - asa_dot0 (x, R+l, i))/R [l+i] ;
            l += m ;
        }
    }

/* equivalent to:
    BLAS_INT M, N ;
    M = (BLAS_INT) m ;
    N = (BLAS_INT) n ;
    if ( w ) CG_DTRSV ("u", "n", "n", &N, R, &M, x, blas_one) ;
    else     CG_DTRSV ("u", "t", "n", &N, R, &M, x, blas_one) ; */

    return ;
}

/* =========================================================================
   ==== asa_scale0 =========================================================
   =========================================================================
   compute y = s*x where s is a scalar
   ========================================================================= */
 void asa_scale0
(
    double *y, /* output vector */
    double *x, /* input vector */
    double  s, /* scalar */
    int     n /* length of vector */
)
{
    int i, n5 ;
    n5 = n % 5 ;
    if ( s == -ONE)
    {
       for (i = 0; i < n5; i++) y [i] = -x [i] ;
       for (; i < n;)
       {
           y [i] = -x [i] ;
           i++ ;
           y [i] = -x [i] ;
           i++ ;
           y [i] = -x [i] ;
           i++ ;
           y [i] = -x [i] ;
           i++ ;
           y [i] = -x [i] ;
           i++ ;
       }
    }
    else
    {
        for (i = 0; i < n5; i++) y [i] = s*x [i] ;
        for (; i < n;)
        {
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
        }
    }
    return ;
}

/* =========================================================================
   ==== asa_scale ==========================================================
   =========================================================================
   compute y = s*x where s is a scalar
   ========================================================================= */
 void asa_scale
(
    double *y, /* output vector */
    double *x, /* input vector */
    double  s, /* scalar */
    ASA_INT     n /* length of vector */
)
{
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    if ( y == x)
    {
#ifdef NOBLAS
        for (i = 0; i < n5; i++) y [i] *= s ;
        for (; i < n;)
        {
            y [i] *= s ;
            i++ ;
            y [i] *= s ;
            i++ ;
            y [i] *= s ;
            i++ ;
            y [i] *= s ;
            i++ ;
            y [i] *= s ;
            i++ ;
        }
#endif
#ifndef NOBLAS
        if ( n < DSCAL_START )
        {
            for (i = 0; i < n5; i++) y [i] *= s ;
            for (; i < n;)
            {
                y [i] *= s ;
                i++ ;
                y [i] *= s ;
                i++ ;
                y [i] *= s ;
                i++ ;
                y [i] *= s ;
                i++ ;
                y [i] *= s ;
                i++ ;
            }
        }
        else
        {
            BLAS_INT N ;
            N = (BLAS_INT) n ;
            CG_DSCAL (&N, &s, x, blas_one) ;
        }
#endif
    }
    else
    {
        for (i = 0; i < n5; i++) y [i] = s*x [i] ;
        for (; i < n;)
        {
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
            y [i] = s*x [i] ;
            i++ ;
        }
    }
    return ;
}

/* =========================================================================
   ==== asa_daxpy0 ==========================================================
   =========================================================================
   Compute x = x + alpha d
   ========================================================================= */
 void asa_daxpy0
(
    double     *x, /* input and output vector */
    double     *d, /* direction */
    double  alpha, /* stepsize */
    int         n  /* length of the vectors */
)
{
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    if (alpha == -ONE)
    {
        for (i = 0; i < n5; i++) x [i] -= d[i] ;
        for (; i < n; i += 5)
        {
            x [i]   -= d [i] ;
            x [i+1] -= d [i+1] ;
            x [i+2] -= d [i+2] ;
            x [i+3] -= d [i+3] ;
            x [i+4] -= d [i+4] ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++) x [i] += alpha*d[i] ;
        for (; i < n; i += 5)
        {
            x [i]   += alpha*d [i] ;
            x [i+1] += alpha*d [i+1] ;
            x [i+2] += alpha*d [i+2] ;
            x [i+3] += alpha*d [i+3] ;
            x [i+4] += alpha*d [i+4] ;
        }
    }
    return ;
}

/* =========================================================================
   ==== asa_daxpy ===========================================================
   =========================================================================
   Compute x = x + alpha d
   ========================================================================= */
 void asa_daxpy
(
    double     *x, /* input and output vector */
    double     *d, /* direction */
    double  alpha, /* stepsize */
    ASA_INT         n  /* length of the vectors */
)
{
#ifdef NOBLAS
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    if (alpha == -ONE)
    {
        for (i = 0; i < n5; i++) x [i] -= d[i] ;
        for (; i < n; i += 5)
        {
            x [i]   -= d [i] ;
            x [i+1] -= d [i+1] ;
            x [i+2] -= d [i+2] ;
            x [i+3] -= d [i+3] ;
            x [i+4] -= d [i+4] ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++) x [i] += alpha*d[i] ;
        for (; i < n; i += 5)
        {
            x [i]   += alpha*d [i] ;
            x [i+1] += alpha*d [i+1] ;
            x [i+2] += alpha*d [i+2] ;
            x [i+3] += alpha*d [i+3] ;
            x [i+4] += alpha*d [i+4] ;
        }
    }
#endif

#ifndef NOBLAS
    ASA_INT i, n5 ;
    BLAS_INT N ;
    if ( n < DAXPY_START )
    {
        n5 = n % 5 ;
        if (alpha == -ONE)
        {
            for (i = 0; i < n5; i++) x [i] -= d[i] ;
            for (; i < n; i += 5)
            {
                x [i]   -= d [i] ;
                x [i+1] -= d [i+1] ;
                x [i+2] -= d [i+2] ;
                x [i+3] -= d [i+3] ;
                x [i+4] -= d [i+4] ;
            }
        }
        else
        {
            for (i = 0; i < n5; i++) x [i] += alpha*d[i] ;
            for (; i < n; i += 5)
            {
                x [i]   += alpha*d [i] ;
                x [i+1] += alpha*d [i+1] ;
                x [i+2] += alpha*d [i+2] ;
                x [i+3] += alpha*d [i+3] ;
                x [i+4] += alpha*d [i+4] ;
            }
        }
    }
    else
    {
        N = (BLAS_INT) n ;
        CG_DAXPY (&N, &alpha, d, blas_one, x, blas_one) ;
    }
#endif

    return ;
}

/* =========================================================================
   ==== asa_dot0 ===========================================================
   =========================================================================
   Compute dot product of x and y, vectors of length n
   ========================================================================= */
 double asa_dot0
(
    double *x, /* first vector */
    double *y, /* second vector */
    int     n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double t ;
    t = ZERO ;
    if ( n <= 0 ) return (t) ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) t += x [i]*y [i] ;
    for (; i < n; i += 5)
    {
        t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
                        + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
    }
    return (t) ;
}

/* =========================================================================
   ==== asa_dot ============================================================
   =========================================================================
   Compute dot product of x and y, vectors of length n
   ========================================================================= */
double asa_dot
(
    double *x, /* first vector */
    double *y, /* second vector */
    ASA_INT     n /* length of vectors */
)
{
#ifdef NOBLAS
    ASA_INT i, n5 ;
    double t ;
    t = ZERO ;
    if ( n <= 0 ) return (t) ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) t += x [i]*y [i] ;
    for (; i < n; i += 5)
    {
        t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
                        + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
    }
    return (t) ;
#endif

#ifndef NOBLAS
    ASA_INT i, n5 ;
    double t ;
    BLAS_INT N ;
    if ( n < DDOT_START )
    {
        t = ZERO ;
        if ( n <= 0 ) return (t) ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) t += x [i]*y [i] ;
        for (; i < n; i += 5)
        {
            t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
                            + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
        }
        return (t) ;
    }
    else
    {
        N = (BLAS_INT) n ;
        return (CG_DDOT (&N, x, blas_one, y, blas_one)) ;
    }
#endif
}

/* =========================================================================
   === asa_copy0 ============================================================
   =========================================================================
   Copy vector x into vector y
   ========================================================================= */
 void asa_copy0
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    int     n  /* length of vectors */
)
{
    int i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) y [i] = x [i] ;
    for (; i < n; )
    {
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
    }

    return ;
}

/* =========================================================================
   === asa_copy =============================================================
   =========================================================================
   Copy vector x into vector y
   ========================================================================= */
void asa_copy
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    ASA_INT     n  /* length of vectors */
)
{
#ifdef NOBLAS
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) y [i] = x [i] ;
    for (; i < n; )
    {
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
        y [i] = x [i] ;
        i++ ;
    }
#endif

#ifndef NOBLAS
    ASA_INT i, n5 ;
    BLAS_INT N ;
    if ( n < DCOPY_START )
    {
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) y [i] = x [i] ;
        for (; i < n; )
        {
            y [i] = x [i] ;
            i++ ;
            y [i] = x [i] ;
            i++ ;
            y [i] = x [i] ;
            i++ ;
            y [i] = x [i] ;
            i++ ;
            y [i] = x [i] ;
            i++ ;
        }
    }
    else
    {
        N = (BLAS_INT) n ;
        CG_DCOPY (&N, x, blas_one, y, blas_one) ;
    }
#endif

    return ;
}

/* =========================================================================
   === asa_update_xy =======================================================
   =========================================================================
   Set x = xnew and y = ynew
   ========================================================================= */
 void asa_update_xy
(
    double *x, /* output of copy */
    double *xnew, /* input of copy */
    double *y, /* output of copy */
    double *ynew, /* input of copy */
    ASA_INT     n  /* length of vectors */
)
{
#ifdef NOBLAS
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++)
    {
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
    }
    for (; i < n; )
    {
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
        i++ ;
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
        i++ ;
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
        i++ ;
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
        i++ ;
        x [i] = xnew [i] ;
        y [i] = ynew [i] ;
        i++ ;
    }
#endif

#ifndef NOBLAS
    ASA_INT i, n5 ;
    BLAS_INT N ;
    if ( n < DCOPY_START )
    {
        n5 = n % 5 ;
        for (i = 0; i < n5; i++)
        {
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
        }
        for (; i < n; )
        {
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
            i++ ;
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
            i++ ;
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
            i++ ;
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
            i++ ;
            x [i] = xnew [i] ;
            y [i] = ynew [i] ;
            i++ ;
        }
    }
    else
    {
        N = (BLAS_INT) n ;
        CG_DCOPY (&N, xnew, blas_one, x, blas_one) ;
        CG_DCOPY (&N, ynew, blas_one, y, blas_one) ;
    }
#endif

    return ;
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       End of routines that could use the BLAS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* =========================================================================
   ==== asa_update_pg ======================================================
   =========================================================================
   compute ginorm and pgnorm of g
   ========================================================================= */
 double asa_update_pg
(
    double      *x, /* x */
    double      *g, /* gradient */
    double *Ginorm, /* ginorm of g*/
    double   bdist, /* distance of x to the boundary is greater than bdist*/
    double     *lo, /* lower bound */
    double     *hi, /* upper bound */
    ASA_INT      nfree, /* dimension of subspace */
    ASA_INT         n, /* length of vectors */
    asa_com   *Com
)
{
    ASA_INT j ;
    double ginorm, pgnorm ;
    double xj, gj, xg, xp, t ;

    ginorm = 0. ;
    for (j = 0; j < nfree; j++)
    {
         xj = x [j] ;
         gj = g [j] ;
         t = fabs (gj) ;
         if ( bdist < t ) break ; /* check for active constraint */
         ginorm = MAX2 (ginorm, t) ;
    }
    pgnorm = ginorm ;
    if ( j < nfree )
    {
         ginorm = t ;
         xg = xj - gj ;
         if      ( xg > hi [j] ) xp = hi [j] - xj ;
         else if ( xg < lo [j] ) xp = xj - lo [j] ;
         else                    xp = t ;
         pgnorm = MAX2 (pgnorm, xp) ;
         for (j++; j < nfree; j++)
         {
             xj = x [j] ;
             gj = g [j] ;
             ginorm = MAX2 (ginorm, fabs (gj)) ;
             xg = xj - gj ;
             if      ( xg > hi [j] ) xp = hi [j] - xj ;
             else if ( xg < lo [j] ) xp = xj - lo [j] ;
             else                    xp = fabs (gj) ;
             pgnorm = MAX2 (pgnorm, xp) ;
         }
    }
    for (; j < n; j++)
    {
         xj = Com->x [j] ;
         gj = g [j] ;
         xg = xj - gj ;
         if      ( xg > hi [j] ) xp = hi [j] - xj ;
         else if ( xg < lo [j] ) xp = xj - lo [j] ;
         else                    xp = fabs (gj) ;
         pgnorm = MAX2 (pgnorm, xp) ;
    }
    *Ginorm = ginorm ;
    return (pgnorm) ;
}

/* =========================================================================
   ==== asa_update_2 =======================================================
   =========================================================================
   Set gold = gnew (if not equal), compute 2-norm^2 of gnew, and optionally
      set d = -gnew
   ========================================================================= */
 double asa_update_2
(
    double *gold, /* old g */
    double *gnew, /* new g */
    double    *d, /* d */
    ASA_INT        n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double s, t ;
    t = ZERO ;
    n5 = n % 5 ;

    if ( d == NULL )
    {
        for (i = 0; i < n5; i++)
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
        }
        for (; i < n; )
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            i++ ;
        }
    }
    else if ( gold != NULL )
    {
        for (i = 0; i < n5; i++)
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
        }
        for (; i < n; )
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++)
        {
            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
        }
        for (; i < n; )
        {
            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;
        }
    }
    return (t) ;
}
/* =========================================================================
   ==== asa_update_dg0 ======================================================
   =========================================================================
   Set gold = gnew, d = -gnew, optionally compute 2-norm of gnew
   ========================================================================= */
 void asa_update_dg0
(
    double   *gold, /* old g */
    double   *gnew, /* new g */
    double      *d, /* d */
    double *gnorm2, /* 2-norm of g */
    ASA_INT          n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double s, t ;
    s = ZERO ;
    n5 = n % 5 ;

    if ( gnorm2 == NULL )
    {
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
        }
        for (; i < n; )
        {
            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
        }
        for (; i < n; )
        {
            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;

            t = gnew [i] ;
            s += t*t ;
            gold [i] = t ;
            d [i] = -t ;
            i++ ;
        }
        *gnorm2 = s ;
    }

}

/* =========================================================================
   ==== asa_update_dg ======================================================
   =========================================================================
   Set d = -gnew + beta*d, gold = gnew, compute 2-norm of d,
   and optionally the 2-norm of gnew
   ========================================================================= */
 double asa_update_dg
(
    double   *gold,
    double   *gnew,
    double      *d,
    double    beta,
    double *gnorm2, /* 2-norm of gnew */
    ASA_INT          n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double dnorm2, s, t ;
    s = ZERO ;
    dnorm2 = ZERO ;
    n5 = n % 5 ;

    if ( gnorm2 == NULL )
    {
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;
        }
    }
    else
    {
        s = ZERO ;
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;

            t = gnew [i] ;
            gold [i] = t ;
            s += t*t ;
            t = -t + beta*d [i] ;
            d [i] = t ;
            dnorm2 += t*t ;
            i++ ;
        }
        *gnorm2 = s ;
    }

    return (dnorm2) ;
}


/* =========================================================================
   ==== asa_update_ykyk ====================================================
   =========================================================================
     compute    ykyk = 2-norm(gnew-gold)^2
                ykgk = (gnew-gold) dot gnew
   ========================================================================= */
void asa_update_ykyk
(
    double *gold, /* old g */
    double *gnew, /* new g */
    double *Ykyk,
    double *Ykgk,
    ASA_INT        n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double t, yk, ykyk, ykgk ;
    ykyk = ZERO ;
    ykgk = ZERO ;
    n5 = n % 5 ;

    for (i = 0; i < n5; i++)
    {
        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
    }
    for (; i < n; )
    {
        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
        i++ ;

        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
        i++ ;

        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
        i++ ;

        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
        i++ ;

        t = gnew [i] ;
        yk = t - gold [i] ;
        ykgk += yk*t ;
        ykyk += yk*yk ;
        i++ ;
    }
    *Ykyk = ykyk ;
    *Ykgk = ykgk ;
}

/* =========================================================================
   ==== asa_update_inf2 ====================================================
   =========================================================================
   Set gold = gnew, compute inf-norm of gnew & 2-norm of gnew, set d = -gnew
   ========================================================================= */
 double asa_update_inf2
(
    double   *gold, /* old g */
    double   *gnew, /* new g */
    double      *d, /* d */
    double *gnorm2, /* 2-norm of g */
    ASA_INT          n /* length of vectors */
)
{
    ASA_INT i, n5 ;
    double gnorm, s, t ;
    gnorm = ZERO ;
    s = ZERO ;
    n5 = n % 5 ;

    for (i = 0; i < n5; i++)
    {
        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
    }
    for (; i < n; )
    {
        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;
    }
    *gnorm2 = s ;
    return (gnorm) ;
}

/* =========================================================================
   === asa_step ============================================================
   =========================================================================
   Compute xtemp = x + alpha d
   ========================================================================= */
 void asa_step
(
    double *xtemp , /*output vector */
    double     *x , /* initial vector */
    double     *d , /* search direction */
    double  alpha , /* stepsize */
    ASA_INT         n   /* length of the vectors */
)
{
    ASA_INT n5, i ;
    n5 = n % 5 ;
    if (alpha == -ONE)
    {
        for (i = 0; i < n5; i++) xtemp [i] = x[i] - d[i] ;
        for (; i < n; i += 5)
        {
            xtemp [i]   = x [i]   - d [i] ;
            xtemp [i+1] = x [i+1] - d [i+1] ;
            xtemp [i+2] = x [i+2] - d [i+2] ;
            xtemp [i+3] = x [i+3] - d [i+3] ;
            xtemp [i+4] = x [i+4] - d [i+4] ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++) xtemp [i] = x[i] + alpha*d[i] ;
        for (; i < n; i += 5)
        {
            xtemp [i]   = x [i]   + alpha*d [i] ;
            xtemp [i+1] = x [i+1] + alpha*d [i+1] ;
            xtemp [i+2] = x [i+2] + alpha*d [i+2] ;
            xtemp [i+3] = x [i+3] + alpha*d [i+3] ;
            xtemp [i+4] = x [i+4] + alpha*d [i+4] ;
        }
    }
    return ;
}

/* =========================================================================
   ==== asa_init ============================================================
   =========================================================================
   initialize x to a given scalar value
   ========================================================================= */
 void asa_init
(
    double *x, /* input and output vector */
    double  s, /* scalar */
    ASA_INT     n /* length of vector */
)
{
    ASA_INT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = s ;
    for (; i < n;)
    {
        x [i] = s ;
        i++ ;
        x [i] = s ;
        i++ ;
        x [i] = s ;
        i++ ;
        x [i] = s ;
        i++ ;
        x [i] = s ;
        i++ ;
    }
    return ;
}

/* =========================================================================
   === asa_max =============================================================
   =========================================================================
   Return max {fabs (x [j]) : 1 <= j < n}
   ========================================================================= */
double asa_max
(
    double *x,
    ASA_INT     n
)
{
    double xnorm ;
    ASA_INT j, n5 ;
    n5 = n % 5 ;
    xnorm = ZERO ;
    for (j = 0; j < n5; j++) if ( xnorm < fabs (x [j]) ) xnorm = fabs (x [j]) ;
    for (; j < n; j += 5)
    {
        if ( xnorm < fabs (x [j]  ) ) xnorm = fabs (x [j]) ;
        if ( xnorm < fabs (x [j+1]) ) xnorm = fabs (x [j+1]) ;
        if ( xnorm < fabs (x [j+2]) ) xnorm = fabs (x [j+2]) ;
        if ( xnorm < fabs (x [j+3]) ) xnorm = fabs (x [j+3]) ;
        if ( xnorm < fabs (x [j+4]) ) xnorm = fabs (x [j+4]) ;
    }
    return (xnorm) ;
}


/* =========================================================================
   === asa_project =========================================================
   =========================================================================
   Project a vector into the feasible set
   ========================================================================= */
void asa_project
(
    double  *xnew,
    double     *x,
    double     *d,
    double  alpha,
    asa_com  *Com   /* cg com structure */
)
{
    ASA_INT j, n ;
    double t, *lo, *hi ;
    lo = Com->lo ;
    hi = Com->hi ;
    n = Com->nfree ;
    for (j = 0; j < n; j++)
    {
        t = x [j] + alpha*d [j] ;
        if      ( t > hi [j] ) t = hi [j] ;
        else if ( t < lo [j] ) t = lo [j] ;
        xnew [j] = t ;
    }
}

/* =========================================================================
   === asa_maxstep =========================================================
   =========================================================================
   Compute maximum step in the search direction until hitting boundary
   ========================================================================= */
 void asa_maxstep
(
    double       *x, /* current iterate */
    double       *d, /* direction */
    asa_com    *Com
)
{
    double bdist, step, minstep, maxstep, xj, t, *lo, *hi ;
    ASA_INT j, n ;

    Com->minflag = TRUE ;
    n = Com->nfree ;
    minstep = ASA_INF ;
    maxstep = ZERO ;
    bdist = ASA_INF ;
    lo = Com->lo ;
    hi = Com->hi ;

    for (j = 0;  j < n; j++)
    {
        xj = x [j] ;
        if ( d [j] > ZERO )
        {
            if ( hi [j] < ASA_INF )
            {
                t = hi [j] - xj ;
                if ( bdist > t ) bdist = t ;
                step = t/d [j] ;
                minstep = MIN2 (minstep, step) ;
                maxstep = MAX2 (maxstep, step) ;
            }
            t = xj - lo [j] ;
            if ( bdist > t ) bdist = t ;
        }
        else if ( d [j] < ZERO )
        {
            if ( lo [j] >-ASA_INF )
            {
                t = xj - lo [j] ;
                if ( bdist > t ) bdist = t ;
                step = -t/d [j] ;
                minstep = MIN2 (minstep, step) ;
                maxstep = MAX2 (maxstep, step) ;
            }
            t = hi [j] - xj ;
            if ( bdist > t ) bdist = t ;
        }
        else
        {
            t = xj - lo [j] ;
            if ( bdist > t ) bdist = t ;
            t = hi [j] - xj ;
            if ( bdist > t ) bdist = t ;
        }
    }
    Com->bdist = bdist ;
    Com->minstep = minstep ;
    Com->maxstep = maxstep ;
    return ;
}

/* =========================================================================
   === asa_grad_proj =======================================================
   =========================================================================
   Nonmonotone gradient projection algorithm modified to combine with cg_descent
   for bound constraints problem.

   Notation:
   fmin  = min {f(x_i),  0 <= i <= k}, k is current iteration (best value)
   fmax  = max {f_{k-i}, i = 0, 1, ... , min(k,m-1)}, m = memory
   fr    = reference function value
   f     = f(x_k), current function value
   fc    = maximum objective function value since the last minimum
          (best) function value was found. In other words, if
          k is the current iteration f(k_1) = fmin, then
          fc = max {f(x_i), k_1 <= i <= k}
   fcomp = min {fr, fmax}
   fr_pert = fr + pert, pert = Parm->eps*|fcomp| or Parm->eps
                        depending on PertRule
   fcomp_pert = fcomp + pert
   ftemp = f(xtemp), temporary (or a trial) function value at xtemp
   ========================================================================= */
 int asa_grad_proj /*return:
                      -1 (give active constraints to cg routine)
                       0 (convergence tolerance satisfied)
                      14 (number of iterations or function evaluations
                          exceed limit)
                      15 (line search fails)
                      16 (search direction in linesearch is not descent)
                      17 (function value became nan) */
(
    asa_com *Com
)
{
    int status, hitbound, getbound, freebound ;
    ASA_INT count, i, ident, index, iter, j, ll, ll0, mcount, mm, n, nf, nf_line,
        ng, nl, np ;
    double alpha, armijo_decay, armijo0, armijo1, f, fmin, fc, fr, sts, gtd,
           fmax, lambda, pgnorm, ginorm, gnorm, xnorm,
           xj, gj, xp, xg, t, th, tl,
           pert_lo, pert_hi, atemp, ftemp, fcomp, sty, yty, s, y, cosine,
           *lo, *hi, *x, *d, *g, *xtemp, *gtemp, *pg, *lastfvalues ;
    double fr_pert, fcomp_pert, dphia, Armijo_hi, AArmijo_hi ;
    asa_parm *Parm ;

    n = Com->n ;
    x = Com->x ;
    lo = Com->lo ;
    hi = Com->hi ;
    d = Com->d ;
    g = Com->g ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    pg = Com->pg ;
    lastfvalues = Com->lastfvalues ;
    pgnorm = Com->pgnorm ;
    ginorm = Com->ginorm ;
    nf = Com->nf ;
    ng = Com->ng ;
    Parm = Com->asaParm ;
    pert_lo = Parm->pert_lo ;
    pert_hi = Parm->pert_hi ;
    armijo_decay = Parm->armijo_decay ;
    armijo0 = Parm->armijo0 ;
    armijo1 = Parm->armijo1 ;
    f = Com->f ;

    iter = 0 ;
    status = 0 ;
    count = 0 ;
    mcount = 2 ;
    ident = FALSE ;
    lambda = Com->alpha ;
    lastfvalues [0] = f ;
    for (i = 1; i < Parm->m; i++) lastfvalues [i] = -ASA_INF ;
    mm = 0 ; /* number of iterations in the current CBB cycle */
    ll = 0 ; /* zero as long as unit steps in line search */
    nl = 0 ; /* number of iterations since fmin decreased in value */
    np = 0 ; /* number of times initial stepsize was accepted in line search */
    fmin = f ;
    fr = f ;
    fc = f ;

    if ( Parm->PrintLevel >= 3 )
    {
        printf ("Initial stepsize in cbb: %14.6e\n", lambda) ;
    }

    while ( TRUE )
    {
        if ( Parm->PrintLevel >= 2 )
        {
            printf ("cbb iter: %5i f: %14.6e pgnorm: %14.6e\n\n",
                     (int) iter, f, pgnorm) ;
        }

        if ( !Parm->GradProjOnly )
        {
            if ( ginorm >= Com->tau1*pgnorm )
            {
               if ( ident || (count >= mcount) )
               {
                   status = -1 ;
                   goto Exit ;
               }
            }
            else
            {
                if ( ident )
                {
                    ident = FALSE ;
                    Com->tau1 *= Parm->tau1_decay ;
                    Com->tau2 *= Parm->tau2_decay ;
                }
            }
        }
        iter++ ;
        hitbound = FALSE ;
        getbound = FALSE ;
        freebound = FALSE ;
        sts = ZERO ;
        gtd = ZERO ;
        for (j = 0; j < n; j++)
        {
            xj = x [j] ;
            gj = g [j] ;
            xp = -lambda*gj ;
            xg = xj + xp ;
            th = hi [j] ;
            tl = lo [j] ;
            if ( xg >= th )
            {
                xp = th - xj ;
                xtemp [j] = th ;
                if ( xp > pert_hi ) getbound = TRUE ;
            }
            else if ( xg <= tl )
            {
                xp = tl - xj ;
                xtemp [j] = tl ;
                if ( -xp > pert_lo ) getbound = TRUE ;
            }
            else
            {
                xtemp [j] = xg ;
                if ( (xj == th) || (xj == tl) ) freebound = TRUE ;
            }
            d [j] = xp ;
            gtd += gj*xp ; /* g'd (derivative in search direction) */
            sts += xp*xp ;
        }
        if ( getbound ) ll++ ;
        nf_line = Com->nf ;

        if (gtd >= ZERO)
        {
            status = 13 ;
            Com->gtd = gtd ;
            goto Exit ;
        }

       /* start of cbb line search */
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("Linesearch in cbb, f: %14.6e gtd: %14.6e\n", f, gtd) ;
        }
        fmax = lastfvalues [0] ;
        for (i = 1; i < Parm->m; i++) fmax = MAX2 (fmax, lastfvalues [i]) ;
        alpha = ONE ;
        ftemp = asa_f (xtemp, Com) ;

        if ( nl == Parm->L )
        {
            fr = fmax ;
            t = (fr-fmin)/(fc-fmin) ;
            if ( t > Parm->gamma1 ) fr = fc ;
            nl = 0 ;
        }

        if ( (np > Parm->P) && (fmax > f) )
        {
           t = (fr-f)/(fmax-f) ;
           if ( t > Parm->gamma2 ) fr = fmax ;
        }

        fcomp = MIN2 (fmax, fr) ;

        if ( Parm->PrintLevel >= 4 )
        {
            printf ("fr: %14.6e fcomp: %14.6e\n", fr, fcomp) ;
        }

        /* Approximate nonmonotone Armijo line search, decrease alpha until:
           phi'(alpha) <= [2(phi_r - phi(0))/alpha] + (2 delta - 1) phi'(0) and
           phi(alpha) <= phi_r, where phi_r = fr_pert or fcomp_pert. */
        if ( Com->AArmijo)
        {
            if ( Parm->PertRule ) t = Parm->eps*fabs(fcomp) ;
            else                  t = Parm->eps ;
            fr_pert = fr + t ;
            fr = fr_pert ;
            fcomp_pert = fcomp + t ;
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Perform approximate Armijo line search\n") ;
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("fr_pert: %14.6e fcomp_pert: %14.6e\n",
                             fr_pert, fcomp_pert) ;
                }
            }

            AArmijo_hi = (TWO*Parm->delta - ONE)*gtd ;
            if ( ftemp != ftemp ) /* function value is nan, reduce stepsize */
            {
                for (i = 0; i < Parm->nshrink; i++)
                {
                    ll++ ;
                    alpha *= Parm->nan_fac ;
                    asa_step (xtemp, x, d, alpha, n) ;
                    ftemp = asa_f (xtemp, Com) ;
                    if ( ftemp == ftemp ) break ;
                }
                if ( (i == Parm->nshrink) || (alpha == ZERO) )
                {
                    status = 14 ;
                    goto exit_with_error ;
                }

                if ( ftemp <= fcomp_pert)
                {
                    asa_g (gtemp, xtemp, Com) ;
                    dphia = asa_dot (gtemp, d, n) ;
                    if (dphia <= TWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ; /* unit step is valid */
                }
            }
            else
            {
                if ( mm == 0 )
                {
                    if ( ftemp <= fr_pert )
                    {
                        asa_g (gtemp, xtemp, Com) ;
                        dphia = asa_dot (gtemp, d, n) ;
                        if (dphia <= TWO*(fr_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
                else
                {
                    if ( ftemp <= fcomp_pert )
                    {
                        asa_g (gtemp, xtemp, Com) ;
                        dphia = asa_dot (gtemp, d, n) ;
                        if (dphia <= TWO*(fcomp_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
            }

            /* backtracking approximate nonmonotone line search */
            ll0 = ll ;
            while ( TRUE )
            {
                /* Modified Raydan's quadratic interpolation line search */
                t = TWO*(ftemp-f-alpha*gtd) ;
                if ( t != ZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                asa_step (xtemp, x, d, alpha, n) ; /* xtemp = x + alpha*d */
                ftemp = asa_f (xtemp, Com) ;
                ll++ ;

                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e ftemp: %14.6e\n", alpha, ftemp) ;
                }

                if ( ftemp <= fcomp_pert )
                {
                    asa_g (gtemp, xtemp, Com) ;
                    dphia = asa_dot (gtemp, d, n) ;
                    if (dphia <= TWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ;
                }

                if ( (alpha <= ZERO) || (ll-ll0 >= Parm->max_backsteps) )
                {
                   status = 15 ;
                   goto exit_with_error ;
                }
            }
            /* End of approximate Armijo line search */
        }

        /* Ordinary nonmonotone Armijo line search, decrease alpha until
           phi(alpha) <= phi_r + alpha * delta * phi'(0)
           where phi_r = fr or fcomp. */
        else
        {
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Perform ordinary Armijo line search\n") ;
            }

            Armijo_hi = Parm->delta*gtd ;
            if ( ftemp != ftemp ) /* function value is nan, reduce stepsize */
            {
                for (i = 0; i < Parm->nshrink; i++)
                {
                    ll++ ;
                    alpha *= Parm->nan_fac ;
                    asa_step (xtemp, x, d, alpha, n) ;
                    ftemp = asa_f (xtemp, Com) ;
                    if ( ftemp == ftemp ) break ;
                }
                if ( (i == Parm->nshrink) || (alpha == ZERO) )
                {
                    status = 17 ;
                    goto exit_with_error ;
                }
                if ( ftemp <= fcomp+alpha*Armijo_hi ) goto exit_cbbls ;
            }
            else
            {
                if ( mm == 0 ) t = fr ;
                else           t = fcomp ;
                if ( ftemp <= t+Armijo_hi )
                {
                    mm++ ;
                    goto exit_cbbls ;
                }
            }

            ll0 = ll ;
            while ( TRUE )
            {
                /* Modified Raydan's quadratic interpolation line search */
                t = TWO*(ftemp-f-alpha*gtd) ;
                if ( t != ZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                asa_step (xtemp, x, d, alpha, n) ; /* xtemp = x + alpha*d */
                ftemp = asa_f (xtemp, Com) ;
                ll++ ;

                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e ftemp: %14.6e\n", alpha, ftemp) ;
                }

                if ( ftemp <= fcomp+alpha*Armijo_hi ) break ;

                if ( (alpha <= ZERO) || (ll-ll0 >= Parm->max_backsteps) )
                {
                    /* try approximate Armijo line search  */
                    if ( Parm->AArmijoFac > ZERO ) fr = fcomp ;
                    else /* line search fails */
                    {
                        status = 15 ;
                        goto exit_with_error ;
                    }
                }
            }
            /* End of ordinary Armijo line search */
        }

        exit_cbbls:

        if ( ftemp <= fmin )
        {
             fmin = ftemp ;
             fc = ftemp ;
             nl = 0 ;
        }
        else nl++ ;
        if ( ftemp > fc ) fc = ftemp ;

        exit_with_error:
        /* end of cbbls */

        if ( getbound && (alpha == ONE) ) hitbound = TRUE ;
        if ( Parm->PrintLevel >= 3 )
        {
            printf ("hitbound = %i freebound = %i alpha = %14.6e\n",
                     hitbound, freebound, alpha) ;
        }

        if ( hitbound || freebound ) count = 0 ;
        else                         count++ ;

        sts *= alpha*alpha ;
        if ( Com->nf == nf_line + 1 ) np++ ;
        else                          np = 0 ;

        if ( !Com->AArmijo)  asa_g (gtemp, xtemp, Com) ;

        /* linesearch fails */
        if ( status > 0 )
        {
            if ( ftemp < f )
            {
                f = ftemp ;
                asa_copy(x, xtemp, n) ;
                asa_copy(g, gtemp, n) ;
            }
            pgnorm = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = x [j] ;
                gj = g [j] ;
                xg = xj - gj ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                else                     xp = -gj ;
                pgnorm = MAX2 (pgnorm, fabs (xp)) ;
                pg [j] = xp ;
            }
            goto Exit ;
        }

        index = 0 ;

        if ( (ll >= 1) || (mm >= Parm->nm) || (iter <= 1) )
        {
            index = 1 ;
            sty = ZERO ;
            pgnorm = ZERO ;
            ginorm = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                if ( (xj - lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                    ginorm = MAX2 (ginorm, fabs (gj));
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = xj - lo [j] ;
                else                     xp = fabs (gj) ;
                pgnorm = MAX2 (pgnorm, xp) ;
                sty += (xj - x [j])*(gj - g [j]) ;
                x [j] = xj ;
                g [j] = gj ;
            }

            if ( asa_tol (pgnorm, Com) )
            {
                f = ftemp ;
                for (j  = 0; j < n; j++)
                {
                    xj = xtemp [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                    else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                    else                     xp = -gj ;
                    pg [j] = xp ;
                }
                status = 0 ;

                goto Exit ;
            }
        }

        else
        {
            pgnorm = ZERO ;
            ginorm = ZERO ;
            gnorm = ZERO ;
            sty = ZERO ;
            yty = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                t = fabs (gj) ;
                gnorm = MAX2 (gnorm, t) ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = xj - lo [j] ;
                else                     xp = t ;
                pgnorm = MAX2 (pgnorm, xp) ;
                if ( (xj - lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                {
                    ginorm = MAX2 (ginorm, t) ;
                }
                s = xj - x [j] ;
                y = gj - g [j] ;
                sty += s*y ;
                yty += y*y ;
                x [j] = xj ;
                g [j] = gj ;
            }
            if ( asa_tol (pgnorm, Com) )
            {
                f = ftemp ;
                for (j  = 0; j < n; j++)
                {
                    xj = xtemp [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                    else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                    else                     xp = -gj ;
                    pg [j] = xp ;
                }
                status = 0 ;
                goto Exit ;
            }
            s = Parm->parm3*fabs (ftemp)/gnorm ;
            t = MAX2 (s, ONE) ;
            if ( sts > t*t ) index = 1 ;
            else
            {
                t = MIN2 (s, ONE) ;
                if ( sts <= t*t )
                {
                    cosine = fabs (sty)/sqrt (sts*yty) ;
                    if ( cosine >= Parm->gamma ) index = 1 ;
                }
            }
        }

        if ( index == 1 )
        {
            ll = 0 ;
            if ( sty <= ZERO)
            {
                if ( mm >= Parm->parm4 )
                {
                    xnorm = asa_max (x, n) ;
                    t = MIN2 (ONE/pgnorm, xnorm/pgnorm) ;
                    lambda = MAX2 (t, lambda) ;
                    mm = 0 ;
                }
            }
            else
            {
                t = MAX2 (Parm->lmin, sts/sty) ;
                lambda = MIN2 (Parm->lmax, t) ;
                mm = 0 ;
            }
        }

        /* If not GradProjOnly, check if the active constraints are identified*/
        if ( !Parm->GradProjOnly &&
              pgnorm < Parm->pgdecay*MAX2 (ONE, Com->pgnorm_start) )
        {
            ident = asa_identify(x, g, pgnorm, Com) ;
        }

        f = ftemp ;
        lastfvalues [iter % Parm->m] = f ;

        /* check for excessive iterations/function evaluations */
        if ( (iter >= Com->pgmaxit) || (Com->nf - nf  >= Com->pgmaxfunc) )
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] ;
                gj = g [j] ;
                xg = xj - gj ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                else                     xp = -gj ;
                pg [j] = xp ;
            }
            status = 14 ;
            goto Exit ;
        }

        if ( !Com->AArmijo )
        {
            if ( fabs(fr - f) <= Parm->AArmijoFac*fabs(fcomp) )
                Com->AArmijo = TRUE ;
        }
    }
    Exit:
    Com->f = f ;
    Com->ginorm = ginorm ;
    Com->pgnorm = pgnorm ;
    Com->cbbiter += iter ;
    Com->cbbfunc += Com->nf - nf ;
    Com->cbbgrad += Com->ng - ng ;
    if ( Parm->PrintLevel >= 2 )
    {
        if(status != -1) printf ("cbb iter: %5i f: %14.6e pgnorm: %14.6e\n\n",
                                  (int) iter, f, pgnorm) ;
    }
    if ( Parm->PrintLevel >= 1 )
    {
        printf ("\nCBB Termination status: %i\n", status) ;
        if ( status == -1 )
            printf ("terminate cbb iteration, branch to cg iteration\n") ;

        printf ("proj gradient max norm: %13.6e\n", pgnorm) ;
        printf ("function value:         %13.6e\n", f) ;
        printf ("cbb iterations:         %13.6e\n", (double) iter) ;
        printf ("function evaluations:   %13.6e\n", (double) Com->nf - nf) ;
        printf ("gradient evaluations:   %13.6e\n", (double) Com->ng - ng) ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_init_bbstep =====================================================
   =========================================================================
   Calculate initial BB stepsize
   ========================================================================= */
 double asa_init_bbstep
(
    asa_com *Com
)
{
    int n ;
    double alpha, lmax, lmin, pgnorm, xnorm, sts, sty, t, *x ;
    x = Com->x ;
    sts = Com->sts ;
    sty = Com->sty ;
    pgnorm = Com->pgnorm ;
    n = Com->n ;
    lmin = Com->asaParm->lmin ;
    lmax = Com->asaParm->lmax ;

    if ( sty > ZERO )
    {
        t = MIN2 (sts/sty, lmax) ;
        alpha = MAX2 (lmin, t) ;
    }
    else
    {
        xnorm = asa_max (x, n) ;
        if ( xnorm > ZERO ) alpha = MIN2 (ONE, xnorm)/pgnorm ;
        else                alpha = ONE/pgnorm ;
    }
    return (alpha) ;
}

/* =========================================================================
   ==== asa_f ==============================================================
   Evaluate the function
   =========================================================================*/
double asa_f
(
    double    *x,
    asa_com *Com
)
{
    double f ;
    asa_objective *user ;
    user = Com->user ;
    user->x = x ;
    Com->nf++ ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate function */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        if(Com->objpoint!=NULL)
        f = Com->objpoint->myvalue(user) ;
        if(Com->avec_objpoint!=NULL)
        f = Com->avec_objpoint->myvalue(user) ;

        /* Shrink x to the reduced space */
        asa_shrinkx (x, Com) ;
    }
    else
    {
        /* Evaluate function */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        if(Com->objpoint!=NULL)
        f = Com->objpoint->myvalue(user) ;
        if(Com->avec_objpoint!=NULL)
        f = Com->avec_objpoint->myvalue(user) ;
    }
    return (f) ;

}

/* =========================================================================
   ==== asa_g ==============================================================
   Evaluate the gradient
   =========================================================================*/
 void asa_g
(
    double    *g,
    double    *x,
    asa_com *Com
)
{
    asa_objective *user ;
    user = Com->user ;
    user->x = x ;
    user->g = g ;
    Com->ng++ ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate gradient */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        if(Com->objpoint!=NULL)
        Com->objpoint->mygrad(user) ;
        if(Com->avec_objpoint!=NULL)
        Com->avec_objpoint->mygrad(user) ;

        /* Shrink x and g to the reduced space */
        asa_shrinkxg (x, g, Com) ;
    }
    else
    {
        /* Evaluate gradient */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        if(Com->objpoint!=NULL)
        Com->objpoint->mygrad(user) ;
        if(Com->avec_objpoint!=NULL)
        Com->avec_objpoint->mygrad(user) ;
    }
}


/* =========================================================================
   ==== asa_fg =============================================================
   Evaluate the function and gradient
   =========================================================================*/
double asa_fg
(
    double    *g,
    double    *x,
    asa_com *Com
)
{
    asa_objective *user ;
    double f ;
    Com->nf++ ;
    Com->ng++ ;
    user = Com->user ;
    user->x = x ;
    user->g = g ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate function and gradient */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        if ( Com->valgrad != NULL )
        {
            //f = Com->valgrad (user) ;
        }
        else
        {
            if(Com->objpoint!=NULL){
            Com->objpoint->mygrad(user) ;
            f = Com->objpoint->myvalue(user) ;
            }
            if(Com->avec_objpoint!=NULL){
            Com->avec_objpoint->mygrad(user) ;
            f = Com->avec_objpoint->myvalue(user) ;
            }
        }

        /* Shrink x and g to the reduced space */
        asa_shrinkxg (x, g, Com) ;
    }
    else
    {
        /* Evaluate function and gradient */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        if ( Com->valgrad != NULL )
        {
            //f = Com->valgrad (user) ;
        }
        else
        {
            if(Com->objpoint!=NULL){
            Com->objpoint->mygrad(user) ;
            f = Com->objpoint->myvalue(user) ;
            }
            if(Com->avec_objpoint!=NULL){
            Com->avec_objpoint->mygrad(user) ;
            f = Com->avec_objpoint->myvalue(user) ;
            }
        }
    }
    return (f) ;
}

/* =========================================================================
   ==== asa_identify =======================================================
   Check whether the bounds with strict complementarity
   are approximately identified
   =========================================================================*/
 int asa_identify
(
   double     *x,
   double     *g,
   double pgnorm,
   asa_com  *Com
)
{
    int ident ;
    ASA_INT j, n ;
    double t, t1, xj, *lo, *hi ;
    n = Com->n ;
    lo = Com->lo ;
    hi = Com->hi ;
    ident = TRUE ;
    t = sqrt (pgnorm) ;
    t1 = t*t*t ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        if ( ((xj - lo [j] >= t1) && (g [j] > t)) ||
             ((hi [j] - xj >= t1) && (g [j] <-t)) ) ident = FALSE ;
    }
    return (ident) ;
}

/* =========================================================================
   === asa_expandx =========================================================
   =========================================================================
   Expand x array from size nfree to full size of dimension n based on
   indices of free variables
   ========================================================================= */
 void asa_expandx
(
    double    *x,
    asa_com *Com
)
{
    ASA_INT i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = nfree-1; j >= 0; j--)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrinkx =========================================================
   =========================================================================
   Compress x array to dimension nfree based on indices of free variables
   ========================================================================= */
void asa_shrinkx
(
    double    *x,
    asa_com *Com
)
{
    ASA_INT i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrinkxg ========================================================
   =========================================================================
   Compress x and g arrays based on indices of free variables
   ========================================================================= */
 void asa_shrinkxg
(
    double    *x,
    double    *g,
    asa_com *Com
)
{
    ASA_INT i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_expand_all ======================================================
   =========================================================================
   Expand vectors x, g, pg, lo and hi from the reduced space (dimension nfree)
   to the full space (dimension n).
   ========================================================================= */
 void asa_expand_all
(
    asa_com *Com
)
{
    ASA_INT i, j, nfree, *ifree ;
    double t, *x, *g, *pg, *lo, *hi ;
    x = Com->x ;
    g = Com->g ;
    pg = Com->pg ;
    lo = Com->lo ;
    hi = Com->hi ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = nfree-1; j >= 0; j--)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;

            t = pg [i] ;
            pg [i] = pg [j] ;
            pg [j] = t ;

            t = lo [i] ;
            lo [i] = lo [j] ;
            lo [j] = t ;

            t = hi [i] ;
            hi [i] = hi [j] ;
            hi [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrink_all ======================================================
   =========================================================================
   Shrink vectors x, g, lo and hi from the full space (dimension n)
   to the reduced space (dimension nfree).
   ========================================================================= */

 void asa_shrink_all
(
    asa_com *Com
)
{
    ASA_INT i, j, nfree, *ifree ;
    double t, *lo, *hi, *g, *x ;
    x = Com->x ;
    g = Com->g ;
    lo = Com->lo ;
    hi = Com->hi ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( i != j )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;

            t = lo [i] ;
            lo [i] = lo [j] ;
            lo [j] = t ;

            t = hi [i] ;
            hi [i] = hi [j] ;
            hi [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_printcgParms ====================================================
   =========================================================================
   Print the contents of the asacg_parm structure
   ========================================================================= */
 void asa_printcgParms
(
    asacg_parm  *Parm
)
{
    printf ("\nCG PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("Wolfe line search parameter ..................... delta: %e\n",
             Parm->delta) ;
    printf ("Wolfe line search parameter ..................... sigma: %e\n",
             Parm->sigma) ;
    printf ("decay factor for bracketing interval ............ gamma: %e\n",
             Parm->gamma) ;
    printf ("growth factor for bracket interval ................ rho: %e\n",
             Parm->rho) ;
    printf ("growth factor for bracket interval after nan .. nan_rho: %e\n",
             Parm->nan_rho) ;
    printf ("decay factor for stepsize after nan ......... nan_decay: %e\n",
             Parm->nan_decay) ;
    printf ("parameter in lower bound for beta ........... BetaLower: %e\n",
             Parm->BetaLower) ;
    printf ("parameter describing cg_descent family .......... theta: %e\n",
             Parm->theta) ;
    printf ("perturbation parameter for function value ......... eps: %e\n",
             Parm->eps) ;
    printf ("factor by which eps grows if necessary .......... egrow: %e\n",
             Parm->egrow) ;
    printf ("factor for computing average cost .............. Qdecay: %e\n",
             Parm->Qdecay) ;
    printf ("relative change in cost to stop QuadStep ... QuadCutOff: %e\n",
             Parm->QuadCutOff) ;
    printf ("maximum factor quadstep reduces stepsize ..... QuadSafe: %e\n",
             Parm->QuadSafe) ;
    printf ("skip quadstep if |f| <= SmallCost*start cost  SmallCost: %e\n",
             Parm->SmallCost) ;
    printf ("relative change in cost to stop cubic step  CubicCutOff: %e\n",
             Parm->CubicCutOff) ;
    printf ("terminate if no improvement over nslow iter ..... nslow: %i\n",
             Parm->nslow) ;
    printf ("cost change factor, approx Wolfe transition . AWolfeFac: %e\n",
             Parm->AWolfeFac) ;
    printf ("restart cg every restart_fac*n iterations . restart_fac: %e\n",
             Parm->restart_fac) ;
    printf ("cost error in quadratic restart is qeps*cost ..... qeps: %e\n",
             Parm->qeps) ;
    printf ("number of quadratic iterations before restart  qrestart: %i\n",
             Parm->qrestart) ;
    printf ("parameter used to decide if cost is quadratic ... qrule: %e\n",
             Parm->qrule) ;
    printf ("stop when cost change <= feps*|f| ................. eps: %e\n",
             Parm->feps) ;
    printf ("starting guess parameter in first iteration ...... psi0: %e\n",
             Parm->psi0) ;
    printf ("lower bound factor in quad step ................ psi_lo: %e\n",
             Parm->psi_lo) ;
    printf ("upper bound factor in quad step ................ psi_hi: %e\n",
             Parm->psi_hi) ;
    printf ("factor multiply starting guess in quad step ...... psi1: %e\n",
             Parm->psi1) ;
    printf ("initial guess factor for general iteration ....... psi2: %e\n",
             Parm->psi2) ;
    printf ("starting step in first iteration if nonzero ...... step: %e\n",
             Parm->step) ;
    printf ("max tries to find non NAN function value ...... nshrink: %i\n",
             Parm->nshrink) ;
    printf ("max expansions in line search .................. ntries: %i\n",
             Parm->ntries) ;
    printf ("maximum growth of secant step in expansion . ExpandSafe: %e\n",
             Parm->ExpandSafe) ;
    printf ("growth factor for secant step during expand . SecantAmp: %e\n",
             Parm->SecantAmp) ;
    printf ("growth factor for rho during expansion phase .. RhoGrow: %e\n",
             Parm->RhoGrow) ;
    printf ("distance threshhold for entering subspace ........ eta0: %e\n",
             Parm->eta0) ;
    printf ("distance threshhold for leaving subspace ......... eta1: %e\n",
             Parm->eta1) ;
    printf ("distance threshhold for invariant space .......... eta2: %e\n",
             Parm->eta2) ;
    printf ("number of vectors stored in memory ............. memory: %i\n",
             Parm->memory) ;
    printf ("check subspace condition mem*SubCheck its .... SubCheck: %i\n",
             Parm->SubCheck) ;
    printf ("skip subspace checking for mem*SubSkip its .... SubSkip: %i\n",
             Parm->SubSkip) ;
    printf ("max number of times that eps is updated .......... neps: %i\n",
             Parm->neps) ;
    printf ("max number of iterations in line search ......... nline: %i\n",
             Parm->nline) ;
    printf ("error tolerance when debugger turned on ..... .debugtol: %e\n",
             Parm->debugtol) ;
    printf ("print level (0 = none, 4 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("\nLogical parameters:\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps*Ck\n") ;
    else
        printf ("    Error estimate for function value is eps\n") ;
    if ( Parm->QuadStep )
        printf ("    Use quadratic interpolation step\n") ;
    else
        printf ("    No quadratic interpolation step\n") ;
    if ( Parm->UseCubic)
        printf ("    Use cubic interpolation step when possible\n") ;
    else
        printf ("    Avoid cubic interpolation steps\n") ;
    if ( Parm->AdaptiveBeta )
        printf ("    Adaptively adjust direction update parameter beta\n") ;
    else
        printf ("    Use fixed parameter theta in direction update\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AWolfe)
        printf ("    Approximate Wolfe line search\n") ;
    else
        printf ("    Wolfe line search") ;
        if ( Parm->AWolfeFac > ZERO )
            printf (" ... switching to approximate Wolfe\n") ;
        else
            printf ("\n") ;
    if ( Parm->debug)
        printf ("    Check for decay of cost, debugger is on\n") ;
    else
        printf ("    Do not check for decay of cost, debugger is off\n") ;
}

/* =========================================================================
   === asa_default =========================================================
   =========================================================================
   Set default parameter values for the ASA routine. The CG default
   parameter values are set by asa_cg_default.  If the parameter argument of
   asa_descent is NULL, this routine is called by asa_cg automatically.
   If the user wishes to set parameter values, then the asa_parameter structure
   should be allocated in the main program. The user could call asa_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to asa_cg.
   =========================================================================*/
void asa_default
(
    asa_parm *Parm
)
{
    double eps, t ;

    /* T => print final statistics
       F => no printout of statistics */
    Parm->PrintFinal = TRUE ;

    /* Level 0  = no printing), ... , Level 4 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parameter values */
    Parm->PrintParms = FALSE ;

    /* T => use approximate nonmonotone Armijo line search
       F => use ordinary nonmonotone Armijo line search, switch to
            approximate Armijo when |f_r-f| < AArmijoFac*|min (f_r, f_{max})| */
    Parm->AArmijo = FALSE ;
    Parm->AArmijoFac = 1.e-8 ;

    /* Stop Rules (these override the corresponding cg parameters):
       T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFac)
       F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
    Parm->StopRule = TRUE ;
    Parm->StopFac = 0.e-12 ;

    /* T => estimated error in function value = eps*|min (f_r, f_{max}) |
       F => estimated error in function value = eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* T => only use gradient projection algorithm
       F => let algorithm decide between grad_proj and cg_descent */
    Parm->GradProjOnly = FALSE ;

    /* maximum number of times the Armijo line search will perform
       backtracking steps */
    Parm->max_backsteps = (int) 50 ;

    /* abort cbb after maxit_fac*n iterations in one pass through cbb */
    Parm->maxit_fac = ASA_INF ;

    /* abort cbb after totit_fac*n iterations in all passes through cbb */
    Parm->totit_fac = ASA_INF ;

    /* abort cbb iteration after maxfunc_fac*n function evaluations */
    Parm->maxfunc_fac = ASA_INF ;

    /* perturbation in bounds based on machine epsilon, which we now compute */
    eps = ONE ;
    t = ONE ;
    while ( t > 0 )
    {
        eps /= TWO ;
        t = ONE + eps ;
        t -= ONE ;
    }
    eps *= 2 ;                   /* machine epsilon */
    Parm->pert_lo = 1.e3*eps ;   /* perturbation of lower bounds */
    Parm->pert_hi = 1.e3*eps ;   /* perturbation of upper bounds */

    /* search for non nan function value by shrinking search interval
       at most nshrink times */
    Parm->nshrink = (int) 50 ;

    /* factor by which interval shrinks when searching for non nan value */
    Parm->nan_fac = 2.e-1 ;

    /* T => Dynamically assign memory in cg when solving subproblem
       F => Enough memory was assigned before solving the problem */
    Parm->DynamicMemory = FALSE ;

    /* In cg_descent, when calculating the initial stepsize,
       we sample the function or gradient at another point
       along the search direction.  This sampling point could
       violate the box constraints.  If HardConstraint is FALSE,
       then we evaluate the function and gradient at this sampling
       point even when it violates the box constraints. If HardConstraint
       is TRUE, then the sampling point is chosen to satisfy the the
       box constraints. In versions of the code that predate version 3.0,
       the sampling point was always chosen to satisfy the box constraints
       (HardConstraint = TRUE).

       T => Hard bound constraints
       F => Soft bound constraints */
    Parm->HardConstraint = TRUE ;

    /* update fr if fmin was not improved after L iterations */
    Parm->L = 3 ;

    /* fmax = max (f_{k-i}, i = 0, 1, ..., min (k, m-1) ) */
    Parm->m = 8 ;

    /* update fr if initial stepsize was accepted in previous P iterations */
    Parm->P = 40 ;

    /* CBB cycle length */
    Parm->nm = 4 ;

    /* Reinitialize BB stepsize, if (s^t y)/(||s|| ||y||) >= gamma
       and ||s|| <= min (parm3*|f_k+1|/||g_k+1||_infty, 1) */
    Parm->gamma = 0.975e0 ;

    /* update reference value fr if (fr-fmin)/(fc-fmin) > gamma1 */
    Parm->gamma1 = (double) Parm->m / (double) Parm->L ;

    /* update fr if (fr-f)/(fmax-f) > gamma2, np > P, and fmax > f */
    Parm->gamma2 = (double) Parm->P / (double) Parm->m ;

    /* terminate Armijo line search when
       phi(alpha) <= phi_r + alpha * delta * phi'(0) where phi_r = fr or fcomp*/
    Parm->delta = 1.0e-4 ;   /* Armijo line search parameter */

    /* stepsize s in the line search must satisfy lmin <= s <= lmax */
    Parm->lmin = 1.0e-20 ;
    Parm->lmax = 1.0e+20 ;

    /* attempt a quadratic interpolation step in cg_descent if the
       provisional stepsize times parm1 <= stepsize to boundary */
    Parm->parm1 = 1.e-1 ;

    /* if quadratic interpolation step is attempted, the provisional step
       is at most parm2*stepsize to boundary */
    Parm->parm2 = 9.e-1 ;

    /* used in the the criterion of reinitializing the BB stepsize */
    Parm->parm3 = 1.e-1 ;

    /* maximum number of previous BB steps used when s^t y <= ZERO */
    Parm->parm4 = 6 ;

    /* if ginorm < tau1*pgnorm, continue gradient projection steps  */
    Parm->tau1 = 1.e-1 ;

    /* decay factor for tau1 */
    Parm->tau1_decay = 5.e-1 ;

    /* ginorm < tau2*pgnorm implies subproblem solved in cgdescent */
    Parm->tau2 = 1.e-1 ;

    /* decay factor for tau2 */
    Parm->tau2_decay = 5.e-1 ;

    /* if pgnorm < pgdecay*MAX (pgnorm0, ONE), check the undecided index set
                                pgnorm0 = pgnorm at starting point */
    Parm->pgdecay = 1.e-4 ;

    /* backtracking decay factor in the Armijo line search */
    Parm->armijo_decay = 5.e-1 ;

    /* use quadratic interpolation to compute Armijo step if it
       lies in the interval [.1 alpha, .9 alpha] */
    Parm->armijo0 = 1.e-1 ;
    Parm->armijo1 = 9.e-1 ;
}

/* =========================================================================
   === asa_cg_default ======================================================
   =========================================================================
   Set default conjugate gradient parameter values. If the parameter argument
   of asa_cg is NULL, this routine is called by asa_cg automatically.
   If the user wishes to set parameter values, then the asa_parameter structure
   should be allocated in the main program. The user could call asa_cg_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to asa_cg.
   =========================================================================*/

 void asa_cg_default
(
    asacg_parm   *Parm
)
{
    /* Level 0 = no printing, ... , Level 4 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parameter values */
    Parm->PrintParms = FALSE ;

    /* T => use LBFGS
       F => only use L-BFGS when memory >= n */
    Parm->LBFGS = FALSE ;

    /* number of vectors stored in memory (code breaks in the Yk update if
       memory = 1 or 2) */
    Parm->memory = 11 ;

    /* SubCheck and SubSkip control the frequency with which the subspace
       condition is checked. It it checked for SubCheck*mem iterations and
       if it is not activated, then it is skipped for Subskip*mem iterations
       and Subskip is doubled. Whenever the subspace condition is satisfied,
       SubSkip is returned to its original value. */
    Parm->SubCheck = 8 ;
    Parm->SubSkip = 4 ;

    /* when relative distance from current gradient to subspace <= eta0,
       enter subspace if subspace dimension = mem (eta0 = 0 means gradient
       inside subspace) */
    Parm->eta0 = 0.001 ;

    /* when relative distance from current gradient to subspace >= eta1,
       leave subspace (eta1 = 1 means gradient orthogonal to subspace) */
    Parm->eta1 = 0.900 ;

    /* when relative distance from current direction to subspace <= eta2,
       always enter subspace (invariant space) */
    Parm->eta2 = 1.e-10 ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost */
    Parm->AWolfe = FALSE ;
    Parm->AWolfeFac = 1.e-3 ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    Parm->Qdecay = 0.7 ;

    /* terminate after 2*n + nslow iterations without strict improvement in
       either function value or gradient */
    Parm->nslow = 1000 ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* factor by which eps grows when line search fails during contraction */
    Parm->egrow = 10. ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutOff
       F => no quadratic interpolation step */
    Parm->QuadStep = TRUE ;
    Parm->QuadCutOff = 1.e-12 ;

    /* maximum factor by which a quad step can reduce the step size */
    Parm->QuadSafe = 1.e-10 ;

    /* T => when possible, use a cubic step in the line search */
    Parm->UseCubic = TRUE ;

    /* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
    Parm->CubicCutOff = 1.e-12 ;

    /* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
    Parm->SmallCost = 1.e-30 ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    Parm->debug = FALSE ;
    Parm->debugtol = 1.e-10 ;

    /* if step is nonzero, it is the initial step of the initial line search */
    Parm->step = ZERO ;

    /* abort cg after maxit iterations */
    Parm->maxit = ASA_INT_INF ;

    /* maximum number of times the bracketing interval grows during expansion */
    Parm->ntries = (int) 50 ;

    /* maximum factor secant step increases stepsize in expansion phase */
    Parm->ExpandSafe = 200. ;

    /* factor by which secant step is amplified during expansion phase
       where minimizer is bracketed */
    Parm->SecantAmp = 1.05 ;

    /* factor by which rho grows during expansion phase where minimizer is
       bracketed */
    Parm->RhoGrow = 2.0 ;

    /* maximum number of times that eps is updated */
    Parm->neps = (int) 5 ;

    /* maximum number of times the bracketing interval shrinks */
    Parm->nshrink = (int) 10 ;

    /* maximum number of secant iterations in line search is nline */
    Parm->nline = (int) 50 ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    Parm->restart_fac = 6.0 ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    Parm->feps = ZERO ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    Parm->nan_rho = 1.3 ;

    /* after encountering nan, decay factor for stepsize */
    Parm->nan_decay = 0.1 ;

    /* Wolfe line search parameter, range [0, .5]
       phi (a) - phi (0) <= delta phi'(0) */
    Parm->delta = .1 ;

    /* Wolfe line search parameter, range [delta, 1]
       phi' (a) >= sigma phi' (0) */
    Parm->sigma = .9 ;

    /* decay factor for bracket interval width in line search, range (0, 1) */
    Parm->gamma = .66 ;

    /* growth factor in search for initial bracket interval */
    Parm->rho = 5. ;

    /* starting guess for line search =
         psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
         psi0 |f(x_0)|/||g_0||_2               otherwise */
    Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

    /* for a QuadStep, function evaluated on interval
       [psi_lo, phi_hi]*psi2*previous step */
    Parm->psi_lo = 0.1 ;
    Parm->psi_hi = 10. ;

    /* when the function is approximately quadratic, use gradient at
       psi1*psi2*previous step for estimating initial stepsize */
    Parm->psi1 = 1.0 ;

    /* when starting a new cg iteration, our initial guess for the line
       search stepsize is psi2*previous step */
    Parm->psi2 = 2. ;

    /* choose theta adaptively if AdaptiveBeta = T */
    Parm->AdaptiveBeta = FALSE ;

    /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
    Parm->BetaLower = 0.4 ;

    /* value of the parameter theta in the cg_descent update formula:
       W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
       methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */
    Parm->theta = 1.0 ;

    /* parameter used in cost error estimate for quadratic restart criterion */
    Parm->qeps = 1.e-12 ;

    /* number of iterations the function is nearly quadratic before a restart */
    Parm->qrestart = 6 ;

    /* treat cost as quadratic if
       |1 - (cost change)/(quadratic cost change)| <= qrule */
    Parm->qrule = 1.e-8 ;
}

/* =========================================================================
   === asa_printParms ====================================================
   =========================================================================
   Print the contents of the asa_parm structure
   ========================================================================= */
 void asa_printParms
(
    asa_parm  *Parm
)
{
    printf ("\nASA PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("update fr if fmin not improved after L iterations.... L: %i\n",
             Parm->L) ;
    printf ("fmax = max (f_{k-i}, i = 0, 1, ..., min (k, m-1) )... m: %i\n",
             Parm->m) ;
    printf ("update fr if P previous initial stepsizes accepted... P: %i\n",
             Parm->P) ;
    printf ("CBB cycle length.................................... nm: %i\n",
             Parm->nm) ;
    printf ("criterion for updating reference value fr....... gamma1: %e\n",
             Parm->gamma1) ;
    printf ("criterion for updating reference value fr....... gamma2: %e\n",
             Parm->gamma2) ;
    printf ("interval decay factor in NAN search ............nan_fac: %e\n",
             Parm->nan_fac) ;
    printf ("perturbation parameter for function value.......... eps: %e\n",
             Parm->eps) ;
    printf ("cost change factor, approx Armijo transition,AArmijoFac: %e\n",
             Parm->AArmijoFac) ;
    printf ("Armijo line search parameter .................... delta: %e\n",
             Parm->delta) ;
    printf ("Armijo decay factor .......................armijo_decay: %e\n",
             Parm->armijo_decay) ;
    printf ("criterion for Q interpolation, cbb line search,.armijo0: %e\n",
             Parm->armijo0) ;
    printf ("criterion for Q interpolation, cbb line search,.armijo1: %e\n",
             Parm->armijo1) ;
    printf ("criterion for reinitializing BB stepsize ........ gamma: %e\n",
             Parm->gamma) ;
    printf ("Lower bound for initial stepsize ................. lmin: %e\n",
             Parm->lmin) ;
    printf ("Upper bound for initial stepsize ................. lmax: %e\n",
             Parm->lmax) ;
    printf ("used when trying a quadratic interpolation step.. parm1: %e\n",
             Parm->parm1) ;
    printf ("used when trying a quadratic interpolation step.. parm2: %e\n",
             Parm->parm2) ;
    printf ("criterion for reinitializing the BB stepsize..... parm3: %e\n",
             Parm->parm3) ;
    printf ("maximum previous BB steps used when s^t y <= 0... parm4: %i\n",
             Parm->parm4) ;
    printf ("if ginorm < tau1*pgnorm, continue grad_proj ...... tau1: %e\n",
             Parm->tau1) ;
    printf ("decay factor for tau1 ...................... tau1_decay: %e\n",
             Parm->tau1_decay) ;
    printf ("ginorm < tau2*pgnorm => subproblem solved in cg... tau2: %e\n",
             Parm->tau2) ;
    printf ("decay factor for tau2 ...................... tau2_decay: %e\n",
             Parm->tau2_decay) ;
    printf ("max number of Armijo backtracking steps . max_backsteps: %i\n",
             Parm->max_backsteps) ;
    printf ("max cbb iterations in 1 pass is n*maxit_fac . maxit_fac: %e\n",
             Parm->maxit_fac) ;
    printf ("max number of contracts in the line search .... nshrink: %i\n",
             Parm->nshrink) ;
    printf ("total number cbb iterations is n*totit_fac .. totit_fac: %e\n",
             Parm->totit_fac) ;
    printf ("max func evals in cbb is n*maxfunc_fac .... maxfunc_fac: %e\n",
             Parm->maxfunc_fac) ;
    printf ("criterion for checking undecided index set..... pgdecay: %e\n",
             Parm->pgdecay) ;
    printf ("perturbation of lower bounds .................. pert_lo: %e\n",
             Parm->pert_lo) ;
    printf ("perturbation of upper bounds .................. pert_hi: %e\n",
             Parm->pert_hi) ;
    printf ("factor multiplying gradient in stop condition . StopFac: %e\n",
             Parm->StopFac) ;
    printf ("print level (0 = none, 4 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("\nLogical parameters:\n") ;
    if ( Parm->DynamicMemory)
        printf ("    dynamically assign cg memory when solving subproblem\n") ;
    else
        printf ("    assign memory before solving the problem\n") ;
    if ( Parm->HardConstraint )
        printf ("    Bound constraints are hard constraints\n") ;
    else
        printf ("    Bound constraints are not hard constraints\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps*|fcomp|\n") ;
    else
        printf ("    Error estimate for function value is eps\n") ;
    if ( Parm->PrintFinal )
        printf ("    Print final cost and statistics\n") ;
    else
        printf ("    Do not print final cost and statistics\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AArmijo)
        printf ("    Approximate nonmonotone Armijo line search\n") ;
    else
        printf ("    Nonmonotone Armijo line search") ;
        if ( Parm->AArmijoFac > ZERO )
            printf (" ... switching to approx nonmonotone Armijo\n") ;
        else
            printf ("\n") ;
    if ( Parm->StopRule )
    {
        if ( Parm->StopFac == ZERO )
        {
            printf ("    Stopping condition based on gradient tolerance\n") ;
        }
        else
        {
            printf ("    Stopping condition uses initial grad tolerance\n") ;
        }
    }
    else
        printf ("    Stopping condition weighted by absolute cost\n") ;
    if ( Parm->GradProjOnly )
        printf ("    Only use the gradient projection algorithm\n") ;
    else
        printf ("    Apply gradient projection algorithm and cg_descent\n") ;
}
/*
Version 1.1 Change:
    1. Pass a structure asa_objective to the user evaluation routines.
       This allows asa_cg to pass more information to the user which
       might be used to speedup his routines to evaluate the objective
       function and its gradient.  Two elements of the structure are
       ifree and nfree.  If ifree is not NULL, then ifree is a pointer
       to an integer array containing the indices of the free variables
       while nfree is the number of free variables.

    2. Halt the Armijo backtracking line search in cbb when the number
       of backtracking steps reaching Parm->max_backsteps

Version 1.2:
    Correct the Armijo line search in cbb by dividing by including the
    factor "/alpha" in the termination condition

Version 1.3:
    In asa_identify, correct the formula for identifying constraints.

Version 2.0:
    Update the cg routine utilizing cg_descent 4.0

Version 2.1:
    Correct several bugs connected with installation of cg_descent 4.0

Version 2.2:
    1. Modify cg_descent line search so that when there are too many
       contractions, the code will increase eps and switch to expansion
       of the search interval.  This fixes some cases where the code
       terminates when eps is too small.  When the estimated error in
       the cost function is too small, the algorithm could fail in
       cases where the slope is negative at both ends of the search
       interval and the objective function value on the right side of the
       interval is larger than the value at the left side (because the true
       objective function value on the right is not greater than the value on
       the left).
    2. Fix bug in asa_lineW

Version 3.0:
    1. Update the conjugate gradient routine using cg_descent version 6.0.
    2. Add the parameter HardConstraint. If HardConstraint is FALSE, then
       cg_descent is allowed to evaluate the function or gradient at a
       point that violates the box constraints. This could improve
       performance by giving the code additional flexibility in the
       initial stepsize routine.
*/
