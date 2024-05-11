
#include "pc-iluk.h"
#include "itsol.h"


/*----------------------------------------------------------------------------
 * ILUK preconditioner
 * incomplete LU factorization with level of fill dropping
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill: all entries with level of fill > lofM are
 *            dropped. Setting lofM = 0 gives BILU(0).
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 *            on format
 * lu       = pointer to a ILUKSpar struct -- see heads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> error in lofC
 *            ierr  = -2  --> zero diagonal found
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *--------------------------------------------------------------------------*/


ITS_ILUSpar *itsol_from_numpy(int n, int lfil, double *data, int *indices, int *indptr)
{
    ITS_SparMat *csmat = NULL;  /* matrix in csr formt             */
    ITS_ILUSpar *lu = NULL;     /* ilu preconditioner structure    */
    ITS_SMat *MAT;              /* Matrix structure for matvecs    */
    ITS_PC *PRE;                /* general precond structure       */
    ITS_PARS io;
    double *x, *rhs, *sol, terr, norm;
    int ierr, i, k, its;
    csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main");
    lu = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "main");
    printf("begin iluk(%d)\n", lfil);
    /*-------------------- conversion from COO to CSR format */
    csmat->n = n;
    csmat->nzcount = (int *) itsol_malloc(sizeof(int)*n, "main");
    csmat->ja = (int **) itsol_malloc(sizeof(int*)*n, "main");
    csmat->ma = (double **) itsol_malloc(sizeof(double*)*n, "main");
    for(i=0; i < n; i++)
    {
        csmat->nzcount[i] = indptr[i+1]-indptr[i];
        csmat->ja[i] = (int *) itsol_malloc(sizeof(int)*csmat->nzcount[i], "main");
        csmat->ma[i] = (double *) itsol_malloc(sizeof(double)*csmat->nzcount[i], "main");
        for(k=0; k < csmat->nzcount[i]; k++)
        {
            (csmat->ja)[i][k] = indices[indptr[i]+k];
            (csmat->ma)[i][k] = data[indptr[i]+k];
        }
         
    }
    ierr = itsol_pc_ilukC(lfil, csmat, lu, stdout);
    


    //printf("%d\n", lu->L->nzcount[5]);
    /*MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");
    itsol_solver_init_pars(&io);
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");
    for( i = 0; i < n; i++ ) x[i] = 0.0;
    for (i = 0; i < n; i++) rhs[i] = i;
    MAT->n = n;
    MAT->CS = csmat;
    MAT->matvec = itsol_matvecCSR;
    PRE->ILU = lu;
    PRE->precon = itsol_preconILU;
    */
    /*-------------------- call itsol_solver_fgmres */
    
    /*
    itsol_solver_fgmres(MAT, PRE, rhs, x, io, &its, NULL);

    printf("solver converged in %d steps...\n\n", its);
    */
    /*-------------------- calculate residual norm */
    //itsol_matvec(csmat, x, sol);

    /* error */
    /*
    terr = 0.0;
    norm = 0.;
    for (i = 0; i < n; i++) {
        terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

        norm += rhs[i] * rhs[i];
    }

    printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));
    */
    return lu;
    //return ierr;
    //printf("%d\n:", n);
    //for(int i=0;i<=n;i++)
    //{
    //    printf("%d, %d, %f\n:", i, indices[i], data[i]);
    //}
    

}


int itsol_pc_ilukC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE * fp)
{
    int ierr;
    int n = csmat->n;
    int *jw, i, j, k, col, jpos, jrow;
    ITS_SparMat *L, *U;
    double *D;
    itsol_setupILU(lu, n);
    //printf("%d", n);//lu->L->nzcount[5]);
    L = lu->L;
    U = lu->U;
    D = lu->D;
    
    /* symbolic factorization to calculate level of fill index arrays */
    if ((ierr = itsol_pc_lofC(lofM, csmat, lu, fp)) != 0) {
        fprintf(fp, "Error: lofC\n");
        return -1;
    }
    //printf("%d\n", ierr);
    //printf("%d\n", lu->L->nzcount[0]);

    jw = lu->work;
    /* set indicator array jw to -1 */
    for (j = 0; j < n; j++)
        jw[j] = -1;

    /* beginning of main loop */
    for (i = 0; i < n; i++) {
        /* set up the i-th row accroding to the nonzero information from
           symbolic factorization */
        itsol_mallocRow(lu, i);

        /* setup array jw[], and initial i-th row */
        for (j = 0; j < L->nzcount[i]; j++) {   /* initialize L part   */
            col = L->ja[i][j];
            jw[col] = j;
            L->ma[i][j] = 0;
        }
        jw[i] = i;
        D[i] = 0;               /* initialize diagonal */
        for (j = 0; j < U->nzcount[i]; j++) {   /* initialize U part   */
            col = U->ja[i][j];
            jw[col] = j;
            U->ma[i][j] = 0;
        }

        /* copy row from csmat into lu */
        for (j = 0; j < csmat->nzcount[i]; j++) {
            col = csmat->ja[i][j];
            jpos = jw[col];
            if (col < i)
                L->ma[i][jpos] = csmat->ma[i][j];
            else if (col == i)
                D[i] = csmat->ma[i][j];
            else
                U->ma[i][jpos] = csmat->ma[i][j];
        }

        /* eliminate previous rows */
        for (j = 0; j < L->nzcount[i]; j++) {
            jrow = L->ja[i][j];
            /* get the multiplier for row to be eliminated (jrow) */
            L->ma[i][j] *= D[jrow];

            /* combine current row and row jrow */
            for (k = 0; k < U->nzcount[jrow]; k++) {
                col = U->ja[jrow][k];
                jpos = jw[col];
                if (jpos == -1)
                    continue;
                if (col < i)
                    L->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
                else if (col == i)
                    D[i] -= L->ma[i][j] * U->ma[jrow][k];
                else
                    U->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
            }
        }

        /* reset double-pointer to -1 ( U-part) */
        for (j = 0; j < L->nzcount[i]; j++) {
            col = L->ja[i][j];
            jw[col] = -1;
        }
        jw[i] = -1;
        for (j = 0; j < U->nzcount[i]; j++) {
            col = U->ja[i][j];
            jw[col] = -1;
        }

        if (D[i] == 0) {
            for (j = i + 1; j < n; j++) {
                L->ma[j] = NULL;
                U->ma[j] = NULL;
            }
            fprintf(fp, "fatal error: Zero diagonal found...\n");
            return -2;
        }
        D[i] = 1.0 / D[i];
    }

    return 0;
}

/*--------------------------------------------------------------------
 * symbolic ilu factorization to calculate structure of ilu matrix
 * for specified level of fill
 *--------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill, lofM >= 0
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 *            on format
 * lu       = pointer to a ILUSpar struct -- see heads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *--------------------------------------------------------------------
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr != 0   --> error
 * lu->n    = dimension of the block matrix
 *   ->L    = L part -- stored in SpaFmt format, patterns only in lofC
 *   ->U    = U part -- stored in SpaFmt format, patterns only in lofC
 *------------------------------------------------------------------*/
int itsol_pc_lofC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE * fp)
{
    int n = csmat->n;
    int *levls = NULL, *jbuf = NULL, *iw = lu->work;
    int **ulvl;                 /*  stores lev-fils for U part of ILU factorization */
    ITS_SparMat *L = lu->L, *U = lu->U;
    /*--------------------------------------------------------------------
     * n        = number of rows or columns in matrix
     * inc      = integer, count of nonzero(fillin) element of each row
     *            after symbolic factorization
     * ju       = entry of U part of each row
     * lvl      = buffer to store levels of each row
     * jbuf     = buffer to store column index of each row
     * iw       = work array
     *------------------------------------------------------------------*/
    int i, j, k, col, ip, it, jpiv;
    int incl, incu, jmin, kmin;

    (void)fp;
    levls = (int *)itsol_malloc(n * sizeof(int), "lofC");
    jbuf = (int *)itsol_malloc(n * sizeof(int), "lofC");
    ulvl = (int **)itsol_malloc(n * sizeof(int *), "lofC");

    /* initilize iw */
    for (j = 0; j < n; j++)
        iw[j] = -1;
    for (i = 0; i < n; i++) {
        incl = 0;
        incu = i;
        /*-------------------- assign lof = 0 for matrix elements */
        for (j = 0; j < csmat->nzcount[i]; j++) {
            col = csmat->ja[i][j];
            if (col < i) {
                /*-------------------- L-part  */
                jbuf[incl] = col;
                levls[incl] = 0;
                iw[col] = incl++;
            }
            else if (col > i) {
                /*-------------------- U-part  */
                jbuf[incu] = col;
                levls[incu] = 0;
                iw[col] = incu++;
            }
        }
        /*-------------------- symbolic k,i,j Gaussian elimination  */
        jpiv = -1;
        while (++jpiv < incl) {
            k = jbuf[jpiv];
            /*-------------------- select leftmost pivot */
            kmin = k;
            jmin = jpiv;
            for (j = jpiv + 1; j < incl; j++) {
                if (jbuf[j] < kmin) {
                    kmin = jbuf[j];
                    jmin = j;
                }
            }
            /*-------------------- swap  */
            if (jmin != jpiv) {
                jbuf[jpiv] = kmin;
                jbuf[jmin] = k;
                iw[kmin] = jpiv;
                iw[k] = jmin;
                j = levls[jpiv];
                levls[jpiv] = levls[jmin];
                levls[jmin] = j;
                k = kmin;
            }
            /*-------------------- symbolic linear combinaiton of rows  */
            for (j = 0; j < U->nzcount[k]; j++) {
                col = U->ja[k][j];
                it = ulvl[k][j] + levls[jpiv] + 1;
                if (it > lofM)
                    continue;
                ip = iw[col];
                if (ip == -1) {
                    if (col < i) {
                        jbuf[incl] = col;
                        levls[incl] = it;
                        iw[col] = incl++;
                    }
                    else if (col > i) {
                        jbuf[incu] = col;
                        levls[incu] = it;
                        iw[col] = incu++;
                    }
                }
                else
                    levls[ip] = its_min(levls[ip], it);
            }
        }                       /* end - while loop */
        /*-------------------- reset iw */
        for (j = 0; j < incl; j++)
            iw[jbuf[j]] = -1;
        for (j = i; j < incu; j++)
            iw[jbuf[j]] = -1;
        /*-------------------- copy L-part */
        L->nzcount[i] = incl;
        if (incl > 0) {
            L->ja[i] = (int *)itsol_malloc(incl * sizeof(int), "lofC");
            memcpy(L->ja[i], jbuf, sizeof(int) * incl);
        }
        /*-------------------- copy U - part        */
        k = incu - i;
        U->nzcount[i] = k;
        if (k > 0) {
            U->ja[i] = (int *)itsol_malloc(sizeof(int) * k, "lofC");
            memcpy(U->ja[i], jbuf + i, sizeof(int) * k);
            /*-------------------- update matrix of levels */
            ulvl[i] = (int *)itsol_malloc(k * sizeof(int), "lofC");
            memcpy(ulvl[i], levls + i, k * sizeof(int));
        }
    }

    /*-------------------- free temp space and leave --*/
    free(levls);
    free(jbuf);
    for (i = 0; i < n - 1; i++) {
        if (U->nzcount[i])
            free(ulvl[i]);
    }
    free(ulvl);

    return 0;
}
