
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


typedef long ptrdiff_t_blas;   
#define BINT ptrdiff_t_blas

extern void   dgemm_(const char*, const char*, const BINT*, const BINT*, const BINT*,
                     const double*, const double*, const BINT*, const double*, const BINT*,
                     const double*, double*, const BINT*);
extern void   dgemv_(const char*, const BINT*, const BINT*, const double*, const double*,
                     const BINT*, const double*, const BINT*, const double*, double*, const BINT*);
extern void   dcopy_(const BINT*, const double*, const BINT*, double*, const BINT*);
extern void   daxpy_(const BINT*, const double*, const double*, const BINT*, double*, const BINT*);
extern void   dscal_(const BINT*, const double*, double*, const BINT*);
extern double ddot_ (const BINT*, const double*, const BINT*, const double*, const BINT*);
extern void   dpotrf_(const char*, const BINT*, double*, const BINT*, BINT*);
extern void   dpotrs_(const char*, const BINT*, const BINT*, const double*, const BINT*,
                      double*, const BINT*, BINT*);

static const BINT   IONE  = 1;
static const double DONE  = 1.0;
static const double DMONE = -1.0;
static const double DZERO = 0.0;

typedef struct {
    BINT   iterations;
    double mu;
    double duality;
    int    status;      /* 0 converged, 1 max-iter, -1 Cholesky failure */
    BINT   chol_info;
} BoxQPInfo;

typedef struct {
    BINT n;
    double *gamma, *theta, *phi, *psi, *grad;
    double *Hbar, *rhs, *r1, *r2;
    double *dz_a, *dg_a, *dt_a, *dp_a, *dq_a;
    double *dz, *dg, *dt, *dp, *dq;
    double *z_shift;
} BoxQPWorkspace;

static BoxQPWorkspace g_boxqp_workspace = {0};

static void boxqp_workspace_free(BoxQPWorkspace *w)
{
    if (!w) return;
    free(w->gamma); free(w->theta); free(w->phi); free(w->psi); free(w->grad);
    free(w->Hbar); free(w->rhs); free(w->r1); free(w->r2);
    free(w->dz_a); free(w->dg_a); free(w->dt_a); free(w->dp_a); free(w->dq_a);
    free(w->dz); free(w->dg); free(w->dt); free(w->dp); free(w->dq);
    free(w->z_shift);
    memset(w, 0, sizeof(*w));
}

static int boxqp_workspace_prepare(BoxQPWorkspace *w, BINT n)
{
    if (!w || n <= 0) return -2;
    if (w->n == n && w->gamma && w->Hbar && w->z_shift) return 0;

    boxqp_workspace_free(w);
    w->n = n;
    w->gamma = calloc((size_t)n, sizeof(double));
    w->theta = calloc((size_t)n, sizeof(double));
    w->phi   = calloc((size_t)n, sizeof(double));
    w->psi   = calloc((size_t)n, sizeof(double));
    w->grad  = calloc((size_t)n, sizeof(double));
    w->Hbar  = calloc((size_t)n*n, sizeof(double));
    w->rhs   = calloc((size_t)n, sizeof(double));
    w->r1    = calloc((size_t)n, sizeof(double));
    w->r2    = calloc((size_t)n, sizeof(double));
    w->dz_a  = calloc((size_t)n, sizeof(double));
    w->dg_a  = calloc((size_t)n, sizeof(double));
    w->dt_a  = calloc((size_t)n, sizeof(double));
    w->dp_a  = calloc((size_t)n, sizeof(double));
    w->dq_a  = calloc((size_t)n, sizeof(double));
    w->dz    = calloc((size_t)n, sizeof(double));
    w->dg    = calloc((size_t)n, sizeof(double));
    w->dt    = calloc((size_t)n, sizeof(double));
    w->dp    = calloc((size_t)n, sizeof(double));
    w->dq    = calloc((size_t)n, sizeof(double));
    w->z_shift = calloc((size_t)n, sizeof(double));

    if (!w->gamma||!w->theta||!w->phi||!w->psi||!w->grad||!w->Hbar||!w->rhs||
        !w->r1||!w->r2||!w->dz_a||!w->dg_a||!w->dt_a||!w->dp_a||!w->dq_a||
        !w->dz||!w->dg||!w->dt||!w->dp||!w->dq||!w->z_shift) {
        boxqp_workspace_free(w);
        return -2;
    }
    return 0;
}

void deepc_boxqp_release_workspace(void)
{
    boxqp_workspace_free(&g_boxqp_workspace);
}

static double clip_box_interior(double v)
{
    const double margin = 1e-8;
    if (v >  1.0 - margin) return  1.0 - margin;
    if (v < -1.0 + margin) return -1.0 + margin;
    return v;
}

/* grad = H*z + h   (H dense n x n, column-major) */
static void compute_grad(const double *H, const double *h, const double *z,
                         BINT n, double *grad)
{
    dcopy_(&n, h, &IONE, grad, &IONE);
    dgemv_("n", &n, &n, &DONE, H, &n, z, &IONE, &DONE, grad, &IONE);
}

static double max_abs_vec(const double *x, BINT n)
{
    BINT i; double out = 0.0;
    for (i = 0; i < n; ++i) { double a = fabs(x[i]); if (a > out) out = a; }
    return out;
}

/* Fraction-to-boundary step length: largest alpha in (0,1] keeping all dual/slack variables strictly positive (0.99 safety factor). */
static double fraction_to_boundary(const double *gamma, const double *theta,
                                   const double *phi,   const double *psi,
                                   const double *dgamma, const double *dtheta,
                                   const double *dphi,   const double *dpsi,
                                   BINT n)
{
    BINT i; double alpha = 1.0;
    for (i = 0; i < n; ++i) {
        if (dgamma[i] < 0.0) { double c = -0.99*gamma[i]/dgamma[i]; if (c < alpha) alpha = c; }
        if (dtheta[i] < 0.0) { double c = -0.99*theta[i]/dtheta[i]; if (c < alpha) alpha = c; }
        if (dphi[i]   < 0.0) { double c = -0.99*phi[i]  /dphi[i];   if (c < alpha) alpha = c; }
        if (dpsi[i]   < 0.0) { double c = -0.99*psi[i]  /dpsi[i];   if (c < alpha) alpha = c; }
    }
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    return alpha;
}

/* The reduced KKT system after eliminating gamma, theta, phi, psi is the dense linear system:(H + diag(gamma./phi + theta./psi)) dz = r2_over_psi - r1_over_phi
 where r1_over_phi and r2_over_psi already carry the appropriate RHS for the predictor (-gamma, -theta) or corrector (centering + 2nd-order term)*/
static int factor_dense_system(const double *H,
                               const double *gamma, const double *theta,
                               const double *phi,   const double *psi,
                               BINT n,
                               double *Hbar,
                               BINT *chol_info)
{
    BINT i, j, info = 0;

    /* Build the lower triangle of Hbar = H + diag(gamma/phi + theta/psi). dpotrf("l",...) and dpotrs("l",...) ignore the upper triangle. */
    for (j = 0; j < n; ++j) {
        for (i = j; i < n; ++i)
            Hbar[i + j*n] = H[i + j*n];
        Hbar[j + j*n] += gamma[j]/phi[j] + theta[j]/psi[j];
    }

    dpotrf_("l", &n, Hbar, &n, &info);
    *chol_info = info;
    if (info != 0) return -1;

    return 0;
}

static int solve_factored_direction(const double *Hbar,
                                    const double *gamma, const double *theta,
                                    const double *phi,   const double *psi,
                                    const double *r1_over_phi,
                                    const double *r2_over_psi,
                                    BINT n,
                                    double *rhs,
                                    double *dz, double *dgamma, double *dtheta,
                                    double *dphi, double *dpsi,
                                    BINT *chol_info)
{
    BINT i, nrhs = 1, info = 0;

    for (i = 0; i < n; ++i)
        rhs[i] = r2_over_psi[i] - r1_over_phi[i];

    dpotrs_("l", &n, &nrhs, Hbar, &n, rhs, &n, &info);
    *chol_info = info;
    if (info != 0) return -1;

    /* dz solved; recover the eliminated directions. */
    dcopy_(&n, rhs, &IONE, dz, &IONE);
    for (i = 0; i < n; ++i) {
        dphi[i]   = -dz[i];
        dpsi[i]   =  dz[i];
        dgamma[i] = (gamma[i]/phi[i])*dz[i] + r1_over_phi[i];
        dtheta[i] = -(theta[i]/psi[i])*dz[i] + r2_over_psi[i];
    }
    return 0;
}

/* dense_boxqp_solve(Algorithm 1):
 H        : n x n dense SPD-ish reduced Hessian (column-major)
 h        : n linear term
 epsilon  : duality-gap tolerance 
 max_iter : iteration cap 
 z_warm   : optional warm start (NULL for cold start at 0); clipped into (-1,1)
 z        : output, length n
 info_out : optional diagnostics (NULL allowed)
 Returns: status (0 converged, 1 max-iter, -1 Cholesky failure, -2 alloc).*/
static int dense_boxqp_solve_workspace(const double *H, const double *h,
                                       BINT n, double epsilon, BINT max_iter,
                                       const double *z_warm, double *z,
                                       BoxQPInfo *info_out,
                                       BoxQPWorkspace *work)
{
    BINT i, iter = 0, chol_info = 0;
    int status = 1;
    double duality = 0.0, mu = INFINITY;

    double *gamma, *theta, *phi, *psi, *grad, *Hbar, *rhs, *r1, *r2;
    double *dz_a, *dg_a, *dt_a, *dp_a, *dq_a, *dz, *dg, *dt, *dp, *dq;

    if (boxqp_workspace_prepare(work, n) != 0) {
        status = -2;
        goto finish;
    }

    gamma = work->gamma; theta = work->theta; phi = work->phi; psi = work->psi;
    grad = work->grad; Hbar = work->Hbar; rhs = work->rhs; r1 = work->r1; r2 = work->r2;
    dz_a = work->dz_a; dg_a = work->dg_a; dt_a = work->dt_a; dp_a = work->dp_a; dq_a = work->dq_a;
    dz = work->dz; dg = work->dg; dt = work->dt; dp = work->dp; dq = work->dq;

    for (i = 0; i < n; ++i) {
        z[i] = z_warm ? clip_box_interior(z_warm[i]) : 0.0;
        phi[i] = 1.0 - z[i];
        psi[i] = 1.0 + z[i];
    }

    /* Strictly-feasible dual start: choose gamma, theta so that H*z + h + gamma - theta = 0 exactly and gamma, theta > 0. */
    compute_grad(H, h, z, n, grad);
    {
        double eta = max_abs_vec(grad, n);
        if (eta < 1.0) eta = 1.0;
        for (i = 0; i < n; ++i) {
            gamma[i] = eta - 0.5*grad[i];
            theta[i] = eta + 0.5*grad[i];
        }
    }

    /* ---- IPM main loop ---- */
    for (iter = 0; iter < max_iter; ++iter) {
        double alpha_aff, alpha, duality_aff = 0.0, mu_aff, sigma, sigma_mu;

        duality = ddot_(&n, gamma, &IONE, phi, &IONE)
                + ddot_(&n, theta, &IONE, psi, &IONE);
        mu = duality / (2.0*(double)n);
        if (mu <= epsilon) { status = 0; break; }

        if (factor_dense_system(H, gamma, theta, phi, psi, n, Hbar,
                                &chol_info) != 0) { status = -1; break; }

        /* Predictor (affine): complementarity RHS = (-gamma.*phi, -theta.*psi),
           pre-divided form r1 = -gamma, r2 = -theta. */
        for (i = 0; i < n; ++i) { r1[i] = -gamma[i]; r2[i] = -theta[i]; }

        if (solve_factored_direction(Hbar, gamma, theta, phi, psi, r1, r2, n,
                                     rhs, dz_a, dg_a, dt_a, dp_a, dq_a,
                                     &chol_info) != 0) { status = -1; break; }

        alpha_aff = fraction_to_boundary(gamma, theta, phi, psi,
                                         dg_a, dt_a, dp_a, dq_a, n);

        for (i = 0; i < n; ++i) {
            double ga = gamma[i] + alpha_aff*dg_a[i];
            double th = theta[i] + alpha_aff*dt_a[i];
            double ph = phi[i]   + alpha_aff*dp_a[i];
            double ps = psi[i]   + alpha_aff*dq_a[i];
            duality_aff += ga*ph + th*ps;
        }
        mu_aff = duality_aff / (2.0*(double)n);

        /* Centering parameter (Mehrotra). */
        sigma = mu_aff / mu;
        sigma = sigma*sigma*sigma;
        if (!(sigma >= 0.0)) sigma = 0.0;     /* NaN-safe */
        if (sigma > 1.0)     sigma = 1.0;
        sigma_mu = sigma*mu;

        /* Corrector: centering + Mehrotra 2nd-order term. */
        for (i = 0; i < n; ++i) {
            r1[i] = -gamma[i] + (sigma_mu - dg_a[i]*dp_a[i]) / phi[i];
            r2[i] = -theta[i] + (sigma_mu - dt_a[i]*dq_a[i]) / psi[i];
        }

        if (solve_factored_direction(Hbar, gamma, theta, phi, psi, r1, r2, n,
                                     rhs, dz, dg, dt, dp, dq,
                                     &chol_info) != 0) { status = -1; break; }

        alpha = fraction_to_boundary(gamma, theta, phi, psi, dg, dt, dp, dq, n);

        /* Commit the step. */
        for (i = 0; i < n; ++i) {
            z[i]     += alpha*dz[i];
            gamma[i] += alpha*dg[i];
            theta[i] += alpha*dt[i];
            phi[i]   += alpha*dp[i];
            psi[i]   += alpha*dq[i];
        }
    }

finish:
    if (info_out) {
        info_out->iterations = iter;
        info_out->mu = mu;
        info_out->duality = duality;
        info_out->status = status;
        info_out->chol_info = chol_info;
    }
    return status;
}

int dense_boxqp_solve(const double *H, const double *h,
                      BINT n, double epsilon, BINT max_iter,
                      const double *z_warm, double *z,
                      BoxQPInfo *info_out)
{
    return dense_boxqp_solve_workspace(H, h, n, epsilon, max_iter,
                                       z_warm, z, info_out,
                                       &g_boxqp_workspace);
}

/* Generic Woodbury/K-form reduction offline */
typedef struct {
    BINT nfixed, nz, ncol, nrows;
    double ridge;
    double *A;              /* [A_fixed; A_decision], nrows x ncol */
    double *Gchol;          /* chol(A A' / ridge + diag(1 / weights)) */
    double *Hphys;          /* physical-variable reduced Hessian */
    double *Fphys;          /* fixed_data -> physical linear term */
    double *Hbox;           /* unit-box Hessian */
    double *Fbox;           /* fixed_data -> unit-box linear term */
    double *center, *scale, *hcenter, *hbox;
    double *swarm, *sout;
} DeePCKForm;

void deepc_kform_free(DeePCKForm *d)
{
    if (!d) return;
    free(d->A); free(d->Gchol); free(d->Hphys); free(d->Fphys);
    free(d->Hbox); free(d->Fbox); free(d->center); free(d->scale);
    free(d->hcenter); free(d->hbox); free(d->swarm); free(d->sout);
    memset(d, 0, sizeof(*d));
}

/*Eliminate g with the Woodbury identity.  */
int deepc_kform_offline(DeePCKForm *d,
                        const double *A_fixed, const double *A_decision,
                        BINT nfixed, BINT nz, BINT ncol,
                        const double *fixed_weights,
                        const double *decision_weights,
                        double ridge, const double *decision_cost,
                        const double *lower, const double *upper)
{
    BINT i, j, info = 0, nrhs;
    BINT nrows = nfixed + nz;
    double inv_ridge;
    double *solve_z = NULL;

    if (!d || !A_fixed || !A_decision || !fixed_weights ||
        !decision_weights || !decision_cost || !lower || !upper ||
        nfixed <= 0 || nz <= 0 || ncol <= 0 || ridge <= 0.0) return -10;

    memset(d, 0, sizeof(*d));
    d->nfixed = nfixed; d->nz = nz; d->ncol = ncol; d->nrows = nrows;
    d->ridge = ridge;
    d->A       = calloc((size_t)nrows*ncol, sizeof(double));
    d->Gchol   = calloc((size_t)nrows*nrows, sizeof(double));
    d->Hphys   = calloc((size_t)nz*nz, sizeof(double));
    d->Fphys   = calloc((size_t)nz*nfixed, sizeof(double));
    d->Hbox    = calloc((size_t)nz*nz, sizeof(double));
    d->Fbox    = calloc((size_t)nz*nfixed, sizeof(double));
    d->center  = calloc((size_t)nz, sizeof(double));
    d->scale   = calloc((size_t)nz, sizeof(double));
    d->hcenter = calloc((size_t)nz, sizeof(double));
    d->hbox    = calloc((size_t)nz, sizeof(double));
    d->swarm   = calloc((size_t)nz, sizeof(double));
    d->sout    = calloc((size_t)nz, sizeof(double));
    solve_z    = calloc((size_t)nrows*nz, sizeof(double));
    if (!d->A || !d->Gchol || !d->Hphys || !d->Fphys || !d->Hbox ||
        !d->Fbox || !d->center || !d->scale || !d->hcenter || !d->hbox ||
        !d->swarm || !d->sout || !solve_z) { info = -1; goto done; }

    for (j = 0; j < ncol; ++j) {
        for (i = 0; i < nfixed; ++i)
            d->A[i + j*nrows] = A_fixed[i + j*nfixed];
        for (i = 0; i < nz; ++i)
            d->A[nfixed + i + j*nrows] = A_decision[i + j*nz];
    }

    inv_ridge = 1.0 / ridge;
    dgemm_("n", "t", &nrows, &nrows, &ncol, &inv_ridge,
           d->A, &nrows, d->A, &nrows, &DZERO, d->Gchol, &nrows);
    for (i = 0; i < nfixed; ++i) {
        if (fixed_weights[i] <= 0.0) { info = -11; goto done; }
        d->Gchol[i + i*nrows] += 1.0 / fixed_weights[i];
    }
    for (i = 0; i < nz; ++i) {
        if (decision_weights[i] <= 0.0 || upper[i] <= lower[i]) {
            info = -12; goto done;
        }
        d->Gchol[(nfixed+i) + (nfixed+i)*nrows] += 1.0 / decision_weights[i];
    }
    dpotrf_("l", &nrows, d->Gchol, &nrows, &info);
    if (info != 0) { info = -2; goto done; }

    /* solve_z = K[:, decision rows], K=(A A'/ridge+W^{-1})^{-1}. */
    for (j = 0; j < nz; ++j) solve_z[(nfixed+j) + j*nrows] = 1.0;
    nrhs = nz;
    dpotrs_("l", &nrows, &nrhs, d->Gchol, &nrows,
            solve_z, &nrows, &info);
    if (info != 0) { info = -3; goto done; }

    for (j = 0; j < nfixed; ++j)
        for (i = 0; i < nz; ++i)
            d->Fphys[i + j*nz] = solve_z[j + i*nrows];
    for (j = 0; j < nz; ++j)
        for (i = 0; i < nz; ++i)
            d->Hphys[i + j*nz] = decision_cost[i + j*nz]
                                   + solve_z[(nfixed+i) + j*nrows];

    for (i = 0; i < nz; ++i) {
        d->center[i] = 0.5 * (lower[i] + upper[i]);
        d->scale[i] = 0.5 * (upper[i] - lower[i]);
    }
    for (j = 0; j < nz; ++j)
        for (i = 0; i < nz; ++i)
            d->Hbox[i + j*nz] = d->scale[i] * d->Hphys[i + j*nz] * d->scale[j];
    for (j = 0; j < nfixed; ++j)
        for (i = 0; i < nz; ++i)
            d->Fbox[i + j*nz] = d->scale[i] * d->Fphys[i + j*nz];
    for (i = 0; i < nz; ++i) {
        double value = 0.0;
        for (j = 0; j < nz; ++j) value += d->Hphys[i + j*nz] * d->center[j];
        d->hcenter[i] = d->scale[i] * value;
    }

    info = 0;
done:
    free(solve_z);
    if (info != 0) deepc_kform_free(d);
    return (int)info;
}

/* Fixed-size online map and unit-box solve; z_warm/z_out use physical units. */
int deepc_kform_online_step(DeePCKForm *d,
                            const double *fixed_data,
                            const double *tracking_linear,
                            const double *z_warm,
                            double epsilon, BINT max_iter,
                            double *z_out, BoxQPInfo *info_out)
{
    BINT i;
    int status;
    if (!d || !fixed_data || !tracking_linear || !z_out) return -10;
    dcopy_(&d->nz, d->hcenter, &IONE, d->hbox, &IONE);
    dgemv_("n", &d->nz, &d->nfixed, &DONE, d->Fbox, &d->nz,
           fixed_data, &IONE, &DONE, d->hbox, &IONE);
    for (i = 0; i < d->nz; ++i)
        d->hbox[i] += d->scale[i] * tracking_linear[i];
    if (z_warm) {
        for (i = 0; i < d->nz; ++i)
            d->swarm[i] = clip_box_interior((z_warm[i] - d->center[i]) / d->scale[i]);
    }
    status = dense_boxqp_solve(d->Hbox, d->hbox, d->nz, epsilon, max_iter,
                               z_warm ? d->swarm : NULL, d->sout, info_out);
    for (i = 0; i < d->nz; ++i)
        z_out[i] = d->center[i] + d->scale[i] * d->sout[i];
    return status;
}

/* Recover the eliminated coordinate only for diagnostics and full-objective evaluation.*/
int deepc_kform_recover_g(DeePCKForm *d,
                          const double *fixed_data,
                          const double *z,
                          double *g_out)
{
    BINT i, info = 0, nrhs = 1;
    double alpha;
    double *rhs;

    if (!d || !fixed_data || !z || !g_out || !d->A || !d->Gchol ||
        d->ridge <= 0.0) return -10;

    rhs = calloc((size_t)d->nrows, sizeof(double));
    if (!rhs) return -1;

    for (i = 0; i < d->nfixed; ++i) rhs[i] = fixed_data[i];
    for (i = 0; i < d->nz; ++i) rhs[d->nfixed + i] = z[i];

    dpotrs_("l", &d->nrows, &nrhs, d->Gchol, &d->nrows,
            rhs, &d->nrows, &info);
    if (info != 0) {
        free(rhs);
        return -2;
    }

    alpha = 1.0 / d->ridge;
    dgemv_("t", &d->nrows, &d->ncol, &alpha,
           d->A, &d->nrows, rhs, &IONE, &DZERO, g_out, &IONE);

    free(rhs);
    return 0;
}

/* DeePC OFFLINE / ONLINE ASSEMBLY   
 m      : number of inputs
 p      : number of outputs
 Tini   : length of the initial (past) window
 N      : prediction horizon
 L      = Tini + N           (Hankel column height multiplier)
 Tdata  : length of the offline data sequences ud, yd
 ncol   = Tdata - L + 1      (number of Hankel columns)
 Up : (m*Tini) x ncol      Uf : (m*N) x ncol
 Yp : (p*Tini) x ncol      Yf : (p*N) x ncol
 Z  = [Up; Yp; Uf]         : (m*Tini + p*Tini + m*N) x ncol
 M  : ncol x ncol          (elimination matrix, Cholesky-factored)
 nz = m*N + p*N            (reduced BoxQP dimension, z = col(u,y))
 Hhat : nz x nz dense reduced Hessian*/

typedef struct {
    BINT m, p, Tini, N, L, ncol, nz, nu, ny;
    /* Hankel blocks (column-major) */
    double *Up, *Yp, *Uf, *Yf;     /* future/past partitions */
    double *Z;                     /* [Up; Yp; Uf] stacked, rows = nZrows */
    BINT    nZrows;                /* m*Tini + p*Tini + m*N */
    /* elimination matrix and its Cholesky factor (lower) */
    double *Mchol;                 /* ncol x ncol, factored in place */
    /* reduced Hessian */
    double *Hhat;                  /* nz x nz */
    /* weights and penalties */
    double *Wu, *Wy;               /* m*N x m*N , p*N x p*N (block-diagonal ok) */
    double  lambda, rho, lambda_g;
    /* T-independent online linear-term maps and small past window */
    double *Su;                    /* m*N x Tini*(m+p) */
    double *Sy;                    /* p*N x Tini*(m+p) */
    double *wpast;                 /* Tini*(m+p) */
    double *hhat;                  /* nz */
} DeePC;

/* Build a block-Hankel matrix of column-height 'block_rows = dim*L' from a signal s of length Tdata with 'dim' channels (column-major signal: s[k*dim + c]
   is channel c at time k).  Output Hank is (dim*L) x ncol, column-major.
   Hankel column j, row-block r (r = 0..L-1), channel c:Hank[(r*dim + c) + j*(dim*L)] = s[(j + r)*dim + c]*/
static void build_hankel(const double *s, BINT dim, BINT Tdata, BINT L,
                         BINT ncol, double *Hank)
{
    BINT j, r, c;
    BINT block_rows = dim*L;
    for (j = 0; j < ncol; ++j)
        for (r = 0; r < L; ++r)
            for (c = 0; c < dim; ++c)
                Hank[(r*dim + c) + j*block_rows] = s[(j + r)*dim + c];
}

/*This API is an older version of offline construction and it is not used in the current implementation.*/
int deepc_offline(DeePC *d,
                  const double *ud, const double *yd, BINT Tdata,
                  BINT m, BINT p, BINT Tini, BINT N,
                  const double *Wu, const double *Wy,
                  double lambda, double rho, double lambda_g)
{
    BINT info = 0, i, j;
    BINT L     = Tini + N;
    BINT ncol  = Tdata - L + 1;
    BINT mTini = m*Tini, pTini = p*Tini, mN = m*N, pN = p*N;
    BINT nPast = mTini + pTini;
    BINT nZrows = mTini + pTini + mN;
    BINT nz = mN + pN;

    double *Uhank = NULL, *Yhank = NULL;
    double *M = NULL, *Minv_Uf = NULL, *Minv_Yf = NULL, *Ppast = NULL;
    BINT nrhs;

    if (ncol <= 0) return -10;

    memset(d, 0, sizeof(*d));
    d->m=m; d->p=p; d->Tini=Tini; d->N=N; d->L=L; d->ncol=ncol;
    d->nz=nz; d->nu=mN; d->ny=pN; d->nZrows=nZrows;
    d->lambda=lambda; d->rho=rho; d->lambda_g=lambda_g;
    
    d->Up=calloc((size_t)mTini*ncol,sizeof(double));
    d->Yp=calloc((size_t)pTini*ncol,sizeof(double));
    d->Uf=calloc((size_t)mN*ncol,sizeof(double));
    d->Yf=calloc((size_t)pN*ncol,sizeof(double));
    d->Z =calloc((size_t)nZrows*ncol,sizeof(double));
    d->Mchol=calloc((size_t)ncol*ncol,sizeof(double));
    d->Hhat =calloc((size_t)nz*nz,sizeof(double));
    d->Wu=calloc((size_t)mN*mN,sizeof(double));
    d->Wy=calloc((size_t)pN*pN,sizeof(double));
    d->Su=calloc((size_t)mN*nPast,sizeof(double));
    d->Sy=calloc((size_t)pN*nPast,sizeof(double));
    d->wpast=calloc((size_t)nPast,sizeof(double));
    d->hhat=calloc((size_t)nz,sizeof(double));

    /* temporaries */
    Uhank   = calloc((size_t)m*L*ncol,sizeof(double));
    Yhank   = calloc((size_t)p*L*ncol,sizeof(double));
    M       = calloc((size_t)ncol*ncol,sizeof(double));
    Minv_Uf = calloc((size_t)ncol*mN,sizeof(double));  /* holds M^{-1} Uf'  (ncol x mN) */
    Minv_Yf = calloc((size_t)ncol*pN,sizeof(double));  /* holds M^{-1} Yf'  (ncol x pN) */
    Ppast    = calloc((size_t)ncol*nPast,sizeof(double)); /* M^{-1} lambda [Up;Yp]' */

    if(!d->Up||!d->Yp||!d->Uf||!d->Yf||!d->Z||!d->Mchol||!d->Hhat||!d->Wu||!d->Wy||
       !d->Su||!d->Sy||!d->wpast||!d->hhat||!Uhank||!Yhank||!M||!Minv_Uf||!Minv_Yf||!Ppast){
        info=-1; goto done;
    }

    memcpy(d->Wu, Wu, sizeof(double)*mN*mN);
    memcpy(d->Wy, Wy, sizeof(double)*pN*pN);

    /* --- Build full Hankel matrices, then split into past/future --- */
    build_hankel(ud, m, Tdata, L, ncol, Uhank);   /* (m*L) x ncol */
    build_hankel(yd, p, Tdata, L, ncol, Yhank);   /* (p*L) x ncol */

    /* Up = first m*Tini rows of Uhank ; Uf = next m*N rows */
    for (j = 0; j < ncol; ++j) {
        for (i = 0; i < mTini; ++i) d->Up[i + j*mTini] = Uhank[i + j*(m*L)];
        for (i = 0; i < mN;    ++i) d->Uf[i + j*mN]    = Uhank[(mTini + i) + j*(m*L)];
    }
    /* Yp = first p*Tini rows of Yhank ; Yf = next p*N rows */
    for (j = 0; j < ncol; ++j) {
        for (i = 0; i < pTini; ++i) d->Yp[i + j*pTini] = Yhank[i + j*(p*L)];
        for (i = 0; i < pN;    ++i) d->Yf[i + j*pN]    = Yhank[(pTini + i) + j*(p*L)];
    }

    /* --- Z = [Up; Yp; Uf]  (nZrows x ncol) --- */
    for (j = 0; j < ncol; ++j) {
        BINT off = 0;
        for (i = 0; i < mTini; ++i) d->Z[off + i + j*nZrows] = d->Up[i + j*mTini];
        off += mTini;
        for (i = 0; i < pTini; ++i) d->Z[off + i + j*nZrows] = d->Yp[i + j*pTini];
        off += pTini;
        for (i = 0; i < mN;    ++i) d->Z[off + i + j*nZrows] = d->Uf[i + j*mN];
    }

    /* --- M = lambda Z'Z + rho Yf'Yf + lambda_g I   (ncol x ncol) --- */
    /* M = lambda * Z'Z */
    dgemm_("t","n", &ncol,&ncol,&nZrows, &lambda, d->Z,&nZrows, d->Z,&nZrows,
           &DZERO, M, &ncol);
    /* M += rho * Yf'Yf */
    dgemm_("t","n", &ncol,&ncol,&pN, &rho, d->Yf,&pN, d->Yf,&pN, &DONE, M, &ncol);
    /* M += lambda_g I */
    for (i = 0; i < ncol; ++i) M[i + i*ncol] += lambda_g;

    /* Cholesky factor M (store factor in d->Mchol; keep M intact for solves) */
    memcpy(d->Mchol, M, sizeof(double)*ncol*ncol);
    dpotrf_("l", &ncol, d->Mchol, &ncol, &info);
    if (info != 0) { info = -2; goto done; }

    /* --- M^{-1} Uf'  and  M^{-1} Yf'  via the stored factor --- */
    /* Minv_Uf = Uf' (ncol x mN), then solve in place */
    for (j = 0; j < mN; ++j)
        for (i = 0; i < ncol; ++i)
            Minv_Uf[i + j*ncol] = d->Uf[j + i*mN];   /* (Uf')[i,j] = Uf[j,i] */
    nrhs = mN;
    dpotrs_("l", &ncol, &nrhs, d->Mchol, &ncol, Minv_Uf, &ncol, &info);
    if (info != 0) { info = -3; goto done; }

    for (j = 0; j < pN; ++j)
        for (i = 0; i < ncol; ++i)
            Minv_Yf[i + j*ncol] = d->Yf[j + i*pN];   /* (Yf')[i,j] = Yf[j,i] */
    nrhs = pN;
    dpotrs_("l", &ncol, &nrhs, d->Mchol, &ncol, Minv_Yf, &ncol, &info);
    if (info != 0) { info = -4; goto done; }

    /* --- Direct online linear-term maps ---
       Only the past rows of wini are nonzero, so precompute

           Ppast = M^{-1} lambda [Up;Yp]'
           Su    = -lambda Uf Ppast
           Sy    = -rho    Yf Ppast.

       This moves every ncol-dependent operation offline. */
    for (j = 0; j < nPast; ++j)
        for (i = 0; i < ncol; ++i)
            Ppast[i + j*ncol] = lambda * d->Z[j + i*nZrows];
    nrhs = nPast;
    dpotrs_("l", &ncol, &nrhs, d->Mchol, &ncol, Ppast, &ncol, &info);
    if (info != 0) { info = -5; goto done; }

    {
        double neg_lambda = -lambda;
        double neg_rho = -rho;
        dgemm_("n","n", &mN,&nPast,&ncol, &neg_lambda,
               d->Uf,&mN, Ppast,&ncol, &DZERO, d->Su,&mN);
        dgemm_("n","n", &pN,&nPast,&ncol, &neg_rho,
               d->Yf,&pN, Ppast,&ncol, &DZERO, d->Sy,&pN);
    }

    /* --- Reduced Hessian blocks ---
       Huu = Wu + lambda I  - lambda^2  Uf (M^{-1} Uf')
       Hyy = Wy + rho    I  - rho^2     Yf (M^{-1} Yf')
       Huy =               - lambda rho  Uf (M^{-1} Yf')
       Assemble directly into Hhat (nz x nz), z = col(u,y), u-block first.   */
    {
        double l2  = -lambda*lambda;
        double r2  = -rho*rho;
        double lr  = -lambda*rho;
        BINT   nz_ = nz;

        /* Huu into top-left (mN x mN): start from Wu + lambda I */
        for (j = 0; j < mN; ++j)
            for (i = 0; i < mN; ++i)
                d->Hhat[i + j*nz_] = d->Wu[i + j*mN] + ((i==j)? lambda : 0.0);
        /* top-left -= lambda^2 Uf*(M^{-1}Uf') :  Uf (mN x ncol) * Minv_Uf (ncol x mN) */
        dgemm_("n","n", &mN,&mN,&ncol, &l2, d->Uf,&mN, Minv_Uf,&ncol,
               &DONE, d->Hhat, &nz_);

        /* Hyy into bottom-right (pN x pN): start from Wy + rho I */
        for (j = 0; j < pN; ++j)
            for (i = 0; i < pN; ++i)
                d->Hhat[(mN+i) + (mN+j)*nz_] = d->Wy[i + j*pN] + ((i==j)? rho : 0.0);
        /* bottom-right -= rho^2 Yf*(M^{-1}Yf') */
        dgemm_("n","n", &pN,&pN,&ncol, &r2, d->Yf,&pN, Minv_Yf,&ncol,
               &DONE, d->Hhat + mN + mN*nz_, &nz_);

        /* Huy into top-right (mN x pN) = -lambda rho Uf*(M^{-1}Yf') */
        dgemm_("n","n", &mN,&pN,&ncol, &lr, d->Uf,&mN, Minv_Yf,&ncol,
               &DZERO, d->Hhat + 0 + mN*nz_, &nz_);

        /* Huy' into bottom-left (pN x mN): transpose of the block above */
        for (j = 0; j < mN; ++j)
            for (i = 0; i < pN; ++i)
                d->Hhat[(mN+i) + j*nz_] = d->Hhat[j + (mN+i)*nz_];
    }

    info = 0;
done:
    free(Uhank); free(Yhank); free(M); free(Minv_Uf); free(Minv_Yf); free(Ppast);
    return (int)info;
}

/* Shift the previous optimizer forward by one MPC step, so [u0,u1,...,uN-1,y0,y1,...,yN-1]becomes [u1,...,uN-1,uN-1,y1,...,yN-1,yN-1]. */
void deepc_shift_warm_start(const DeePC *d, const double *z_prev, double *z_shift)
{
    BINT k, c;
    BINT m = d->m, p = d->p, N = d->N, mN = d->nu;
    if (!d || !z_prev || !z_shift) return;

    for (k = 0; k < N; ++k) {
        BINT src_k = (k + 1 < N) ? (k + 1) : (N - 1);
        for (c = 0; c < m; ++c)
            z_shift[k*m + c] = clip_box_interior(z_prev[src_k*m + c]);
        for (c = 0; c < p; ++c)
            z_shift[mN + k*p + c] = clip_box_interior(z_prev[mN + src_k*p + c]);
    }
}

/*deepc_online_step that solve reduced BoxQP  min 0.5 z'Hhat z + hhat' z  s.t. -1<=z<=1
 wpast = [uini; yini]
 hu    = Su wpast
 hy    = -Wy rt + Sy wpast
 hhat = [hu; hy]

 uini : length m*Tini   (most recent Tini inputs, stacked col-major)
 yini : length p*Tini   (most recent Tini outputs)
 rt   : length p*N      (output reference over the horizon)
 z_warm : optional warm start (length nz) or NULL
 z_out  : output, length nz = m*N + p*N, = col(u_opt, y_opt)
 The first optimal input u0 is z_out[0 .. m-1, and returns BoxQP status (0 converged, 1 max-iter, -1 chol fail, -2 alloc).*/
int deepc_online_step(DeePC *d,
                      const double *uini, const double *yini, const double *rt,
                      const double *z_warm,
                      double epsilon, BINT max_iter,
                      double *z_out, BoxQPInfo *info_out)
{
    BINT i;
    BINT mTini = d->m*d->Tini, pTini = d->p*d->Tini, mN = d->nu, pN = d->ny;
    BINT nPast = mTini + pTini;

    /* wpast = [uini; yini] */
    for (i = 0; i < mTini; ++i) d->wpast[i] = uini[i];
    for (i = 0; i < pTini; ++i) d->wpast[mTini + i] = yini[i];

    /* hu = Su wpast: dimensions depend only on N,m,p,Tini. */
    dgemv_("n", &mN, &nPast, &DONE, d->Su, &mN,
           d->wpast, &IONE, &DZERO, d->hhat, &IONE);

    /* hy = -Wy rt + Sy wpast  (write into hhat + mN) */
    {
        /* hy = -Wy rt */
        dgemv_("n", &pN, &pN, &DMONE, d->Wy, &pN, rt, &IONE, &DZERO,
               d->hhat + mN, &IONE);
        /* hy += Sy wpast */
        dgemv_("n", &pN, &nPast, &DONE, d->Sy, &pN, d->wpast, &IONE, &DONE,
               d->hhat + mN, &IONE);
    }

    /* Solve the reduced BoxQP */
    return dense_boxqp_solve(d->Hhat, d->hhat, d->nz, epsilon, max_iter,
                             z_warm, z_out, info_out);
}

void deepc_free(DeePC *d)
{
    if(!d) return;
    free(d->Up); free(d->Yp); free(d->Uf); free(d->Yf); free(d->Z);
    free(d->Mchol); free(d->Hhat); free(d->Wu); free(d->Wy);
    free(d->Su); free(d->Sy); free(d->wpast); free(d->hhat);
    deepc_boxqp_release_workspace();
    memset(d, 0, sizeof(*d));
}
