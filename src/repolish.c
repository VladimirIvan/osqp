#include "repolish.h"
#include "polish.h"
#include "lin_alg.h"
#include "util.h"
#include "auxil.h"
#include "lin_sys.h"
#include "kkt.h"
#include "proj.h"
#include "error.h"

/**
 * Allocates reduced matrix A, that contains only rows that are active at the
 * solution. Also allocates corresponding right hand side, reduced right
 * hand side, linear solver and polished solution vector.
 * Ared = vstack[Alow, Aupp]
 * Active constraints are guessed from the primal and dual solution returned by
 * the ADMM.
 * @param  work Workspace
 * @param  data Repolish data
 * @return      Number of rows in Ared, negative if error
 */
static c_int form_Ared(OSQPWorkspace *work, OSQPRepolish *data) {
  c_int j, ptr, n;
  c_int Ared_nnz = 0;

  // Initialize counters for active constraints
  work->pol->n_low = 0;
  work->pol->n_upp = 0;

  /* Guess which linear constraints are lower-active, upper-active and free
   *    A_to_Alow[j] = -1    (if j-th row of A is not inserted in Alow)
   *    A_to_Alow[j] =  i    (if j-th row of A is inserted at i-th row of Alow)
   * Aupp is formed in the equivalent way.
   * Ared is formed by stacking vertically Alow and Aupp.
   */
  for (j = 0; j < work->data->m; j++) {
    if (work->z[j] - work->data->l[j] < -work->y[j]) { // lower-active
      work->pol->Alow_to_A[work->pol->n_low] = j;
      work->pol->A_to_Alow[j]                = work->pol->n_low++;
    } else {
      work->pol->A_to_Alow[j] = -1;
    }
  }

  for (j = 0; j < work->data->m; j++) {
    if (work->data->u[j] - work->z[j] < work->y[j]) { // upper-active
      work->pol->Aupp_to_A[work->pol->n_upp] = j;
      work->pol->A_to_Aupp[j]                = work->pol->n_upp++;
    } else {
      work->pol->A_to_Aupp[j] = -1;
    }
  }

  // Check if there are no active constraints
  if (work->pol->n_low + work->pol->n_upp == 0) {
    // Form empty Ared
    work->pol->Ared = csc_spalloc(0, work->data->n, 0, 1, 0);
    if (!(work->pol->Ared)) return -1;
    int_vec_set_scalar(work->pol->Ared->p, 0, work->data->n + 1);
    return 0; // mred = 0
  }

  // Count number of elements in Ared
  for (j = 0; j < work->data->A->p[work->data->A->n]; j++) {
    if ((work->pol->A_to_Alow[work->data->A->i[j]] != -1) ||
        (work->pol->A_to_Aupp[work->data->A->i[j]] != -1)) Ared_nnz++;
  }

  // Form Ared
  // Ared = vstack[Alow, Aupp]
  work->pol->Ared = csc_spalloc(work->pol->n_low + work->pol->n_upp,
                                work->data->n, Ared_nnz, 1, 0);
  if (!(work->pol->Ared)) return -1;

  n = work->data->n + work->pol->Ared->m;

  // Allocate rhs vector
  data->rhs = (c_float *)c_malloc(sizeof(c_float) * n);

  if (!data->rhs) {
    csc_spfree(work->pol->Ared);
    return osqp_error(OSQP_MEM_ALLOC_ERROR);
  }

  // Return number of rows in Ared
  return work->pol->n_low + work->pol->n_upp;
}

/**
 * Perform iterative refinement on the polished solution:
 *    (repeat)
 *    1. (K + dK) * dz = b - K*z
 *    2. z <- z + dz
 * @param  work Solver workspace
 * @param  data Repolish data
 * @param  z    Initial z value
 * @param  b    RHS of the linear system
 * @return      Exitflag
 */
static c_int iterative_refinement(OSQPWorkspace *work,
                                  OSQPRepolish  *data,
                                  c_float       *z,
                                  c_float       *b) {
  c_int i, j, n;
  c_float *rhs = data->rhs;
  LinSysSolver  *p = data->plsh;

  if (work->settings->polish_refine_iter > 0) {

    // Assign dimension n
    n = work->data->n + work->pol->Ared->m;

    for (i = 0; i < work->settings->polish_refine_iter; i++) {
      // Form the RHS for the iterative refinement:  b - K*z
      prea_vec_copy(b, rhs, n);

      // Upper Part: R^{n}
      // -= Px (upper triang)
      mat_vec(work->data->P, z, rhs, -1);

      // -= Px (lower triang)
      mat_tpose_vec(work->data->P, z, rhs, -1, 1);

      // -= Ared'*y_red
      mat_tpose_vec(work->pol->Ared, z + work->data->n, rhs, -1, 0);

      // Lower Part: R^{m}
      mat_vec(work->pol->Ared, z, rhs + work->data->n, -1);

      // Solve linear system. Store solution in rhs
      p->solve(p, rhs);

      // Update solution
      for (j = 0; j < n; j++) {
        z[j] += rhs[j];
      }
    }
  }
  return 0;
}

/**
 * @brief Fills Ared data from workspace
 * 
 * @param  work Workspace
 * @param  data Repolish data
 */
static void fill_Ared(OSQPWorkspace *work, OSQPRepolish *data) {
  c_int Ared_nnz = 0; // counter
  c_int j, ptr;

  for (j = 0; j < work->data->n; j++) { // Cycle over columns of A
    work->pol->Ared->p[j] = Ared_nnz;

    for (ptr = work->data->A->p[j]; ptr < work->data->A->p[j + 1]; ptr++) {
      // Cycle over elements in j-th column
      if (work->pol->A_to_Alow[work->data->A->i[ptr]] != -1) {
        // Lower-active rows of A
        work->pol->Ared->i[Ared_nnz] =
          work->pol->A_to_Alow[work->data->A->i[ptr]];
        work->pol->Ared->x[Ared_nnz++] = work->data->A->x[ptr];
      } else if (work->pol->A_to_Aupp[work->data->A->i[ptr]] != -1) {
        // Upper-active rows of A
        work->pol->Ared->i[Ared_nnz] = work->pol->A_to_Aupp[work->data->A->i[ptr]] \
                                       + work->pol->n_low;
        work->pol->Ared->x[Ared_nnz++] = work->data->A->x[ptr];
      }
    }
  }

  // Update the last element in Ared->p
  work->pol->Ared->p[work->data->n] = Ared_nnz;
}

/**
 * Form reduced right-hand side rhs_red = vstack[-q, l_low, u_upp]
 * @param  work Workspace
 * @param  rhs  right-hand-side
 * @return      reduced rhs
 */
static void form_rhs_red(OSQPWorkspace *work, c_float *rhs) {
  c_int j;

  // Form the rhs of the reduced KKT linear system
  for (j = 0; j < work->data->n; j++) { // -q
    rhs[j] = -work->data->q[j];
  }

  for (j = 0; j < work->pol->n_low; j++) { // l_low
    rhs[work->data->n + j] = work->data->l[work->pol->Alow_to_A[j]];
  }

  for (j = 0; j < work->pol->n_upp; j++) { // u_upp
    rhs[work->data->n + work->pol->n_low + j] =
      work->data->u[work->pol->Aupp_to_A[j]];
  }
}

/**
 * Compute dual variable y from yred
 * @param work Workspace
 * @param yred Dual variables associated to active constraints
 */
static void get_ypol_from_yred(OSQPWorkspace *work, c_float *yred) {
  c_int j;

  // If there are no active constraints
  if (work->pol->n_low + work->pol->n_upp == 0) {
    vec_set_scalar(work->pol->y, 0., work->data->m);
    return;
  }

  // NB: yred = vstack[ylow, yupp]
  for (j = 0; j < work->data->m; j++) {
    if (work->pol->A_to_Alow[j] != -1) {
      // lower-active
      work->pol->y[j] = yred[work->pol->A_to_Alow[j]];
    } else if (work->pol->A_to_Aupp[j] != -1) {
      // upper-active
      work->pol->y[j] = yred[work->pol->A_to_Aupp[j] + work->pol->n_low];
    } else {
      // inactive
      work->pol->y[j] = 0.0;
    }
  }
}

c_int update_active(OSQPWorkspace *work, OSQPRepolish *data) {
  c_int mred, exitflag;

  if(data->mred >= 0) cleanup(work, data);

  // Form Ared by assuming the active constraints and store in work->pol->Ared
  mred = form_Ared(work, data);

  if (mred < 0) { // work->pol->red = OSQP_NULL
    // Polishing failed
    work->info->status_polish = -1;
    data->mred = -1;

    return -1;
  }
  data->mred = mred;

  fill_Ared(work, data);

  // Form and factorize reduced KKT
  exitflag = init_linsys_solver(&data->plsh, work->data->P, work->pol->Ared,
                                work->settings->delta, OSQP_NULL,
                                work->settings->linsys_solver, 1);

  if (exitflag) {
    // Polishing failed
    work->info->status_polish = -1;

    // Memory clean-up
    if (work->pol->Ared) csc_spfree(work->pol->Ared);
    work->pol->Ared = NULL;
    if (data->rhs) c_free(data->rhs);
    data->rhs = NULL;
    data->mred = -1;

    return 1;
  }

  // Form reduced right-hand side rhs_red
  data->rhs_red = c_malloc(sizeof(c_float) * (work->data->n + mred));

  if (!data->rhs_red) {
    // Polishing failed
    work->info->status_polish = -1;

    // Memory clean-up
    if (work->pol->Ared) csc_spfree(work->pol->Ared);
    work->pol->Ared = NULL;
    if (data->rhs) c_free(data->rhs);
    data->rhs = NULL;
    if (data->plsh) data->plsh->free(data->plsh);
    data->plsh = NULL;
    data->mred = -1;

    return -1;
  }

  // Allocate polished solution
  data->pol_sol = c_malloc(sizeof(c_float) * (work->data->n + mred));

  if (!data->pol_sol) {
    // Polishing failed
    work->info->status_polish = -1;

    // Memory clean-up
    if (work->pol->Ared) csc_spfree(work->pol->Ared);
    work->pol->Ared = NULL;
    if (data->rhs) c_free(data->rhs);
    data->rhs = NULL;
    if (data->plsh) data->plsh->free(data->plsh);
    data->plsh = NULL;
    if (data->rhs_red) c_free(data->rhs_red);
    data->rhs_red = NULL;
    data->mred = -1;

    return -1;
  }

  return 0;
}

c_int repolish(OSQPWorkspace *work, OSQPRepolish *data) {
  c_int exitflag, polish_successful;

  #ifdef PROFILING
    osqp_tic(work->timer); // Start timer
  #endif /* ifdef PROFILING */

  if(data->mred < 0) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  fill_Ared(work, data);
  form_rhs_red(work, data->rhs_red);
  prea_vec_copy(data->rhs_red, data->pol_sol, work->data->n + data->mred);
  
  // Solve the reduced KKT system
  data->plsh->solve(data->plsh, data->pol_sol);

  // Perform iterative refinement to compensate for the regularization error
  exitflag = iterative_refinement(work, data, data->pol_sol, data->rhs_red);

  if (exitflag) {
    // Polishing failed
    work->info->status_polish = -1;  
    return -1;
  }

  // Store the polished solution (x,z,y)
  prea_vec_copy(data->pol_sol, work->pol->x, work->data->n);   // pol->x
  mat_vec(work->data->A, work->pol->x, work->pol->z, 0); // pol->z
  get_ypol_from_yred(work, data->pol_sol + work->data->n);     // pol->y

  // Ensure (z,y) satisfies normal cone constraint
  project_normalcone(work, work->pol->z, work->pol->y);

  // Compute primal and dual residuals at the polished solution
  update_info(work, 0, 1, 1);

    // Check if polish was successful
  polish_successful = (work->pol->pri_res < work->info->pri_res &&
                       work->pol->dua_res < work->info->dua_res) || // Residuals
                                                                    // are
                                                                    // reduced
                      (work->pol->pri_res < work->info->pri_res &&
                       work->info->dua_res < 1e-10) ||              // Dual
                                                                    // residual
                                                                    // already
                                                                    // tiny
                      (work->pol->dua_res < work->info->dua_res &&
                       work->info->pri_res < 1e-10);                // Primal
                                                                    // residual
                                                                    // already
                                                                    // tiny

  // Update solver information
  work->info->status_polish = polish_successful ? 1 : -1;
  work->info->obj_val       = work->pol->obj_val;
  work->info->pri_res       = work->pol->pri_res;
  work->info->dua_res       = work->pol->dua_res;

  // Update (x, z, y) in ADMM iterations
  // NB: z needed for warm starting
  prea_vec_copy(work->pol->x, work->x, work->data->n);
  prea_vec_copy(work->pol->z, work->z, work->data->m);
  prea_vec_copy(work->pol->y, work->y, work->data->m);

  // Store solution
  store_solution(work);

  #ifdef PROFILING
    work->info->polish_time = osqp_toc(work->timer);
  #endif /* ifdef PROFILING */
}

void cleanup(OSQPWorkspace *work, OSQPRepolish *data) {
  if (work->pol->Ared) csc_spfree(work->pol->Ared);
  work->pol->Ared = NULL;
  clean(data);
}

void clean(OSQPRepolish *data) {
  if (data->rhs) c_free(data->rhs);
  data->rhs = NULL;
  if (data->plsh) data->plsh->free(data->plsh);
  data->plsh = NULL;
  if (data->rhs_red) c_free(data->rhs_red);
  data->rhs_red = NULL;
  if (data->pol_sol) c_free(data->pol_sol);
  data->pol_sol = NULL;
  data->mred = -1;
}

void init_repolish(OSQPRepolish *data) {
  data->rhs = NULL;
  data->plsh = NULL;
  data->rhs_red = NULL;
  data->pol_sol = NULL;
  data->mred = -1;
}