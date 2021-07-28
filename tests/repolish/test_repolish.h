#include "osqp.h"    // OSQP API
#include "repolish.h"
#include "minunit.h" // Basic testing script header

static const char* test_repolish_solve()
{
  // Load problem data
  c_float P_x[3] = { 4.0, 1.0, 2.0, };
  c_int   P_nnz  = 3;
  c_int   P_i[3] = { 0, 0, 1, };
  c_int   P_p[3] = { 0, 1, 3, };
  c_float q[2]   = { 1.0, 1.0, };
  c_float A_x[4] = { 1.0, 1.0, 1.0, 1.0, };
  c_int   A_nnz  = 4;
  c_int   A_i[4] = { 0, 1, 0, 2, };
  c_int   A_p[3] = { 0, 2, 4, };
  c_float l[3]   = { 1.0, 0.0, 0.0, };
  c_float u[3]   = { 1.0, 0.7, 0.7, };
  c_int n = 2;
  c_int m = 3;

  c_int exitflag;

  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // OSQPData

  // Populate data
  data = (OSQPData *)c_malloc(sizeof(OSQPData));
  data->n = n;
  data->m = m;
  data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
  data->q = q;
  data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
  data->l = l;
  data->u = u;

  // Define solver settings as default
  osqp_set_default_settings(settings);

  // Setup workspace
  exitflag = osqp_setup(&work, data, settings);

  // Setup correct
  mu_assert("Repolish test solve: Setup error!", exitflag == 0);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Repolish test solve: Error in solver status!",
	work->info->status_val == OSQP_SOLVED);

  OSQPRepolish *repolish_data = (OSQPRepolish *)c_malloc(sizeof(OSQPRepolish));
  init_repolish(repolish_data);

  mu_assert("Repolish test solve: Error in repolish data!",
	repolish_data->mred == -1);

  update_active(work, repolish_data);

  mu_assert("Repolish test solve: Error in active set update!",
	repolish_data->mred >= 0);

  repolish(work, repolish_data);

  mu_assert("Repolish test solve: Error in repolish!",
	work->info->status_polish == 1);

  cleanup(work, repolish_data);

  mu_assert("Repolish test solve: Error in cleanup!",
	repolish_data->mred == -1);

  // Clean workspace
  osqp_cleanup(work);
  c_free(data->A);
  c_free(data->P);
  c_free(data);
  c_free(settings);

  return 0;
}


static const char* test_repolish()
{
  mu_run_test(test_repolish_solve);

  return 0;
}
