/* Solution polish based on assuming the active set */
#ifndef REPOLISH_H
# define REPOLISH_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


# include "types.h"

# ifndef EMBEDDED

/**
 * Polish structure
 */
typedef struct {
  c_float *rhs_red;     ///< Reduced right hand side
  LinSysSolver *plsh;   ///< Linear solver
  c_float *pol_sol;     ///< Polished solution
  c_float *rhs;         ///< Right hand side for iterative refinement
  c_int mred;           ///< Number of reduced rows
} OSQPRepolish;
# endif // ifndef EMBEDDED

/**
 * @brief Update the set of active constrints
 * 
 * @param  work Workspace
 * @param  data Repolish data
 * @return      Exitflag
 */
c_int update_active(OSQPWorkspace *work, OSQPRepolish *data);

/**
 * Solution repolish: Solve equality constrained QP with assumed active
 *constraints stored in data
 * @param  work Workspace
 * @param  data Repolish data
 * @return      Exitflag
 */
c_int repolish(OSQPWorkspace *work, OSQPRepolish *data);

/**
 * @brief Cleanup repolish data and free memory
 * 
 * @param  work Workspace
 * @param  data Repolish data
 */
void cleanup(OSQPWorkspace *work, OSQPRepolish *data);
void clean(OSQPRepolish *data);

/**
 * @brief Initialises the repolish data
 * 
 * @param  data Repolish data
 */
void init_repolish(OSQPRepolish *data);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef REPOLISH_H
