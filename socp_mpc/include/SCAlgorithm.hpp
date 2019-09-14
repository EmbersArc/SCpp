#include "activeModel.hpp"
#include "ecosWrapper.hpp"
#include "discretization.hpp"
#include "SCProblem.hpp"
#include "parameterServer.hpp"

class SCAlgorithm
{
public:
  /**
     * @brief Construct a new free Final Time Algorithm solver.
     * 
     * @param model     The system model.
     */
  explicit SCAlgorithm(Model::ptr_t model);

  /**
     * @brief Initializes the algorithm. Has to be called before solving the problem.
     * 
     */
  void initialize();

  /**
     * @brief Solves the system.
     * 
     * @param warm_start    Whether to reuse the last computed trajectory.
     */
  void solve(bool warm_start = false);

  /**
     * @brief Get the solution variables object.
     * 
     * @param X     The state trajectory.
     * @param U     The input trajectory.
     * @param t     The final time.
     */
  void getSolution(Model::dynamic_matrix_t &X, Model::dynamic_matrix_t &U, double &t);

private:
  /**
     * @brief Saves solution indices for performance.
     * 
     */
  void cacheIndices();

  /**
     * @brief Reads the solution variables X, U and sigma.
     * 
     */
  void readSolution();

  /**
     * @brief Loads the parameters from the configuration file.
     * 
     */
  void loadParameters();

  /**
     * @brief Performs a Successive Convexification iteration.
     * 
     * @return true     If converged.
     * @return false    If not converged.
     */
  bool iterate();

  ParameterServer param;

  size_t K;

  Model::ptr_t model;

  bool free_final_time;

  bool nondimensionalize;
  double weight_trust_region_time;
  double weight_trust_region_trajectory;
  double weight_virtual_control;
  double trust_region_factor;
  double nu_tol;
  double delta_tol;
  size_t max_iterations;

  Model::state_matrix_v_t A_bar;
  Model::control_matrix_v_t B_bar;
  Model::control_matrix_v_t C_bar;
  Model::state_vector_v_t S_bar;
  Model::state_vector_v_t z_bar;

  double sigma;
  Eigen::MatrixXd X;
  Eigen::MatrixXd U;

  size_t sigma_index;
  Eigen::MatrixXi X_indices;
  Eigen::MatrixXi U_indices;

  op::SecondOrderConeProgram socp;

  std::unique_ptr<EcosWrapper> solver;
};