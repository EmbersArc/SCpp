#include <fstream>
#include <filesystem>

#include <fmt/format.h>

#include "activeModel.hpp"
#include "ecosWrapper.hpp"
#include "discretization.hpp"
#include "successiveConvexificationProblem.hpp"
#include "parameterServer.hpp"

class freeFinalTimeAlgorithm
{
  public:
    explicit freeFinalTimeAlgorithm(std::shared_ptr<Model> model);
    void initialize();
    void solve(bool warm_start = false);
    void getSolution(Model::dynamic_matrix_t &X, Model::dynamic_matrix_t &U, double &t);

  private:
    void cacheIndices();
    void readSolution();
    void loadParameters();
    bool iterate();

    ParameterServer param;

    size_t K;

    std::shared_ptr<Model> model;

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