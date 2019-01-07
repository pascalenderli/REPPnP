/**
 * @file reppnp.cpp
 * @author Pascal Enderli
 *
 * @Reference Paper: [1] Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection [Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer]
 * @date 10.10.2018
 */

#include <algorithm>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct> //unsupported by official eigen
#include <iostream>
#include <float.h>

#include "reppnp.h"

/// @brief Namespace for public function.
namespace REPPnP
{
  using namespace internal;


/**
 * @details
 * @n This function is based on the following paper:
 * @n Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection [Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer]
 * @n REPPnP rejects outliers by nullspace estimation rather than RANSAC. This is significantly faster than using Ransac with similar accuracy.
 * @n
 * @n ### Nomenclature Convention for the Transformations.
 * @n #### Upper Case Letters:
 * @n P represents Datapoints in World Coordinate Frame.
 * @n T represent Translations between Coordinate Frames.
 * @n R represent Rotations between Coordinate Frames.
 * @n K Camera calibration Matrix. Tacitly assumed as K_ic.
 * @n 
 * @n #### Lower Case Letters represent coordinate frames:
 * @n c : Camera Frame (3D)
 * @n w : World Frame (3D)
 * @n i : Image Frame in pixel Coordinates (2D)
 * @n n : Normalized coordinates (3D)
 * @n
 * @n #### Exemples
 * @n T_cw : Translation from w to c Frame.
 * @n c_P : 3D Points represented in Camera Frame.
 * @n
 * @n #### Mathematical Background
 * @n i_P = K_ic * c_P  with c_P = [R|t] * w_P
 * @n i_P = K_ic *[R_cw|t_cw] * w_P
 * @n if i_P, w_P correspondences and K_ic are given, this function finds R_cw and t_cw in a robust and efficient way.
 * @n For more details read the paper mentioned above.
 * @n
 * @param i_P[in] 2D pixel coordinates in image coordinate frame. Point features stacked in a 2xn Matrix.
 * @param w_P[in] Spacial 3D coordinates in the World coordinate Frame. Point features stacked in a 3xn Matrix.
 * @param K[in] the Camera intrinsic Calibration Matrix.
 * @param opt[in] Explained here: REPPnP_Options
 * @return Tuple of two values. The estimated Transform and a vector containing the indeces of the inliers. @n The Transform represents the extrinsics of the camera (T_cw; R_cw).
 */
  std::tuple<TF, std::vector<idx_t>> REPPnP(Coordinates2D_t i_P, Coordinates3D_t w_P, const CameraCalibration_t K, const REPPnP_Options& opt)
  {
    // Initialization of variables
    ControlPointBase c_w = opt.control_point_base_w;
    std::vector<idx_t> inliers_idx;
    nullspace_t Estimated_Nullspace;
    TF T_cw;

    try {
      if (i_P.cols() != w_P.cols()) {
        throw REPPnPException("Number of 3D points and 2D points as imput parameter have to be equal.");
      }
      if (i_P.cols() < 6) {
        throw REPPnPException("Not enough correspondences provided. Min number of correspondences is 6.");
      }



      // Transform Pixel coordinates to coordinates in the camera coordinate Frame.
      // REPPnP function input is i_P and w_P
      // For finding R_cw and t_cw, c_P and w_P is needed. => K * c_P = i_P => solve for c_P
      Coordinates2D_t c_P = Image2CameraCoordinates(i_P, K);

      // Normalize 3D Point Data
      ScaledTF T_wn = Normalize(w_P);

      // Apply REPPnP to find R_cw and t_cw
      M_t M = ComputeM(c_P, w_P, c_w);

      robust_nullspace_estimation(M, opt.dim_kernel, opt.max_algebraic_error, Estimated_Nullspace, inliers_idx);
      if (inliers_idx.size() < 6) {
        throw REPPnPException("Not enough inliers left to determine a solution.  Min number of inliers is 6.");
      }

      TF T_cn = kernelPnP(c_w, Estimated_Nullspace, opt.dim_kernel, opt.refine);

      // Denormalize
      // The following three lines are correct for normalization scale == 1;
      //Translation_t t_cw = T_cn.GetRotmatrix() * T_nw.GetTranslation() + T_cn.GetTranslation();
      //RotMat_t r_cw = T_cn.GetRotmatrix();
      //TF T_cw = TF(r_cw, t_cw);

      // Denormalize with normalization scale != 1 could be numerically unstable if the scale is large!
      TF T_wc = T_wn * T_cn.Inv();
      T_cw = T_wc.Inv();
    }

    catch(const REPPnPException& e)
    {
      throw e;
    }

    catch(const std::exception& e)
    {
      throw REPPnP::REPPnPException("A std exception was caught: " + std::string(e.what()) );
    }

    catch(...)
    {
      throw REPPnP::REPPnPException("Unknown exception.");
    }

    return make_tuple(T_cw, inliers_idx);
  }

  REPPnPException::REPPnPException(const std::string& message)
  : message_(message)
  {}

  const char* REPPnPException::what() const throw()
  {
    return this->message_.c_str();
  }

  /// @brief Internal functions used for REPPnP implementation. Changing those can cause unexpected behaviour.
  namespace internal
  {


    /**
     * @brief Compute the internal Matrix M
     *
     * @param i_P[in] 2D feature point corespondences.
     * @param w_P[in] 3D point corespondences.
     * @param c_w[in] The basis to compute the controlpoints alpha.
     * @return Matrix M which is a kronecker product of the parameters alpha and a pattern containing the 2D features.
     */
    M_t ComputeM(Coordinates2D_t i_P, const Coordinates3D_t w_P, const ControlPointBase c_w)
    {
      size_t n_matches = i_P.cols();
      alphas_t alphas = ComputeAlphas(w_P, c_w);

      Eigen::Map<Eigen::VectorXd> i_P_vectorized(i_P.data(), i_P.size());
      Eigen::Matrix<double, 2, 3> pattern;
      pattern << 1, 0, -1, 0, 1, -1;

      M_t M;
      M.resize(n_matches, Eigen::NoChange);
      M = Eigen::kroneckerProduct<alphas_t, Eigen::Matrix<double, 2, 3>>(alphas, pattern);

      for(int i=2; i<12; i+=3)
      {
        M.col(i) =  M.col(i).cwiseProduct(i_P_vectorized);
      }
      return M;
    }


    /**
     * @brief Compute the Controlpoints alpha of the 3D Points
     * @param w_P Coordinates of Point 3D corespondences.
     * @param c_w The basis to compute the controlpoints alpha.
     * @return 3D Point coordinates in Controlpoint base.
     */
    alphas_t ComputeAlphas(const Coordinates3D_t& w_P, const ControlPointBase& c_w)
    {
      size_t n_matches = w_P.cols();

      Eigen::Matrix<double, 4, 4> c_w_extended;
      c_w_extended << c_w, 1, 1, 1, 1;

      Eigen::Matrix<double, 4, Eigen::Dynamic> w_P_homogenious;
      Eigen::Matrix<double, 1, Eigen::Dynamic> ones_vector;

      ones_vector = Eigen::MatrixXd::Ones(1, n_matches);
      w_P_homogenious.resize(Eigen::NoChange, n_matches);
      w_P_homogenious << w_P, ones_vector;

      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> alphas;
      alphas.resize(Eigen::NoChange, n_matches);
      alphas = c_w_extended.inverse() * w_P_homogenious;

      alphas.transposeInPlace();
      return alphas;
    }


    /**
     * @brief Nullspace etimation.
     * @param M Internal[in] matrix.
     * @param dim_kernel[in] Dimension for the estimated kernel.
     * @param max_algebraic_error[in] Threthold of the maximal accepted algebraic error such that a point is taken as inlier.
     * @param Nullspace_estimation_out[out] Return of the estimated Nullspace of M.
     * @param inliers_idx[out] Vector filled with the indeces of the inliers.
     */
    void robust_nullspace_estimation(M_t M, short dim_kernel, double max_algebraic_error, nullspace_t & Nullspace_estimation_out, std::vector<idx_t>& inliers_idx)
    {
      size_t n_rows_M = M.rows();
      size_t n_measurements = n_rows_M/2;

      std::vector<idx_t> idx_M(n_rows_M);
      unsigned int it = 0;
      generate(idx_M.begin(), idx_M.end(), [&it]()mutable{unsigned int tmp=it; ++it; return tmp; });

      M_t M_of_inliers;
      M_of_inliers.resize(n_rows_M, Eigen::NoChange);

      Eigen::Matrix<double, 12, 12> MTM;
      Eigen::Matrix<double, 12, 12> V;
      Eigen::Matrix<double, 12, 12> result_V;

      std::vector<double> algebraic_errors;
      std::vector<double> algebraic_errors_sorted;
      std::vector<idx_t> algebraic_errors_idx;
      std::vector<idx_t> result_algebraic_errors_idx;
      double Q25;
      unsigned int n_inliers = n_measurements;
      double previous_Q25 = DBL_MAX;

      for(int i=0; i<20; ++i)
      {
        M_of_inliers = SelectRows(M, idx_M);
        MTM = M_of_inliers.transpose() * M_of_inliers;

        Eigen::BDCSVD<Eigen::Matrix<double, 12, 12>> svd(MTM, Eigen::ComputeFullV );
        V = svd.matrixV();
        algebraic_errors = GetAlgebraicErrors(M, V);

        algebraic_errors_sorted = algebraic_errors;
        std::sort(algebraic_errors_sorted.begin(), algebraic_errors_sorted.end());
        algebraic_errors_idx = argsort(algebraic_errors);

        //algebraic error of the correspondence that is at the boundary of the lowest 25% quartile.
        Q25 = algebraic_errors_sorted[(int)(n_measurements / 4 - 1)];
        //pprint::pprint(Q25,"Q25 Error");
        n_inliers = std::lower_bound(algebraic_errors_sorted.begin(), algebraic_errors_sorted.end(), std::max(Q25, max_algebraic_error)) - algebraic_errors_sorted.begin();
        //pprint::pprint(algebraic_errors_sorted,"errors sorted; size"+to_string(algebraic_errors_sorted.size())+"; n_inliers="+to_string(n_inliers) );
        if(Q25 >=  previous_Q25)
        {
          break;
        }
        else
        {
          previous_Q25 = Q25;
          result_V = V;
          result_algebraic_errors_idx = algebraic_errors_idx;
          result_algebraic_errors_idx.resize(n_inliers);
        }
        idx_M = UpdateInlierIdxM(algebraic_errors_idx, n_inliers);

      }
      std::vector<idx_t> idx_kernel(dim_kernel);
      unsigned int kernel_startid = 12-dim_kernel;
      generate(idx_kernel.begin(), idx_kernel.end(), [&kernel_startid]()mutable{ unsigned int tmp=kernel_startid; ++kernel_startid; return tmp; });
      Nullspace_estimation_out = SelectCols(result_V, idx_kernel);
      inliers_idx = result_algebraic_errors_idx;
    }


    /**
     * @details Form a new Matrix by stacking a selection of rows of a given matrix.
     * @param A The Full Matrix.
     * @param idx Vector of Row indeces which will be selected.
     * @return Stacked matrix of the selected rows.
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SelectRows(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, std::vector<idx_t> idx)
    {
      int n_selected_rows = idx.size();
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sub_matrix;
      sub_matrix.resize(n_selected_rows, A.cols());

      for(int i=0; i<n_selected_rows; ++i)
      {
        sub_matrix.row(i) = A.row(idx.at(i));
      }
      return sub_matrix;
    }


    /**
     * @details Form a new Matrix by stacking a selection of columns of a given matrix.
     * @param A The Full Matrix.
     * @param idx vector of column indeces which will be selected.
     * @return Stacked matrix of the selected Columns.
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SelectCols(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, std::vector<idx_t> idx)
    {
      int n_selected_cols = idx.size();
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sub_matrix;
      sub_matrix.resize(A.rows(), n_selected_cols);

      for(int i=0; i<n_selected_cols; ++i)
      {
        sub_matrix.col(i) = A.col(idx.at(i));
      }
      return sub_matrix;
    }


    /**
     * @param M internal Matrix.
     * @param V Eigenvectors of M.
     * @return Vector containing the algebraic error of every point2D/point3D pair.
     */
    std::vector<double> GetAlgebraicErrors(M_t M, Eigen::Matrix<double, 12, 12> V)
    {
      if (M.rows()%2)
      {
        throw std::invalid_argument("M must have a even number of rows!");
      }

      const unsigned int size = M.rows()/2;

      std::vector<idx_t> odd_idx(size);
      std::vector<idx_t> even_idx(size);
      unsigned int odd_id = 1;
      unsigned int even_id = 0;
      std::generate(odd_idx.begin(), odd_idx.end(), [&odd_id]()mutable{ unsigned int tmp=odd_id; odd_id += 2; return tmp; });
      std::generate(even_idx.begin(), even_idx.end(), [&even_id]()mutable{ unsigned int tmp=even_id; even_id += 2; return tmp; });

      M_t odd_M; odd_M.resize(size, Eigen::NoChange);
      M_t even_M; even_M.resize(size, Eigen::NoChange);
      odd_M = SelectRows(M, odd_idx);
      even_M = SelectRows(M, even_idx);

      Eigen::Matrix<double, Eigen::Dynamic, 1> error21; error21.resize(size, Eigen::NoChange);
      Eigen::Matrix<double, Eigen::Dynamic, 1> error22; error22.resize(size, Eigen::NoChange);
      Eigen::Matrix<double, Eigen::Dynamic, 1> error2; error2.resize(size, Eigen::NoChange);

      error21 = odd_M * V.col(11);
      error22 = even_M * V.col(11);
      error2 = error21.cwiseProduct(error21)  + error22.cwiseProduct(error22);
      error2 = error2.cwiseSqrt();

      std::vector<double> vec_error2(error2.data(), error2.data() + error2.size());

      return vec_error2;
    }


    /**
     * @param v[in] Eigen vector which has to be sorted.
     * @return The indices that would sort an array in a vector.
     */
    std::vector<idx_t> argsort(std::vector<double> v)
    {
      std::vector<idx_t> idxs(v.size());
      unsigned int it = 0;
      std::generate(idxs.begin(), idxs.end(), [&it]()mutable{unsigned int tmp=it; ++it; return tmp; });

      std::sort( idxs.begin(), idxs.end(), [&v](int a, int b) {return v[a] < v[b];} );
      return idxs;
    }


    /**
     *
     * @param algebraic_errors_idx Indeces which would sort the algebraic error in ascending order.
     * @param n_inliers Nr of inliers detected in this iteration.
     * @return Indeces of Matrix m which are inliers.
     */
    std::vector<idx_t> UpdateInlierIdxM(std::vector<idx_t>  algebraic_errors_idx, unsigned int n_inliers)
    {
      std::vector<idx_t> result(2 * n_inliers);
      for(unsigned int i=0; i<n_inliers; i++)
      {
        result[2*i] = 2 * algebraic_errors_idx[i];
        result[2*i+1] = 2 * algebraic_errors_idx[i] + 1;
      }
      return result;
    }


    // ControlPoint Base.
    ControlPoint::ControlPoint(){}

    ControlPoint::ControlPoint(ControlPointBase base)
    { this->basis_ = base; }

    ControlPoint::ControlPoint(const ControlPoint&) = default;
    ControlPoint::ControlPoint(ControlPoint&&) = default;

    ControlPoint::~ControlPoint()
    {}

    ControlPoint& ControlPoint::operator=(const ControlPoint&)=default;
    ControlPoint& ControlPoint::operator=(ControlPoint&&)=default;

    ControlPointBase ControlPoint::GetControlPointsBase() const
    {
      return this->basis_;
    }

    Eigen::Matrix<double, 3, 1> ControlPoint::GetMeanControlPoint() const
    {
      return this->basis_.rowwise().mean();
    }

    ControlPointBase  ControlPoint::GetCenteredControlPoints() const
    {
      return this->basis_.colwise() - GetMeanControlPoint();
    }

    ControlPointBase ControlPoint::GetNormalizedCenteredControlPoints() const
    {
      return GetCenteredControlPoints() / GetNorm();
    }

    float ControlPoint::GetNorm() const
    {
      return GetCenteredControlPoints().norm();
    }

    ControlPoint ControlPoint::ProjectToNullspace(nullspace_t nullspace_base) const
    {
      Eigen::Map<const Eigen::VectorXd> this_base_vectorized(this->basis_.data(), 12);

      // Project previous solution into nullspace. Where:  Kernel * a = previous solution.
      Eigen::Matrix<double, 4, 1> a;
      a = nullspace_base.colPivHouseholderQr().solve(this_base_vectorized);

      Eigen::Matrix<double, 12, 1> a_in_nullspace = nullspace_base * a;
      Eigen::Map<ControlPointBase> this_base_in_nullspace(a_in_nullspace.data(), 3, nullspace_base.cols());
      return ControlPoint(this_base_in_nullspace);
    }


    /**
     * @param c_w The basis to compute the controlpoints alpha.
     * @param Estimated_Nullspace Nullspace of Matrix M.
     * @param dim_kernel Dimension of the estimated kernel.
     * @param refine if the xode should optimize the result.
     * @return Transform aka: Rotation and translation of the registered image.
     */
    TF kernelPnP(ControlPointBase c_w, nullspace_t Kernel, short dim_kernel, bool refine)
    {

      // Controlpoint Basis from estimated Kernel.
      Eigen::Matrix<double, Eigen::Dynamic, 1> Kernel_last_col = Kernel.col(dim_kernel-1);
      Eigen::Map<ControlPointBase> v_k(Kernel_last_col.data(), 3, 4); //Reshape
      if(v_k.row(2).mean() < 0)
      { v_k = -v_k; }

      // chosen Control Point Base
      ControlPoint x(c_w);

      // base of Controlpoints of the 3D Points in camera Coordinate Frame.
      ControlPoint y(v_k);

      // Find Transform between base x and y
      RotMat_t r;
      double scaler;
      Eigen::Matrix<double, 3, 1> t;
      ProcrustesSolver(x, y, r, scaler, t);

      // refinement
      ControlPoint previous_solution;
      ControlPoint previous_solution_in_nullspace_base;
      if(refine)
      {
        double current_error;
        double previous_error = DBL_MAX;
        ControlPointBase error_matrix;

        for(int i=0; i<20; ++i)
        {
          previous_solution = ControlPoint(r * (x.GetControlPointsBase() - t.replicate(1, 4)));
          previous_solution_in_nullspace_base = ControlPoint(previous_solution.ProjectToNullspace(Kernel));

          error_matrix = r.transpose() * previous_solution_in_nullspace_base.GetControlPointsBase() + t.replicate(1, 4) - x.GetControlPointsBase();
          current_error = error_matrix.norm();
          if(current_error > previous_error && i>1)
          {
            break;
          }
          else
          {
            ProcrustesSolver(x, previous_solution_in_nullspace_base, r, scaler, t);
            previous_error = current_error;
          }
        }
      }
      Coordinate3D_t translation = (scaler * previous_solution_in_nullspace_base.GetMeanControlPoint()) - (r * x.GetMeanControlPoint());
      return TF(r, translation);
    }


    /**
     * @details The Procrustes problem determines the closest orthogonal transform by minimizing the Frobenius norm.
     * @param x[in] Base of the Controlpoints
     * @param y[in] Base of control points of 3D Points in the camera coordinate system.
     * @param r_out[out] Rotation Matrix.
     * @param b_out[out] Scaler Component.
     * @param mc_out[out] Translation Component
     */
    void ProcrustesSolver(ControlPoint x, ControlPoint y, RotMat_t& r_out, double& scaler_out, Eigen::Matrix<double, 3, 1>& t_out)
    {
      // Cross correlation between a and y
      Eigen::Matrix3d correlation = x.GetNormalizedCenteredControlPoints() * y.GetNormalizedCenteredControlPoints().transpose();

      // svd
      Eigen::JacobiSVD<Eigen::Matrix3d> svd(correlation, Eigen::ComputeFullU | Eigen::ComputeFullV);

      // find R
      Eigen::Matrix3d temp = svd.matrixV()*svd.matrixU().transpose();
      double sigma_z = sign(temp.determinant());
      Eigen::DiagonalMatrix<double, 3> S_artificial(1, 1, sigma_z);

      r_out = svd.matrixV() * S_artificial * svd.matrixU().transpose();

      // find scaler
      scaler_out = svd.singularValues().sum() * x.GetNorm()/y.GetNorm();

      // find translation
      t_out = x.GetMeanControlPoint() - scaler_out * r_out.transpose() * y.GetMeanControlPoint();
    }

    int sign(double val)
    { return (0 < val) - (val < 0); }

    /**
     * @details
     * @n In order PnP works numerically properly the 3D datapoints need to be normalized.
     * @n The centroid of the 3D points are set to the origin and the average distance of all the points are set to 1.
     * @n !!! If the scale is big (some 100) the Denormalization can become numerically unstable, because of error upscaling.
     * @param w_P[in,out] 3xn stacked list of 3D datapoints will be normalized after execution.
     * @param use_scaling If Data has too much scale involved you may try to not scale the data in normalization step. (Defult is set to also scale the data. Use_scaling = true).
     * @return Transform which contains the applied scale and translation for normalization of the datapoints w_P[in] into normalized coordinate frame w_P[out] (which is then actually n_P) \n and the transform parameters s_wn and t_wn.
     */
    ScaledTF Normalize(Coordinates3D_t& w_P, bool use_scaling)
    {
      // Initialize Transforms
      double s_wn;
      Translation_t t_wn;

      // Find mean distance from 3D datapoints to the origin.
      t_wn = w_P.rowwise().mean();

      // Center the Centroid of the datapoints around the origin. (inverse transform for t_wn)
      Eigen::Matrix<double, 3, Eigen::Dynamic> w_P_centered = w_P.colwise() - t_wn;

      // Find the mean distance of all the datapoints to their centroid.
      Eigen::Matrix<double, 1, Eigen::Dynamic> norm = w_P_centered.colwise().norm();
      double mean_distance_to_centroid = norm.mean();

      // Get the scaler which would transform the data to have its mean distance to one.
      s_wn = mean_distance_to_centroid;

      if(use_scaling == false)
      {
        // Do not scale
        s_wn = 1;
      }

      // Apply inverse normalization transform to data.
      w_P = 1/s_wn * w_P_centered;

      return ScaledTF(t_wn, s_wn);
    }

    /**
     * @details
     * @n Projection of 3D coordinates P to image pixel coordinates i_P Projection rule:
     * @n i_P = K * c_P
     * @n where c_P = [R|t] * P
     * @n REPPnP function input is i_P and P.
     * @n But for finding R and t, c_P and P is needed. => K*c_P=i_P => solve for c_P
     * @param i_P[in] 2D feature points in image pixel coordinate Frame.
     * @param K[in] Camera Calibration Matix [intrinsics]
     * @return 2D Feature coordinates in camera coordinate frame.
     */
    Coordinates2D_t Image2CameraCoordinates(Coordinates2D_t i_P, CameraCalibration_t K)
    {
      Coordinates3D_t p_homogeneous = K.fullPivLu().solve(i_P.colwise().homogeneous());
      return p_homogeneous.colwise().hnormalized();
    }

  } //namespace internal
} // namespace REPPnP
