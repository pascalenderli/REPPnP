/**
 * @file reppnp.h
 * @author Pascal Enderli
 * @brief Robust Efficient Procrustes Perspective-n-Point (REPPnP)
 * @details
 * @n
 * @n This is c++ program of REPPnP [1].
 * @n
 * @n
 * @n Quote of the abstract from the referenced paper by L.Ferraz et al. [1]:
 * @n
 * @n "We propose a real-time, robust to outliers and accurate
 * @n solution to the Perspective-n-Point (PnP) problem. The
 * @n main advantages of our solution are twofold: first, it in-
 * @n tegrates the outlier rejection within the pose estimation
 * @n pipeline with a negligible computational overhead; and sec-
 * @n ond, its scalability to arbitrarily large number of correspon-
 * @n dences. Given a set of 3D-to-2D matches, we formulate
 * @n pose estimation problem as a low-rank homogeneous sys-
 * @n tem where the solution lies on its 1D null space. Outlier
 * @n correspondences are those rows of the linear system which
 * @n perturb the null space and are progressively detected by
 * @n projecting them on an iteratively estimated solution of the
 * @n null space. Since our outlier removal process is based on
 * @n an algebraic criterion which does not require computing the
 * @n full-pose and reprojecting back all 3D points on the image
 * @n plane at each step, we achieve speed gains of more than
 * @n 100× compared to RANSAC strategies. An extensive exper-
 * @n imental evaluation will show that our solution yields accu-
 * @n rate results in situations with up to 50% of outliers, and can
 * @n process more than 1000 correspondences in less than 5ms." (L.Ferraz et al. [1])
 * @n
 * @Reference Paper: [1] Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection [Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer]
 * @Reference Paper: [2] EPnP: An Accurate O(n) Solution to the PnP Problem [Vincent LepetitFrancesc, Moreno-Noguer, Pascal Fua]
 * @date 10.10.2018
 */

#ifndef ABSOLUTE_SLAM_REPPNP_H
#define ABSOLUTE_SLAM_REPPNP_H

#include <tuple>

#include "Eigen/Dense"
#include "scaled_transform.h"

namespace REPPnP
{
  //using namespace std;
  using idx_t = unsigned int;
  using Coordinates3D_t = Eigen::Matrix<double, 3, Eigen::Dynamic>;
  using Coordinates2D_t = Eigen::Matrix<double, 2, Eigen::Dynamic>;
  using ControlPointBase = Eigen::Matrix<double, 3, 4>;
  using CameraCalibration_t = Eigen::Matrix<double, 3, 3>;

  /// @brief Options for REPPnP function.
  struct REPPnP_Options
  {
    //Fill struct with default vaues.
    REPPnP_Options()
    {
      dim_kernel = 4;
      refine = true;
      control_point_base_w << 1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0, 1, 0;

      max_algebraic_error = 0.0175;
    }
    /** @brief Dimension of the kernel.
     *  @details Default Value = 4.
     */
    short dim_kernel;

    /**
     * @brief refine initial solution
     * @details Default Value = true.
     */
    bool refine;

    /** @brief threshold for max algebraic_error.
     *  @details
     *  @n Max accepted algebraic error for solving Mx=0 such that a ddatapoint is taken as inlier.
     *  @n This parameter is proportional to image_noise/f
     *  @n Default Value = 0.0175
     */
    float max_algebraic_error;

    /** @brief Control Point Basis representing the world coordinate Frame. (EPnP [Lepetit et all.])
     *  @details
     *  @n Default Value 3x4 Matrix  =
     *  @n [1, 0, 0, 0]
     *  @n [0, 1, 0, 0]
     *  @n [0, 0, 1, 0]
     */
     ControlPointBase control_point_base_w;
  };

  class REPPnPException: public std::exception
  {
  public:
    explicit REPPnPException(const std::string& message);

    virtual const char* what() const throw();

  protected:
    std::string message_;
  };




  // The REPPnP main Function.
  /**
   * @brief Solve PnP Problem with Algebraic Outlier Rejection. (Main Function)
   */
  std::tuple<TF, std::vector<idx_t>> REPPnP(Coordinates2D_t i_P, Coordinates3D_t w_P, const CameraCalibration_t K, const REPPnP_Options& opt = REPPnP_Options());

  namespace internal
  {
    using alphas_t = Eigen::Matrix<double, Eigen::Dynamic, 4 >;
    using M_t = Eigen::Matrix<double, Eigen::Dynamic, 12>;
    using nullspace_t = Eigen::Matrix<double, 12, Eigen::Dynamic>;

    /**
     * @brief Control Points describe a special 3D space basis.
     * @details
     * @n Control Points describe a special 3D space basis. The dimension of this special basis vectors will be four.
     * @n The 3D Points are then expressed as a weighted sum of four virtual controlpoints (aka homogeneous barycentric coordinates). These Coordinates are by design independent of the coordinate frame.
     * @n Explained more in detail in paper [2] EPnP: An Accurate O(n) Solution to the PnP Problem [Vincent LepetitFrancesc, Moreno-Noguer, Pascal Fua]
     * @param[in] c_w Basevectors as matrix of the used control points alpha.
    */
    class ControlPoint
    {
    public:
      /// @brief Default constructor.
      ControlPoint();

      /// @brief Construct from given base matrix.
      ControlPoint(const ControlPoint&);

      /// @brief Move constructor.
      ControlPoint(ControlPoint&&);

      /// @brief Copy constructor.
      ControlPoint(ControlPointBase base);

      /// @brief Destructor.
      virtual ~ControlPoint();

      // Operators
      /// @brief Copy assignment.
      ControlPoint& operator=(const ControlPoint& p);

      /// @brief Move assignment.
      ControlPoint& operator=(ControlPoint&& p);

      // Getters
      /// @brief Get the basis as matrix.
      ControlPointBase GetControlPointsBase() const;

      /// @brief Mean of the three Basisvectors.
      Eigen::Matrix<double, 3, 1> GetMeanControlPoint() const;

      /// @brief Basis centered in origin.
      ControlPointBase GetNormalizedCenteredControlPoints() const;

      /// @brief Project the basis to a Nullspace.
      ControlPoint ProjectToNullspace(nullspace_t nullspace_base) const;

      /// @brief  Frobenius norm of the basis.
      float GetNorm() const;

      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
      Eigen::Matrix<double, 3, 4>  GetCenteredControlPoints() const;
      ControlPointBase basis_;
    };

    /// @brief  Find the special Matrix M.
    M_t ComputeM(Coordinates2D_t i_P, const Coordinates3D_t w_P, const ControlPointBase c_w);

    /// @brief Compute the homogeneous barycentric coordinates of the 3D points in the control point basis.
    alphas_t ComputeAlphas(const Coordinates3D_t& w_P, const ControlPointBase& c_w);

    /// @brief Find the nullspace of matrix M.
    void robust_nullspace_estimation(M_t M, short dim_kernel, double max_algebraic_error, nullspace_t & Nullspace_estimation_out, std::vector<idx_t>& inliers_idx);

    /// @brief  Select certain rows of a eigen matrix.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SelectRows(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, std::vector<idx_t> idx);

    /// @brief Select certain columns of a eigen matrix.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SelectCols(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, std::vector<idx_t> idx);

    /// @brief Calculate the Algebraic errors of every data pair 2D/3D.
    std::vector<double> GetAlgebraicErrors(M_t M, Eigen::Matrix<double, 12, 12> V);

    /// @brief Get the indeces of Matrix m which will be inliers according to the algebraic error.
    std::vector<idx_t> UpdateInlierIdxM(std::vector<idx_t>  algebraic_errors_idx, unsigned int n_inliers);

    /// @brief Calculate the estimatet Rotation and Translation of the registered image from the estimated Nullspace of M.
    TF kernelPnP(ControlPointBase c_w, nullspace_t Kernel, short dim_kernel, bool refine);

    /// @brief Find the indeces which would sort an array. (similar to numpy argsort)
    std::vector<idx_t> argsort(std::vector<double> v);

    /// @brief Solves the Orthogonal Procrustes problem in closed-form as a minimization problem.
    void ProcrustesSolver(ControlPoint x, ControlPoint y, RotMat_t& r_out, double& scaler_out, Eigen::Matrix<double, 3, 1>& t_out);

    /// @brief Get the signum of a number.
    int sign(double val);

    /// @brief Get Parametres which would center the points w_P in origin and scale them to a mean distance to origin of 1.
    ScaledTF Normalize(Coordinates3D_t& w_P, bool use_scaling = true);

    /// @brief Transform Pixel coordinates to coordinates in the camera coordinate Frame.
    Coordinates2D_t Image2CameraCoordinates(Coordinates2D_t i_P, CameraCalibration_t K);
  }
}


#endif //ABSOLUTE_SLAM_REPPNP_H
