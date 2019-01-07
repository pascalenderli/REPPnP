/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file transform.h
 * @brief Function prototypes for transformations
 *
 * @details
 * ### I recommend strongly to use the following Nomenclature Convention for Transformations to avoid Confusion when working with Transforms.
 * @n #### Upper Case Letters Represent the name of the Variable:
 * @n e.g.
 * @n P represents Datapoints.
 * @n T represents a Translation between two Coordinate Frames.
 * @n R represents a Rotation between two Coordinate Frames.
 * @n
 * @n #### Lower Case Letters (indeces) represent coordinate frames:
 * @n A Transform occurs <b> always </b> in pairs. One describes the transform from Frame a to frame b, the other is the inverse, and describes the Transform from frame b to a.
 * @n e.g.
 * @n c : Camera Frame (3D)
 * @n w : World Frame (3D)
 * @n
 * @n #### Transformation Exemples
 * @n T_cw : Translation from w to c Frame.
 * @n c_P : 3D Points represented in Camera Frame.
 * @n c_P = TF_cw * w_P; : Point P in w frame coordinates is transformed to its representation in c coordinate frame
 * @n TF_ca = TF_cw * TF_wa; Transform from a to c frame (TF_ca) is a concatenation of TF_cw and TF_wa.
 *
 */


#ifndef ABSOLUTE_SLAM_TRANSFORM_INTERFACE_H
#define ABSOLUTE_SLAM_TRANSFORM_INTERFACE_H


#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <type_traits>

#include <iostream>

/// 4x4 [Homogenious Transformation Matrix](https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html")
using HomogeniousTF_t = Eigen::Matrix<double, 4, 4>;
/// 3x3 [Rotation Matrix](https://en.wikipedia.org/wiki/Rotation_matrix)
using RotMat_t        = Eigen::Matrix<double, 3, 3>;
/// 3x1 Vector [x,y,z]<SUP>T</SUP>
using Translation_t   = Eigen::Matrix<double, 3, 1>;
/// 4x1 Vector [w,x,y,z]<SUP>T</SUP>
using Quaternion_t    = Eigen::Matrix<double, 4, 1>;
/// 4x1 Vector [x,y,z,theta]<SUP>T</SUP>
using AxisAngle_t    = Eigen::Matrix<double, 4, 1>;
/// 2x1 Vector [x,y]<SUP>T</SUP>
using Coordinate2D_t  = Eigen::Matrix<double, 2, 1>;
/// 3x1 Vector [x,y,z]<SUP>T</SUP>
using Coordinate3D_t  = Eigen::Matrix<double, 3, 1>;
/// 3x1 Vector [lon,lat,h]<SUP>T</SUP>
using geodethic3D_t   = Eigen::Matrix<double, 3, 1>;

/// Pi constant
const double pi = 3.1415926535897932384626433832795028841971;

template<typename T>
class TransformInterface
{
public:

  /// @brief Set this Transform to Identity.
  virtual void SetIdentity() = 0;

  /// @brief Multiply two Transforms of same Type
  virtual T operator*(const T& other) const = 0;

  /// @brief Multiplication.
  virtual T& operator*=(const T& other) = 0;

  /// @brief Copy assignment.
  virtual T& operator=(const T&) = 0;

  /// @brief Move assignment.
  virtual T& operator=(T&&) = 0;

  /// @brief Eqiality comparison
  virtual bool operator==(const T& other) = 0;

  /// @brief ineqiality comparison
  virtual bool operator!=(const T& other) = 0;

  friend std::ostream& operator<<(std::ostream& os, const T& transform)
  {
    transform.Print(os);
    return os;
  };

private:
  virtual void Print( std::ostream& out ) const = 0;

};


/// @brief Multiplication if the right hand side is a vector of Transforms/Derivatives.
template<typename TL, typename TR>
std::vector<TR> operator*(const TL lhs, const std::vector<TR>& rhs)
{
  static_assert(std::is_base_of<TransformInterface<TL>, TL>::value, "Left hand side is not derived from TransformInterface!!");
  static_assert(std::is_base_of<TransformInterface<TR>, TR>::value, "Right hand side Vector Elements are not derived from TransformInterface!!");

  size_t n_elements = rhs.size();
  std::vector<TR> result(n_elements);
  for(unsigned int i=0; i<n_elements; ++i)
  {
    result[i] = lhs * rhs[i];
  }
  return result;
}
#endif //ABSOLUTE_SLAM_TRANSFORM_INTERFACE_H