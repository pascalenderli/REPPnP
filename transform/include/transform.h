/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file transform.h
 * @brief Function prototypes for transformations
 */

#ifndef ABSOLUTE_SLAM_TRANSFORM_H
#define ABSOLUTE_SLAM_TRANSFORM_H


#include "transform_interface.h"

/**
@brief Tools for geometrical transformations in 3D space.

@details
A combination of a translation and a rotation forms a transform (TF). @n
A TF can represent Transformations between Coordinate frames. @n
Beside of that they can be useful to Store a position and an orientation of objects in 3D space. @n
@n
@n All coordinates (parameters and returns) are stored in Eigen Matrices.
#### Notation
Represents a Transformation from B to A Frame. @n
T<SUB>AB</SUB> @n
@n
P is a point which is represented in C frame. @n
<SUB>C</SUB>P @n
@n
Point coordinates (and its orientation) of point P are transformed from B to A Frame representation. @n
<SUB>A</SUB>P = T<SUB>AB</SUB> * <SUB>B</SUB>P  @n
*/
class TF : public TransformInterface<TF>
{

public:

  /// @brief Default Constructor.
  TF();

  /// @brief Copy Constructor.
  TF(const TF&);

  /// @brief Move Constructor.
  TF(TF&&);

  /// @brief Construct TF with 4x1 Quaternion and 3x1 Translation Vector.
  TF(Quaternion_t q, Translation_t t);

  /// @brief Construct TF with 3x3 Roration Matrix and 3x1 Translation Vector.
  TF(const RotMat_t& r, const Translation_t& t);

  /// @brief Construct a Point with identity rotation.
  TF(Coordinate3D_t t);

  /// @brief Construct a Rotation with zero Translation.
  TF(RotMat_t r);

  /// @brief Construct a Rotation with zero Translation.
  TF(Quaternion_t q);

  /// @brief Default Destructor.
  virtual ~TF();

  // Getters
  RotMat_t GetRotmatrix() const;
  Quaternion_t GetQuaternion() const;
  AxisAngle_t GetAxisAngle() const;
  Translation_t GetTranslation() const;
  HomogeniousTF_t GetHomogeniousTF() const;

  // Setters
  /// @brief Set this Transform to identity (Zero translation / identity rotation).
  void SetIdentity() override;

  // Functionality
  /// @brief Get the inverse of a TF.
  TF Inv() const;

  // Operators
  /// @brief Copy assignment.
  TF& operator=(const TF&) override;

  /// @brief Move assignment.
  TF& operator=(TF&&) override;

  /// @brief Multiplication.
  TF operator*(const TF& other) const override;

  /// @brief Multiplication.
  TF& operator*=(const TF& other) override;

  /// @brief Eqiality comparison
  bool operator==(const TF& other) override;

  /// @brief ineqiality comparison
  bool operator!=(const TF& other) override;

  /// @brief formatting class members for output stream in yaml format.
  void Print( std::ostream& out )const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  /// @brief Set all members of the Transform.
  virtual void SetTF(TF other);
  virtual TF GetTF() const;

private:
  RotMat_t r_;
  Translation_t t_;
};

/// @brief Get the Identity Transform.
TF TFIdentity();

/// @brief Get Rotation around x axis.
TF RotX(double rad);

/// @brief Get Rotation around y axis.
TF RotY(double rad);

/// @brief Get Rotation around z axis.
TF RotZ(double rad);

/// @brief Transforms Geodethic coordinates.
Coordinate3D_t GeodethicToECEF(const geodethic3D_t& lon_lat_h);

/// @brief Transforms Geodethic coordinates.
geodethic3D_t ECEFToGeodethic(const Coordinate3D_t& xyz);

/// @brief Calculate the euclidean distance between two points.
double GetDistanceTranslation(const Translation_t& t1, const Translation_t& t2);

/// @brief Calculate error between two rotations.
double GetDistanceRotation(const Quaternion_t& q1, const Quaternion_t & q2);

/// @brief Calculate angle between two Vectors.
double GetAngleBetweenVectors(const Translation_t t1, const Translation_t t2);

/// @brief Checks if a rotation matrix is a proper rotation or not.
bool IsProperRotation(RotMat_t r, float tol=1e-6);

/// @brief Throws an exception if r is not a proper rotation.
void CheckProperRotation(RotMat_t r, float tol=1e-6);

#endif //ABSOLUTE_SLAM_TRANSFORM_H
