/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file transform.h
 * @brief Function prototypes for transformations
 */

#ifndef ABSOLUTE_SLAM_SCALED_TRANSFORM_H
#define ABSOLUTE_SLAM_SCALED_TRANSFORM_H

#include "transform.h"

/// @brief Orthogonal Transformation with scaling proberty (Non Rigid but non Skewed Transforms). Inherited form TF.
class ScaledTF: public TF, TransformInterface<ScaledTF>
{

public:
  /// @brief Default Constructor.
  ScaledTF();

  /// @brief Copy Constructor.
  ScaledTF(const ScaledTF&);

  /// @brief Move Constructor.
  ScaledTF(ScaledTF&&);

  /// @brief Construct ScaledTF from Rotation Matrix, Translation and scale.
  ScaledTF(RotMat_t r, Translation_t t, double scale);

  /// @brief Construct ScaledTF from Quaternion, Translation and scale.
  ScaledTF(Quaternion_t q, Translation_t t, double scale);

  /// @brief Construct ScaledTF from Quaternion, Zero Translation and scale.
  ScaledTF(Quaternion_t q, double scale);

  /// @brief Construct ScaledTF from Rotation Matrix, Zero Translation and scale.
  ScaledTF(RotMat_t r, double scale);

  /// @brief Construct ScaledTF from Unit Rotation, Translation and scale.
  ScaledTF(Translation_t t, double scale);

  /// @brief Construct ScaledTF from Baseclass object and scale.
  ScaledTF(TF tf, double scale);

  /// @brief Construct ScaledTF with Identity Rotation/Translation but with a scale.
  ScaledTF(double scale);

  // Setters
  /// @brief Set this Transform to identity (Zero translation / identity rotation / identity scale).
  void SetIdentity() override;

  // Functionality
  /// @brief Get the inverse of a ScaledTF. such that ScaledTF*inv(ScaledTF) = Identity
  ScaledTF Inv() const;

  // Operators
  /// @brief Copy assignment.
  ScaledTF& operator=(const ScaledTF&) override;

  /// @brief Move assignment.
  ScaledTF& operator=(ScaledTF&&) override;

  /// @brief Multiplication.
  ScaledTF operator*(const ScaledTF& other) const override;

  /// @brief Multiplication.
  ScaledTF& operator*=(const ScaledTF& other) override;

  /// @brief Multiplication ScaledTF * TF (!! not in the interface so vector multiplication does not work)
  TF operator*(const TF& other) const override;

  /// @brief Eqiality comparison
  bool operator==(const ScaledTF& other) override;

  /// @brief ineqiality comparison
  bool operator!=(const ScaledTF& other) override;

  double GetScale() const;

  /// @brief formatting class members for output stream in yaml format.
  void Print( std::ostream& out ) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
protected:
  /// @brief Returns a copy of a scaled Transform.
  ScaledTF GetScaledTF();

private:
  double scale_;
};

// Free Functions
/// @brief Get the Identity Transform.
ScaledTF ScaledTFIdentity();

#endif //ABSOLUTE_SLAM_SCALED_TRANSFORM_H
