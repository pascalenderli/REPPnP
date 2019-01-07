/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file transform.cpp
 * @brief Implementations for transformations.
 */

#include "scaled_transform.h"
#include <iostream>

// Derived Class ScaledTF:
// Constructors
ScaledTF::ScaledTF()
{this->scale_ = 1; }

ScaledTF::ScaledTF(const ScaledTF&) = default;

ScaledTF::ScaledTF::ScaledTF(ScaledTF&&) = default;

ScaledTF::ScaledTF(RotMat_t r, Translation_t t, double scale) : TF(r, t), scale_(scale){}

ScaledTF::ScaledTF(Quaternion_t q, Translation_t t, double scale) : TF(q, t), scale_(scale){}

ScaledTF::ScaledTF(Quaternion_t q, double scale) : TF(q), scale_(scale){}

ScaledTF::ScaledTF(RotMat_t r, double scale) : TF(r), scale_(scale){}

ScaledTF::ScaledTF(Translation_t t, double scale) : TF(t), scale_(scale){}

ScaledTF::ScaledTF(TF tf, double scale) : TF(tf), scale_(scale){}

ScaledTF::ScaledTF(double scale) : TF(TFIdentity()), scale_(scale) {}

double ScaledTF::GetScale() const
{ return this->scale_; }

// Setters
void ScaledTF::SetIdentity()
{
  TF::SetIdentity();
  this->scale_ = 1.0;
}

// Functionality
ScaledTF ScaledTF::Inv() const
{
  Eigen::Matrix<double, 3, 3> I;
  double s_inv = 1/(this->scale_);
  RotMat_t r_inv = this->GetRotmatrix().transpose();
  Translation_t t_inv = -1 * s_inv * r_inv * this->GetTranslation();
  return ScaledTF(r_inv, t_inv, s_inv);
}

// Operators
ScaledTF& ScaledTF::operator=(const ScaledTF&) = default;

ScaledTF& ScaledTF::operator=(ScaledTF&&) = default;

ScaledTF& ScaledTF::operator*=(const ScaledTF& other)
{
  this->SetTF( this->GetTF() * other.GetTF());
  this->scale_ *= other.scale_;
  return *this;
}

ScaledTF ScaledTF::operator*(const ScaledTF& other) const
{
  ScaledTF scaled_tf_result = *this;
  scaled_tf_result *= other;
  return scaled_tf_result;
}

TF ScaledTF::operator*(const TF& other) const
{
  RotMat_t r = this->GetRotmatrix() * other.GetRotmatrix();
  Translation_t t = (this->GetRotmatrix() * other.GetTranslation());
  t = (this->scale_ * t);
  t = t + this->GetTranslation();

  return TF(r, t);
}

ScaledTF ScaledTF::GetScaledTF()
{
  return *this;
}

bool ScaledTF::operator==(const ScaledTF& other)
{
  return this->GetRotmatrix() == other.GetRotmatrix() && this->GetTranslation() == other.GetTranslation() && this->GetScale() == other.GetScale();
}

bool ScaledTF::operator!=(const ScaledTF& other)
{
  return !(*this == other);
}

void ScaledTF::Print( std::ostream& out ) const
{
  out<<this->GetTF();
  out<<"s : "<<this->scale_<<std::endl;
  out<<std::defaultfloat;
}

/// @brief Get the Identity Transform.
ScaledTF ScaledTFIdentity()
{
  ScaledTF result;
  result.SetIdentity();
  return result;
}
