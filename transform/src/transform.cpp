/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file transform.cpp
 * @brief Implementations for transformations.
 */

#include "transform.h"
#include <iostream>
#include "string"

// Constructors
TF::TF()
{
  this->r_ <<
  1,0,0,
  0,1,0,
  0,0,1;

  this->t_ <<
  0,0,0;
}

/**
 * @details
 * This Function is adapted from googles ceres library.
 * @param q quaternion
 * @param t translation
 */
TF::TF(Quaternion_t q, Translation_t t)
{
  // This Function is apded from googles ceres library
  // Make convenient names for elements of q.
  double a = q[0];
  double b = q[1];
  double c = q[2];
  double d = q[3];

  // This is not to eliminate common sub-expression, but to
  // make the lines shorter so that they fit in 80 columns!
  double aa = a * a;
  double ab = a * b;
  double ac = a * c;
  double ad = a * d;
  double bb = b * b;
  double bc = b * c;
  double bd = b * d;
  double cc = c * c;
  double cd = c * d;
  double dd = d * d;

  double R[3][3] = {0};

  R[0][0] = aa + bb - cc - dd; R[0][1] = double(2) * (bc - ad);  R[0][2] = double(2) * (ac + bd);  // NOLINT
  R[1][0] = double(2) * (ad + bc);  R[1][1] = aa - bb + cc - dd; R[1][2] = double(2) * (cd - ab);  // NOLINT
  R[2][0] = double(2) * (bd - ac);  R[2][1] = double(2) * (ab + cd);  R[2][2] = aa - bb - cc + dd; // NOLINT

  // Normalize Rotation Matrix
  double normalizer = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
  normalizer = double(1) / normalizer;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      R[i][j] *= normalizer;
    }
  }

  // Copy C array to Eigen Matrix.
  for(int i = 0; i<3 ; ++i)
  {
    for(int j = 0; j<3 ; ++j)
    {
      this->r_(i,j) = R[i][j];
    }
  }

  this->t_ = t;
}


TF::TF(const RotMat_t& r, const Translation_t& t)
{
  CheckProperRotation(r, 1e-5);

  this->r_ = r;
  this->t_ = t;
}

TF::TF(Translation_t t) : TF(TFIdentity().GetRotmatrix(), t){}

TF::TF(RotMat_t r) : TF(r, TFIdentity().GetTranslation()) {}

TF::TF(Quaternion_t q) : TF(q, TFIdentity().GetTranslation()) {}



TF::TF(const TF&) = default;

TF::TF(TF&&) = default;

TF::~TF()
{}

// Getter
RotMat_t TF::GetRotmatrix() const
{
  return this->r_;
}

/**
 *
 * @return Quaternion in Hamiltoneon Convention. [q<SUB>w</SUB>,q<SUB>x</SUB>,q<SUB>y</SUB>,q<SUB>z</SUB>]<SUP>T</SUP>
 */
Quaternion_t TF::GetQuaternion() const
{
  Quaternion_t q;
  double trace = r_(0, 0) + r_(1, 1) + r_(2,  2);

  if( trace > 0 )
  {
    double s = 0.5f / sqrtf(trace+ 1.0f);
    q << 0.25f / s,
         ( r_(2, 1) - r_(1, 2) ) * s,
         ( r_(0, 2) - r_(2, 0) ) * s,
         ( r_(1, 0) - r_(0, 1) ) * s;
  }
  else {
    if ( r_(0, 0) > r_(1, 1) && r_(0, 0) > r_(2, 2) )
    {
      double s = 2.0f * sqrtf( 1.0f + r_(0, 0) - r_(1, 1) - r_(2, 2));
      q << (r_(2, 1) - r_(1, 2) ) / s,
           0.25f * s,
           (r_(0, 1) + r_(1, 0) ) / s,
           (r_(0, 2) + r_(2, 0) ) / s;

    } else if (r_(1, 1) > r_(2, 2))
    {
      double s = 2.0f * sqrtf( 1.0f + r_(1, 1) - r_(0, 0) - r_(2, 2));
      q << (r_(0, 2) - r_(2, 0) ) / s,
           (r_(0, 1) + r_(1, 0) ) / s,
           0.25f * s,
           (r_(1, 2) + r_(2, 1) ) / s;

    } else
    {
      double s = 2.0f * sqrtf( 1.0f + r_(2, 2) - r_(0, 0) - r_(1, 1) );
      q <<  (r_(1, 0) - r_(0, 1) ) / s,
            (r_(0, 2) + r_(2, 0) ) / s,
            (r_(1, 2) + r_(2, 1) ) / s,
            0.25f * s;
    }
  }
  return q;
}

/**
 *
 * @return AxisAngle  [x, y, z]<SUP>T</SUP>
 */
AxisAngle_t TF::GetAxisAngle() const
{
  Quaternion_t q = this->GetQuaternion();
  Eigen::Matrix<double, 3, 1> q_imag;
  q_imag << q[1], q[2], q[3];
  double q_real = q[0];

  Eigen::Matrix<double, 3, 1> axis = q_imag/q_imag.norm();
  double theta = 2 * atan2(q_imag.norm(),q_real);

  AxisAngle_t result;

  result << axis, theta;
  return result;
}

Translation_t TF::GetTranslation() const
{ return t_; }


HomogeniousTF_t TF::GetHomogeniousTF() const
{
  HomogeniousTF_t HTF;
  HTF << this->r_, this->t_,
         0, 0, 0, 1;
  return HTF;
}

void TF::SetIdentity()
{
  *this = TFIdentity();
}

TF TF::Inv() const
{
  TF result;
  result.r_ << this->r_.transpose();
  result.t_ << -1 * result.r_ * this->t_;
  return result;
}

// Operators
TF& TF::operator=(const TF&)=default;

TF& TF::operator=(TF&&)=default;

TF& TF::operator*=(const TF& other)
{
  Translation_t t_cache = this->r_ * other.t_;
  this->t_  = t_cache + this->t_;
  this->r_ *= other.r_;
  return *this;
}

TF TF::operator*(const TF& other) const
{
  TF tf_result = *this;
  tf_result *= other;
  return tf_result;
}

void TF::SetTF(TF other)
{
  this->r_ = other.GetRotmatrix();
  this->t_ = other.GetTranslation();
}

TF TF::GetTF() const
{
  return *this;
}

void TF::Print( std::ostream& out )const
{
  out<<std::scientific;
  out<<"q : ["<<this->GetQuaternion()[0]<<", "<<this->GetQuaternion()[1]<<", "<<this->GetQuaternion()[2]<<", "<<this->GetQuaternion()[3]<<"]"<<std::endl;
  out<<"t : ["<<this->t_.x()<<", "<<this->t_.y()<<", "<<this->t_.z()<<"]"<<std::endl;
  out<<std::defaultfloat;
}

bool TF::operator==(const TF& other)
{
  return this->GetQuaternion() == other.GetQuaternion() && this->GetTranslation() == other.GetTranslation();
}

bool TF::operator!=(const TF& other)
{
  return !(*this == other);
}

/**
@details Returns a TF with zero translation and identity rotation.
@return TF
*/
TF TFIdentity()
{
  RotMat_t r;
  r << 1,0,0,
       0,1,0,
       0,0,1;

  Translation_t t;
  t << 0,0,0;

  return TF(r, t);
}


TF RotX(double rad)
{
  RotMat_t r;
  r<<
  1,        0,         0,
  0, cos(rad), -sin(rad),
  0, sin(rad),  cos(rad);
  return TF(r);
}

TF RotY(double rad)
{
  RotMat_t r;
  r<<
    cos(rad), 0, sin(rad),
           0, 1,        0,
   -sin(rad), 0, cos(rad);
  return TF(r);
}

TF RotZ(double rad)
{
  RotMat_t r;
  r<<
  cos(rad), -sin(rad), 0,
  sin(rad),  cos(rad), 0,
         0,         0, 1;
  return TF(r);
}

/**
@details
@n Formulas by geodethic toolbox [Matlab file Exchange](https://ch.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox)
@n Peter Wasmeier, Technical University of Munich
@n p.wasmeier@bv.tum.de
@n Jan 18, 2006
@n
@n Implemented Elibsoid: besseldhdn

@param[in] lon_lat_h vector with geodethic coordinates.
@return 3D vector with cartesian coordinates in Earth centered Earth Fixed (ECEF) Coordinate Frame.
*/
Coordinate3D_t GeodethicToECEF(const geodethic3D_t& lon_lat_h)
{
  // a : semimajor axis
  // b : semiminor axis
  const double a = 6.377397155000000e+06;
  const double b = 6.356078963003470e+06;

  double lon = lon_lat_h[0];
  double lat = lon_lat_h[1];
  double h = lon_lat_h[2];

  double rho = 180/pi;
  double B = lat/rho;
  double L = lon/rho;

  double e2 = (a*a - b*b)/(a*a);
  double N = a/sqrt(1-e2*sin(B)*sin(B));

  Coordinate3D_t cartesian_ECEF_coordinates;
  cartesian_ECEF_coordinates[0] = (N + h)*cos(B)*cos(L);
  cartesian_ECEF_coordinates[1] = (N + h)*cos(B)*sin(L);
  cartesian_ECEF_coordinates[2] = (N * (1 - e2) + h)*sin(B);

  return cartesian_ECEF_coordinates;
}


/**
@details
@n Formulas by geodethic toolbox [Matlab file Exchange](https://ch.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox)
@n Peter Wasmeier, Technical University of Munich
@n p.wasmeier@bv.tum.de
@n Jan 18, 2006
@n
@n Implemented Elibsoid: besseldhdn

@param[in] xyz 3D vector with cartesian coordinates in Earth centered Earth Fixed (ECEF) Coordinate Frame.
@return Vector with geodethic coordinates.
*/
geodethic3D_t ECEFToGeodethic(const Coordinate3D_t& xyz)
{
  const double a = 6.377397155000000e+06;
  const double b = 6.356078963003470e+06;

  // lon [-180 180]
  double lon = atan2(xyz.y(), xyz.x()) * 180/pi;

  double B0 = atan2(xyz.z(), sqrt(xyz.x() * xyz.x() + xyz.y() * xyz.y()));
  double B = 100.01;
  double e2 = (a*a-b*b)/(a*a);
  double h = 0;
  double N = 0;
  while(fabs(B-B0) > 1e-10)
  {
    N = a/sqrt(1 - e2 * sin(B0) * sin(B0));
    h = sqrt(xyz.x() * xyz.x() + xyz.y() * xyz.y()) / cos(B0) - N;
    B = B0;
    B0 = atan((xyz.z()/sqrt(xyz.x() * xyz.x() + xyz.y() * xyz.y())) * 1/(1-e2*N/(N+h)));
  }
  double lat = B * 180 / pi;
  geodethic3D_t lon_lat_h; lon_lat_h << lon, lat, h;

  return  lon_lat_h;
}


/**
@param tf1 Expressed in same coordinate frame than tf2.
@param tf2 Expressed in same coordinate frame than tf1.
@return Euclidian distance between the two points given as argument.
*/
double GetDistanceTranslation(const Translation_t & t1, const Translation_t & t2)
{
  Coordinate3D_t diff = t1 - t2;
  return diff.norm();
}

/**
 * Distance of two rotations is calculated. The distance is defined as the minimal absolute angle in radian between two orientations r1 and r2.
 * The rotation axis of the difference rotation r_error is chosen such that the rotation angle is minimized.
 *
 * @param r1[in] Rotation Matrix 1
 * @param r2[in] Rotation Matrix 2
 * @return Distance in Radian.
 */
double GetDistanceRotation(const Quaternion_t& q1, const Quaternion_t & q2)
{
  TF tf1(q1);
  TF tf2(q2);

  TF rot_distance = tf1.Inv() * tf2;

  AxisAngle_t axis_angle = rot_distance.GetAxisAngle();
  return axis_angle[3];
}

bool IsProperRotation(RotMat_t r, float tol)
{
  RotMat_t eye;
  eye <<
  1,0,0,
  0,1,0,
  0,0,1;

  RotMat_t check_orthogonality1 = r.transpose() * r - eye;
  RotMat_t check_orthogonality2 = r * r.transpose() - eye;
  double check_scaling = r.determinant() - 1;

 return (check_orthogonality1.norm() < tol) && (check_orthogonality2.norm() < tol) && (std::fabs(check_scaling) < tol);
}

/**
 *
 * @param t1 First Vector
 * @param t2 Second Vector
 * @return angle between the two vectors in radian.
 */
double GetAngleBetweenVectors(const Translation_t t1, const Translation_t t2)
{
  double dot_product = t1[0] * t2[0] + t1[1] * t2[1] + t1[2] * t2[2];
  double squared_len1 = t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]; double len1 = sqrt(squared_len1);
  double squared_len2 = t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]; double len2 = sqrt(squared_len2);
  return acos(dot_product / (len1 * len2));
}

void CheckProperRotation(RotMat_t r, float tol)
{
  double check_scaling = r.determinant() - 1;
  if (!IsProperRotation(r, tol))
  {
    if (std::fabs(check_scaling) > tol)
    {
      throw std::invalid_argument("Error: Rotation Matrix does not describe a proper Rotation. There is scaling or a flipped axis involved. ( fabs(r.determinant() - 1) > tol ).");
    }
    else
    {
      throw std::invalid_argument("Error: Rotation Matrix does not describe a proper Rotation. The matrix is not orthogonal. ( norm(r.transpose() * r) - eye > tol or norm(r * r.transpose()) - eye > tol )");
    }
  }
}


