/**
 * @author Pascal Enderli
 * @date 04.10.2018
 * @file test_transform.h
 * @brief Test Cases For Transformation Functionalities.
 *
 *  //Install cmake:
 *  sudo apt-get install cmake
 *
 *  Install Google Test:
 *  sudo apt-get install libgtest-dev
 *  cd /usr/src/gtest
 *  sudo mkdir build
 *  cd build
 *  sudo cmake ..
 *  sudo make
 *  sudo ln -s *.a /usr/lib
 *
 *
 *  //Install Google Mock: (Not Needed for TF and derivatives.)
 *  sudo apt-get install google-mock
 *  cd /usr/src/gmock
 *  sudo mkdir build
 *  cd build
 *  sudo cmake ..
 *  sudo make
 *  sudo ln -s *.a /usr/lib
 *
 *  in
 *
 *
 *  //Native GTest:
 *
 *  ASSERT_TRUE(condition); 	     EXPECT_TRUE(condition); 	     condition is true
 *  ASSERT_FALSE(condition); 	     EXPECT_FALSE(condition); 	   condition is false
 *
 *  ASSERT_EQ(val1, val2);       	 EXPECT_EQ(val1, val2);  	      val1 == val2
 *  ASSERT_NE(val1, val2);       	 EXPECT_NE(val1, val2);  	      val1 != val2
 *  ASSERT_LT(val1, val2);       	 EXPECT_LT(val1, val2);  	      val1 < val2
 *  ASSERT_LE(val1, val2);       	 EXPECT_LE(val1, val2);  	      val1 <= val2
 *  ASSERT_GT(val1, val2);       	 EXPECT_GT(val1, val2);  	      val1 > val2
 *  ASSERT_GE(val1, val2);       	 EXPECT_GE(val1, val2);  	      val1 >= val2
 *
 *  ASSERT_STREQ(str1, str2); 	   EXPECT_STREQ(str1, str2); 	    the two C strings have the same content
 *  ASSERT_STRNE(str1, str2); 	   EXPECT_STRNE(str1, str2);    	the two C strings have different contents
 *  ASSERT_STRCASEEQ(str1, str2);  EXPECT_STRCASEEQ(str1, str2); 	the two C strings have the same content, ignoring case
 *  ASSERT_STRCASENE(str1, str2);  EXPECT_STRCASENE(str1, str2); 	the two C strings have different contents, ignoring case
 *
 *  ASSERT_FLOAT_EQ(val1, val2); 	       EXPECT_FLOAT_EQ(val1,val2); 	        the two float values are almost equal
 *  ASSERT_DOUBLE_EQ(val1, val2); 	     EXPECT_DOUBLE_EQ(val1, val2); 	      the two double values are almost equal
 *  ASSERT_NEAR(val1, val2, abs_error);  EXPECT_NEAR(val1, val2, abs_error);  the difference between val1 and val2 doesn't exceed the given absolute error
 *
 *
 *  //From eigen-checks.h (by ASL, ETHZ):
 *
 *  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(MatrixA, MatrixB))
 *  Succeeds if two matrices are binary equal.
 *
 *  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(MatrixA, MatrixB))
 *  Succeeds if two matrices are equal to floating-point precision
 *
 *  EXPECT_TRUE(EIGEN_MATRIX_NEAR(MatrixA, MatrixB, Precision))
 *  Succeeds if two matrices are equal to a user-specified precision.
 *
 *  EXPECT_TRUE(EIGEN_MATRIX_ZERO(MatrixA, Precision))
 *  Succeeds if a matrix is equal to zero to a user-specified precision.
 *
 *  ASSERT_TRUE(EIGEN_MATRIX_EQUAL(MatrixA, MatrixB));
 *  ASSERT_FALSE(EIGEN_MATRIX_EQUAL(MatrixA, MatrixB));
 *
*/

#include <eigen-checks/gtest.h>
#include <gtest/gtest.h>
#include <iostream>
#include "transform_interface.h"
#include "scaled_transform.h"
#include "string.h"
#include <vector>



TEST(RigidTransformation, QuaternionConstructor)
{

  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
             0.9523810, 0.2380952,-0.1904762,
            -0.1904762, 0.9523810, 0.2380952;
  t_init << 1, 1, 1;
  TF tf(r_init , t_init );

  Quaternion_t q_expected;
  q_expected << 1.5, 1, 1, 1;
  q_expected = q_expected.normalized();
  Quaternion_t q_test = tf.GetQuaternion();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(q_expected, q_test, 10e-8));
}


TEST(RigidTransformation, RotationMatrixConstructor)
{
  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
             0.9523810, 0.2380952,-0.1904762,
            -0.1904762, 0.9523810, 0.2380952;

  t_init << 1, 1, 1;
  TF tf(r_init , t_init);

  RotMat_t r_test = tf.GetRotmatrix();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(r_init, r_test, 10e-8));
}


TEST(RigidTransformation, TranslationGetter)
{
  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
             0.9523810, 0.2380952,-0.1904762,
            -0.1904762, 0.9523810, 0.2380952;

  t_init << 1, 1, 1;
  TF tf(r_init , t_init);

  Translation_t t_test = tf.GetTranslation();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(t_init, t_test, 10e-8));
}


TEST(RigidTransformation, HomogeniousGetter1)
{
  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
             0.9523810, 0.2380952,-0.1904762,
            -0.1904762, 0.9523810, 0.2380952;

  t_init << 1, 1, 1;
  TF tf(r_init , t_init);

  HomogeniousTF_t tf_expected;
  tf_expected << r_init , t_init,
                 0, 0, 0, 1;

  HomogeniousTF_t tf_test = tf.GetHomogeniousTF();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(tf_expected, tf_test, 10e-8));
}


TEST(RigidTransformation, QuaternionGetter)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init  << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf(q_init , t_init );

  q_init = q_init.normalized();
  Quaternion_t q_test = tf.GetQuaternion();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(q_init, q_test, 10e-8));
}


TEST(RigidTransformation, RotationMatrixGetter)
{
  Quaternion_t q_init;
  Translation_t t_init;
  q_init  << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf(q_init , t_init);

  RotMat_t r_expected;
  r_expected << 1, 0, 0,
                0, 0,-1,
                0, 1, 0;

  RotMat_t r_test = tf.GetRotmatrix();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(r_expected, r_test, 10e-8));
}


TEST(RigidTransformation, TranslationGetter2)
{
  Quaternion_t q_init;
  Translation_t t_init;
  q_init  << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf(q_init , t_init);

  Translation_t t_test = tf.GetTranslation();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(t_init, t_test, 10e-8));
}


TEST(RigidTransformation, HomogeniousGetter2)
{
  Quaternion_t q_init;
  Translation_t t_init;
  q_init  << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf(q_init , t_init);

  HomogeniousTF_t tf_expected;
  tf_expected << 1, 0, 0, 1,
                 0, 0,-1, 1,
                 0, 1, 0, 1,
                 0, 0, 0, 1;

  HomogeniousTF_t tf_test = tf.GetHomogeniousTF();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(tf_expected, tf_test, 10e-8));
}


TEST(RigidTransformation, GetIdentity)
{
  TF tf_test = TFIdentity();
  Translation_t t_expected;
  Quaternion_t q_expected;
  RotMat_t r_expected;
  HomogeniousTF_t tf_homogen_expected;
  q_expected << 1, 0, 0, 0;
  r_expected << 1,0,0,
                0,1,0,
                0,0,1;
  t_expected << 0,0,0;
  tf_homogen_expected << r_expected, t_expected,
                         0,0,0,1;

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(q_expected, tf_test.GetQuaternion()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(r_expected, tf_test.GetRotmatrix()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(t_expected, tf_test.GetTranslation()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(tf_homogen_expected, tf_test.GetHomogeniousTF()));
}

TEST(RigidTransformation, InverseGetter)
{
  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
             0.9523810, 0.2380952,-0.1904762,
            -0.1904762, 0.9523810, 0.2380952;

  t_init << 5, -1, 10;
  TF tf_init(r_init , t_init);

  Quaternion_t q_expected;
  RotMat_t r_expected;
  Translation_t t_expected;
  HomogeniousTF_t tf_homogen_expected;

  q_expected << 0.654653676942776, -0.436435777354585, -0.436435777354585, -0.436435777354585;
  r_expected << 0.238095210884357,  0.952380929251700, -0.190476140136057,
               -0.190476140136057,  0.238095210884357,  0.952380929251700,
                0.952380929251700, -0.190476140136057,  0.238095210884357;
  t_expected << 1.666666276190484, -8.333333380952354, -7.333332895238130;

  tf_homogen_expected << r_expected, t_expected,
                         0, 0, 0, 1;

  TF tf_test = tf_init.Inv();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(q_expected.normalized(), tf_test.GetQuaternion()  , 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(r_expected, tf_test.GetRotmatrix()  , 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(t_expected, tf_test.GetTranslation(), 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(tf_homogen_expected, tf_test.GetHomogeniousTF(), 10e-7));
}


// Free Functions
TEST(FreeFunction, TFIdentityGetter)
{
  Quaternion_t q_expected;
  Translation_t t_expected;
  q_expected << 1, 0, 0, 0;
  t_expected << 0, 0, 0;

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(q_expected, TFIdentity().GetQuaternion()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(t_expected, TFIdentity().GetTranslation()));
}

TEST(FreeFunction, ConvertGeodethic)
{
  geodethic3D_t zurich_lon_lat_h_init; zurich_lon_lat_h_init << 8.5391825, 47.3686498, 429;
  Coordinate3D_t zurich_xyz = GeodethicToECEF(zurich_lon_lat_h_init);
  geodethic3D_t zurich_lon_lat_h_test = ECEFToGeodethic(zurich_xyz);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(zurich_lon_lat_h_init, zurich_lon_lat_h_test  , 10e-6));
}

TEST(FreeFunction, CalculateVectorDistance)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf_1(q_init , t_init);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  TF tf_2(q_init , t_init);

  EXPECT_DOUBLE_EQ( sqrt(5), GetDistanceTranslation(tf_1.GetTranslation(), tf_2.GetTranslation()));
}

TEST(FreeFunction, BasicRotations)
{
  TF result = RotX(1) * RotY(2) * RotZ(3);
  RotMat_t r_expect;
  r_expect<<
    0.4120,  0.0587,  0.9093,
   -0.6812, -0.6429,  0.3502,
    0.6051, -0.7637, -0.2248;
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(r_expect, result.GetRotmatrix(), 10e-4));
}

TEST(RigidTransformation, Multiply)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf_1(q_init , t_init);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  TF tf_2(q_init , t_init);

  TF tf_result = tf_1 * tf_2;

  HomogeniousTF_t expected_homogeneous = tf_1.GetHomogeniousTF() * tf_2.GetHomogeniousTF();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(expected_homogeneous, tf_result.GetHomogeniousTF(), 10e-4));
}

TEST(RigidTransformation, Equality)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf_1(q_init , t_init);
  TF tf_2(q_init , t_init);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  TF tf_3(q_init , t_init);


  EXPECT_TRUE(tf_1 == tf_2);
  EXPECT_TRUE(tf_1 != tf_3);
}

TEST(FreeFunction, BoxMinus)
{
  Quaternion_t q_init ;
  Translation_t t_init;

  double expected_angle = 0.5;

  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF tf_1(q_init , t_init);

  q_init << (1 + cos(expected_angle/2)), (1 * sin(expected_angle/2)), (0 * sin(expected_angle/2)), (0 * sin(expected_angle/2));
  t_init << 1, 2, 3;
  TF tf_2(q_init , t_init);

  float test_angle = GetDistanceRotation(tf_1.GetQuaternion(), tf_2.GetQuaternion());

  EXPECT_EQ(test_angle, test_angle);
}

TEST(RigidTransformation, SetIdentity)
{
  Quaternion_t q_init ;
  Translation_t t_init;

  q_init << 7, 4, 5, 6;
  t_init << 1, 2, 1;
  TF tf(q_init , t_init);

  tf.SetIdentity();

  EXPECT_TRUE(tf == TFIdentity());
}

TEST(Check, ProperRotation)
{
  RotMat_t r_wrong;
  r_wrong<<
   1, 0, 0.1,
   0, 1,   0,
   0, 0,   1;

  EXPECT_TRUE( IsProperRotation(RotX(0.5).GetRotmatrix()) );
  EXPECT_FALSE(IsProperRotation(r_wrong));
}


// Test ScaledTF
TEST(ScaledTF, ScaleGetter)
{
  ScaledTF s_tf(2);
  EXPECT_EQ(2, s_tf.GetScale());

  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
      0.9523810, 0.2380952,-0.1904762,
      -0.1904762, 0.9523810, 0.2380952;
  t_init << 1, 1, 1;
  TF tf(r_init , t_init );

  ScaledTF s_tf2(tf, 2.1);
  EXPECT_EQ(2.1, s_tf2.GetScale());
}


TEST(ScaledTF, SetIdentity)
{
  Translation_t t_init;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
      0.9523810, 0.2380952,-0.1904762,
      -0.1904762, 0.9523810, 0.2380952;
  t_init << 1, 1, 1;
  TF tf(r_init , t_init );
  ScaledTF s_tf(tf, 2.1);

  s_tf.SetIdentity();

  Translation_t t_expected;
  Quaternion_t q_expected;
  RotMat_t r_expected;
  HomogeniousTF_t tf_homogen_expected;
  q_expected << 1, 0, 0, 0;
  r_expected << 1,0,0,
      0,1,0,
      0,0,1;
  t_expected << 0,0,0;
  tf_homogen_expected << r_expected, t_expected,
      0,0,0,1;


  EXPECT_EQ(1, s_tf.GetScale());
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(q_expected, s_tf.GetQuaternion()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(r_expected, s_tf.GetRotmatrix()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(t_expected, s_tf.GetTranslation()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(tf_homogen_expected, s_tf.GetHomogeniousTF()));
}

TEST(ScaledTF, Inverse)
{
  Translation_t t_init;
  float scale_init = 2.1;
  RotMat_t r_init;
  r_init <<  0.2380952,-0.1904762, 0.9523810,
      0.9523810, 0.2380952,-0.1904762,
      -0.1904762, 0.9523810, 0.2380952;

  t_init << 5, -1, 10;
  TF tf_init(r_init , t_init);

  Quaternion_t q_expected;
  RotMat_t r_expected;
  Translation_t t_expected;
  HomogeniousTF_t tf_homogen_expected;

  q_expected << 0.654653676942776, -0.436435777354585, -0.436435777354585, -0.436435777354585;
  r_expected << 0.238095210884357,  0.952380929251700, -0.190476140136057,
      -0.190476140136057,  0.238095210884357,  0.952380929251700,
      0.952380929251700, -0.190476140136057,  0.238095210884357;
  t_expected << -1/scale_init * r_expected * t_init;

  tf_homogen_expected << r_expected, t_expected,
      0, 0, 0, 1;

  ScaledTF s_tf_test = ScaledTF(tf_init, scale_init).Inv();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(q_expected.normalized(), s_tf_test.GetQuaternion()  , 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(r_expected, s_tf_test.GetRotmatrix()  , 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(t_expected, s_tf_test.GetTranslation(), 10e-7));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(tf_homogen_expected, s_tf_test.GetHomogeniousTF(), 10e-7));

  EXPECT_FLOAT_EQ(1/scale_init, s_tf_test.GetScale());
}

TEST(ScaledTF, Multiply_ScaledTFxScaledTF)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  ScaledTF s_tf_1(TF(q_init , t_init), 2.1);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  ScaledTF s_tf_2(TF(q_init , t_init), 100);

  ScaledTF s_tf_test = s_tf_1 * s_tf_2;

  HomogeniousTF_t expected_homogeneous = s_tf_1.GetHomogeniousTF() * s_tf_2.GetHomogeniousTF();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(expected_homogeneous, s_tf_test.GetHomogeniousTF(), 10e-4));

  EXPECT_EQ(2.1 * 100, s_tf_test.GetScale());
}

TEST(ScaledTF, Multiply_ScaledTFxRigidScaledTF)
{
  double scaler = 100;

  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  TF rigid_tf_1(q_init , t_init);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  TF rigid_tf_2(q_init , t_init);
  ScaledTF scaled_tf_2(TF(q_init , t_init), scaler);

  TF saling_test_result = scaled_tf_2 * rigid_tf_1;
  TF rigid_test_result  =  rigid_tf_2 * rigid_tf_1;

  Translation_t expected_translation = scaler * rigid_tf_2.GetRotmatrix() * rigid_tf_1.GetTranslation() + rigid_tf_2.GetTranslation();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(rigid_test_result.GetRotmatrix(), saling_test_result.GetRotmatrix(), 10e-7));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(expected_translation, saling_test_result.GetTranslation(), 10e-7));
}

TEST(ScaledTF, Equality)
{
  Quaternion_t q_init ;
  Translation_t t_init;
  q_init << 1, 1, 0, 0;
  t_init << 1, 1, 1;
  ScaledTF s_tf_1(TF(q_init , t_init), 1);
  ScaledTF s_tf_2(TF(q_init , t_init), 1);

  q_init << 1.5, 8, 9, 1;
  t_init << 1, 2, 3;
  ScaledTF s_tf_3(TF(q_init , t_init), 2.2);


  EXPECT_TRUE(s_tf_1 == s_tf_2);
  EXPECT_TRUE(s_tf_1 != s_tf_3);
}

TEST(FreeFunction, ScaledTFIdentityGetter)
{
  Quaternion_t q_expected;
  Translation_t t_expected;
  q_expected << 1, 0, 0, 0;
  t_expected << 0, 0, 0;

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(q_expected, ScaledTFIdentity().GetQuaternion()));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(t_expected, ScaledTFIdentity().GetTranslation()));
  EXPECT_EQ(1, ScaledTFIdentity().GetScale());
}

TEST(FreeInterfaceFunction, Transform_times_vecotor)
{
  Quaternion_t q;
  Translation_t t;
  q << 1, 0, 0, 0;
  t << 0, 0, 0;

  TF tf1(q, t);
  q << 1, 1, 2, 2;
  t << 1, 5, 7;
  Quaternion_t q_expected = q;
  Translation_t t_expected = t;

  TF tf2(q, t);
  q << 8, 5, 4, 1;
  t << 0, 0, 0;
  TF tf3(q, t);

  TF mytfs[] = {tf1, tf2, tf3};
  std::vector<TF> tf_vector(mytfs, mytfs + 3);
  std::vector<TF> tf_vector_test_result = tf1 * tf_vector;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(q_expected.normalized(), tf_vector_test_result[1].GetQuaternion(), 1e-7));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL(t_expected, tf_vector_test_result[1].GetTranslation()));
}

TEST(FreeFunction, AngleBetweenVectors)
{
  Translation_t t1;
  t1 << 1, 0, 0;

  Translation_t t2;
  t2 << 0, 0.5, 0;


  EXPECT_DOUBLE_EQ(pi/2, GetAngleBetweenVectors( t1, t2));
}



int main(int argc, char **argv)
{

  Quaternion_t q;
  Translation_t t;
  q << 1, 0, 0, 0;
  t << 0, 0, 0;
  std::cout<<"Print Transforms to Console:"<<std::endl;
  std::cout<<"TF: "<<std::endl<<TF(q,t) <<std::endl<<std::endl;
  std::cout<<"Scaled TF:"<<std::endl<<ScaledTF(q,t,1) <<std::endl;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
