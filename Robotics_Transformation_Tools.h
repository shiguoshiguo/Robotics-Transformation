#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
// quaternion rotation etc.
#include <Eigen/Geometry>

using namespace std;

class QuaternionPose;
class EulerPose;
class LeicaPose;
class KUKAPose;


class QuaternionPose
{
public:
	Eigen::Vector3d position;
	Eigen::Quaterniond gesture;
public:
	QuaternionPose();
	QuaternionPose(const Eigen::Vector3d& pos, const Eigen::Quaterniond& q);
	QuaternionPose(double px, double py, double pz, double qw, double qx, double qy, double qz);
	QuaternionPose(double px, double py, double pz);
	QuaternionPose(double qw, double qx, double qy, double qz);
	QuaternionPose(const Eigen::Vector3d& pos);
	QuaternionPose(const Eigen::Quaterniond& q);	
	QuaternionPose(const Eigen::Matrix3d& m);
	QuaternionPose(const Eigen::Isometry3d& isom);
	QuaternionPose(const Eigen::Vector3d& pos, const Eigen::Vector3d& ges);
	QuaternionPose(const Eigen::Vector3d& pos, const Eigen::Matrix3d& m);
	QuaternionPose(const EulerPose& ep);
	Eigen::Isometry3d IsometryMatrix(void) const;
	Eigen::Vector3d operator*(const Eigen::Vector3d& vec);
	friend std::ostream& operator<< (std::ostream& ostr, const QuaternionPose qp);
};

class EulerPose
{
public:
	Eigen::Vector3d position;
	Eigen::Vector3d gesture;
	//const string roration_seq = "ZYX";
public:
	EulerPose();
	EulerPose(const Eigen::Vector3d& pvec, const Eigen::Vector3d& gvec);
	EulerPose(double x, double y, double z, double a, double b, double c);
	EulerPose(double n1, double n2, double n3, int tag);
	EulerPose(const Eigen::Vector3d& v, int tag);
	EulerPose(const Eigen::Matrix3d& m);
	EulerPose(const Eigen::Isometry3d& isom);
	EulerPose(const Eigen::Vector3d& pos, const Eigen::Matrix3d& m);
	EulerPose(const Eigen::Vector3d& pos, const Eigen::Quaterniond& q);
	EulerPose(const QuaternionPose& qp);
	Eigen::Isometry3d IsometryMatrix(void) const;
	Eigen::Vector3d operator*(const Eigen::Vector3d& vec);
	friend std::ostream& operator<< (std::ostream& ostr, const EulerPose ep);
};



class LeicaPose : public QuaternionPose
{
public:
	using QuaternionPose::QuaternionPose;
	LeicaPose();
	LeicaPose(double h, double v, double d, const Eigen::Quaterniond& q);
	LeicaPose(double h, double v, double d, double w, double x, double y, double z);
	double H() const;
	double V() const;
	double D() const;
	double q0() const;
	double q1() const;
	double q2() const;
	double q3() const;
	friend std::ostream& operator<< (std::ostream& ostr, const QuaternionPose lp);
		
};


class KUKAPose : public EulerPose
{
public:
	using EulerPose::EulerPose;
	KUKAPose();
	KUKAPose(const Eigen::Vector3d& pvec, const Eigen::Vector3d& gvec);
	KUKAPose(double x, double y, double z, double a, double b, double c);
	KUKAPose(double n1, double n2, double n3, int tag);
	KUKAPose(const Eigen::Vector3d& v, int tag);
	double x() const;
	double y() const;
	double z() const;
	double a() const;
	double b() const;
	double c() const;
	friend std::ostream& operator<< (std::ostream& ostr, const KUKAPose ep);
};


// basic transform functions
// left-hand ball coor and rectangle coor
void xyz2lhvd(const double x, const double y, const double z, double& H, double & V, double& D);
void lhvd2xyz(const double h, const double v, const double d, double& x, double & y, double& z);
// qunternion and rotation matrix
Eigen::Matrix3d q2r(const Eigen::Quaterniond& q);
Eigen::Quaterniond r2q(const Eigen::Matrix3d& m, const double epsstd = 1.0e-3);
// zyx Euler angles and rotation matrix
Eigen::Matrix3d abc2r(const double a, const double b, const double c);
void r2abc(const Eigen::Matrix3d&m, double& a, double& b, double& c);
// quaternion and zyx Euler angles
void q2abc(const Eigen::Quaterniond& q);
Eigen::Quaterniond abc2q(const double a, const double b, const double c);



