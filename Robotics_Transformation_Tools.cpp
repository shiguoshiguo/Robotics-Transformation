#include "Robotics_Transformation_Tools.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;


///
/// quaternion based leica tracker transformation
///

LeicaPose::LeicaPose(void){
	
	position = Eigen::Vector3d(0.0, 0.0, 0.0);
	gesture = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);
}


LeicaPose::LeicaPose(double h, double v, double d, const Eigen::Quaterniond& q)
{
	// the Leica camera coordination is left-handed, 
	// so the x,y functions are different from standard ones
	// this means H is the angle by y+ axis
	double px = d * sin(v * M_PI / 180.0) * sin(h * M_PI / 180.0);
	double py = d * sin(v * M_PI / 180.0) * cos(h * M_PI / 180.0);
	double pz = d * cos(v * M_PI / 180.0);
	position = Eigen::Vector3d(px, py, pz);
	gesture = q.normalized();
}

LeicaPose::LeicaPose(double h, double v, double d, double w, double x, double y, double z)
{
	*this = LeicaPose(h, v, d, Eigen::Quaterniond(w, x, y, z));
}

double LeicaPose::H()  const { return atan2(position[0], position[1])  * 180.0 / M_PI; }
double LeicaPose::V()  const { return acos(position[2] / this->D())  * 180.0 / M_PI; }
double LeicaPose::D()  const { return position.norm(); }
double LeicaPose::q0()  const { return gesture.w(); }
double LeicaPose::q1()  const { return gesture.x(); }
double LeicaPose::q2()  const { return gesture.y(); }
double LeicaPose::q3()  const { return gesture.z(); }

///
/// eulaer based kuka robot trasnformation
///

KUKAPose::KUKAPose()
{
	position = Eigen::Vector3d(0.0, 0.0, 0.0);
	gesture = Eigen::Vector3d(0.0, 0.0, 0.0);
}

KUKAPose::KUKAPose(const Eigen::Vector3d& pvec, const Eigen::Vector3d& gvec)
{
	position = pvec;
	gesture = gvec * M_PI / 180.0;
}


KUKAPose::KUKAPose(double x, double y, double z, double a, double b, double c) :
	EulerPose(x, y, z, a * M_PI / 180.0, b * M_PI / 180.0, c * M_PI / 180.0) {}

KUKAPose::KUKAPose(double n1, double n2, double n3, int tag): EulerPose(n1, n2, n3, tag)
{
	double amp;
	if (tag)
	{
		amp = M_PI / 180.0;
	}
	else
	{
		amp = 1.0;
	}
	gesture *= amp;
}

KUKAPose::KUKAPose(const Eigen::Vector3d& v, int tag) : EulerPose(v,tag)
{
	double amp;
	if (tag)
	{
		amp = M_PI / 180.0;
	}
	else
	{
		amp = 1.0;
	}
	gesture *= amp;
}

double KUKAPose::x()  const { return position.x(); }
double KUKAPose::y()  const { return position.y(); }
double KUKAPose::z()  const { return position.z(); }
double KUKAPose::a()  const { return gesture.x() * 180.0 /M_PI; }
double KUKAPose::b()  const { return gesture.y() * 180.0 / M_PI; }
double KUKAPose::c()  const { return gesture.z() * 180.0 / M_PI; }

std::ostream& operator<< (std::ostream& ostr, const KUKAPose kp)
{
	ostr << " " << kp.x() << " " << kp.y() << " " << kp.z() << ", ";
	ostr << kp.a() << " " << kp.b() << " " << kp.c() << " ";
	return ostr;
}

/////////////////////////////////

///
/// Quaternion represented rotation based pose transformation class
///
QuaternionPose::QuaternionPose()
{
	position = Eigen::Vector3d::Zero();
	gesture = Eigen::Quaterniond::Identity();
}

QuaternionPose::QuaternionPose(const Eigen::Vector3d& vec, const Eigen::Quaterniond& q)
{
	position = vec;
	gesture = q.normalized();
}

QuaternionPose::QuaternionPose(double px, double py, double pz, double qw, double qx, double qy, double qz)
{
	*this = QuaternionPose(Eigen::Vector3d(px, py, pz), Eigen::Quaterniond(qw, qx, qy, qz));
}

QuaternionPose::QuaternionPose(double px, double py, double pz)
{
	*this = QuaternionPose(Eigen::Vector3d(px, py, pz), Eigen::Quaterniond::Identity());
}

QuaternionPose::QuaternionPose(double qw, double qx, double qy, double qz)
{
	*this = QuaternionPose(Eigen::Vector3d::Zero(), Eigen::Quaterniond(qw, qx, qy, qz));
}

QuaternionPose::QuaternionPose(const Eigen::Vector3d& pos)
{
	*this = QuaternionPose(pos, Eigen::Quaterniond::Identity());
}

QuaternionPose::QuaternionPose(const Eigen::Quaterniond& q)
{
	*this = QuaternionPose(Eigen::Vector3d::Zero(), q);
}

QuaternionPose::QuaternionPose(const Eigen::Vector3d& pos, const Eigen::Matrix3d& m)
{
	double w, x, y, z;
	double r11 = m(0, 0), r12 = m(0, 1), r13 = m(0, 2);
	double r21 = m(1, 0), r22 = m(1, 1), r23 = m(1, 2);
	double r31 = m(2, 0), r32 = m(2, 1), r33 = m(2, 2);
	double judge = sqrt(1.0 + r11 + r22 + r33) / 2.0;
	if (abs(judge) > 1.0e-3) {
		w = judge;
		x = (r32 - r23) / (4.0 * w);
		y = (r13 - r31) / (4.0 * w);
		z = (r21 - r12) / (4.0 * w);
	}
	else
	{
		int _c = 0;
		if (r11 > r22 && r11 > r33) {
			_c = 1;
		}
		else if (r22 > r33)
		{
			_c = 2;
		}
		else
		{
			_c = 3;
		}
		double t;
		switch (_c)
		{
		case 1:
			t = 2.0 * sqrt(1.0 + r11 - r22 - r33);
			w = (r32 - r23) / t;
			x = t / 4;
			y = (r12 + r21) / t;
			z = (r13 + r31) / t;
			break;
		case 2:
			t = 2.0 * sqrt(1.0 - r11 + r22 - r33);
			w = (r13 - r31) / t;
			x = (r12 + r21) / t;
			y = t / 4;
			z = (r23 + r32) / t;
			break;
		case 3:
			t = 2.0 * sqrt(1.0 - r11 - r22 + r33);
			w = (r21 - r12) / t;
			x = (r13 + r31) / t;
			y = (r23 + r32) / t;
			z = t / 4;
			break;
		default:
			w = 1.0; x = 0.0; y = 0.0; z = 0.0;
			break;
		}
	} // end of judge
	Eigen::Quaterniond q(w, x, y, z);
	*this = QuaternionPose(pos, q);
}

QuaternionPose::QuaternionPose(const Eigen::Matrix3d& m)
{
	*this = QuaternionPose(Eigen::Vector3d::Zero(), m);
}

QuaternionPose::QuaternionPose(const Eigen::Isometry3d& isom)
{
	*this = QuaternionPose(isom.translation(), isom.rotation());
}

QuaternionPose::QuaternionPose(const Eigen::Vector3d& pos, const Eigen::Vector3d& ges)
{
	double a = ges.x(), b = ges.y(), c = ges.z();
	double w = cos(c / 2.0) * cos(b / 2.0) * cos(a / 2.0) + sin(c / 2.0) * sin(b / 2.0) * sin(a / 2.0);
	double x = sin(c / 2.0) * cos(b / 2.0) * cos(a / 2.0) - cos(c / 2.0) * sin(b / 2.0) * sin(a / 2.0);
	double y = cos(c / 2.0) * sin(b / 2.0) * cos(a / 2.0) + sin(c / 2.0) * cos(b / 2.0) * sin(a / 2.0);
	double z = cos(c / 2.0) * cos(b / 2.0) * sin(a / 2.0) - sin(c / 2.0) * sin(b / 2.0) * cos(a / 2.0);
	Eigen::Quaterniond q(w, x, y, z);
	*this = QuaternionPose(pos, q);
}

QuaternionPose::QuaternionPose(const EulerPose& ep)
{
	*this = QuaternionPose(ep.position, ep.gesture);
}

Eigen::Isometry3d QuaternionPose::IsometryMatrix(void) const
{
	Eigen::Isometry3d isot(gesture);
	isot.pretranslate(position);
	return isot;
}


Eigen::Vector3d QuaternionPose::operator*(const Eigen::Vector3d& vec)
{
	return this->gesture * vec + this->position;
}

std::ostream& operator<< (std::ostream& ostr, const QuaternionPose qp)
{
	//"attention : the quaternion is defaulty printed in x-y-z-w sequence ..."
	// we change the sequence by hand
	
	string s1 = "+", s2 = "+", s3 = "+";
	if (qp.gesture.coeffs().x() < 0.0) s1 = "-";
	if (qp.gesture.coeffs().y() < 0.0) s2 = "-";
	if (qp.gesture.coeffs().z() < 0.0) s3 = "-";
	double ax = abs(qp.gesture.coeffs().x());
	double ay = abs(qp.gesture.coeffs().y());
	double az = abs(qp.gesture.coeffs().z());
	ostr << " " << qp.position.transpose() << ", ";
	ostr << qp.gesture.coeffs().w() << " ";
	ostr << s1 << ax << "i ";
	ostr << s2 << ay << "j ";
	ostr << s3 << az << "k";
	return ostr;
}

///
/// Euler angles represented rotation based pose transformation class
///

EulerPose::EulerPose()
{
	position = Eigen::Vector3d::Zero();
	gesture = Eigen::Vector3d::Zero();
}

EulerPose::EulerPose(const Eigen::Vector3d& pvec, const Eigen::Vector3d& gvec)
{
	position = pvec;
	gesture = gvec;
}


EulerPose::EulerPose(double x, double y, double z, double a, double b, double c)
{
	*this = EulerPose(Eigen::Vector3d(x, y, z), Eigen::Vector3d(a, b, c));
}

EulerPose::EulerPose(const Eigen::Vector3d& v, int tag)
{
	// tag == 0, pose; tag == 1 or other values, gesture
	Eigen::Vector3d vp, vg;
	if (tag) {
		vp = Eigen::Vector3d::Zero();
		vg = v;
	}
	else
	{
		vp = v;
		vg = Eigen::Vector3d::Zero();
	}
	*this = EulerPose(vp, vg);
}

EulerPose::EulerPose(double n1, double n2, double n3, int tag)
{
	// tag == 0, pose; tag == 1 or other values, gesture
	Eigen::Vector3d vp, vg;
	if (tag)
	{
		vp = Eigen::Vector3d::Zero();
		vg = Eigen::Vector3d(n1, n2, n3);
	}
	else
	{
		vp = Eigen::Vector3d(n1, n2, n3);
		vg = Eigen::Vector3d::Zero();
	}
	*this = EulerPose(vp, vg);
}


EulerPose::EulerPose(const Eigen::Vector3d& pos, const Eigen::Matrix3d& m)
{
	double m10 = m(1, 0), m00 = m(0, 0);
	double m20 = m(2, 0), m21 = m(2, 1), m22 = m(2, 2);
	double a = atan2(m(1, 0), m(0, 0));
	double b1 = atan2(-m(2, 0), sqrt(pow(m(2, 1), 2) + pow(m(2, 2), 2)));
	double b2 = atan2(-m(2, 0), -sqrt(pow(m(2, 1), 2) + pow(m(2, 2), 2)));
	double b = b1;
	double c = atan2(m(2, 1), m(2, 2));
	*this = EulerPose(pos.x(), pos.y(), pos.z(), a, b, c);
}

EulerPose::EulerPose(const Eigen::Matrix3d& m)
{
	*this = EulerPose(Eigen::Vector3d::Zero(), m);
}

EulerPose::EulerPose(const Eigen::Isometry3d& isom)
{
	*this = EulerPose(isom.translation(), isom.rotation());
}


EulerPose::EulerPose(const Eigen::Vector3d& pos, const Eigen::Quaterniond& q)
{
	double a = atan2(2.0 * (q.w() * q.x() + q.y() * q.z()), 1 - 2.0 * (q.x() * q.x() + q.y() * q.y()));
	double b = asin(2.0 * (q.w() * q.y() - q.x() * q.z()));
	double c = atan2(2.0 * (q.w() * q.z() + q.x() * q.y()), 1 - 2.0 * (q.y() * q.y() + q.z() * q.z()));
	gesture = Eigen::Vector3d(a, b, c);

	*this = EulerPose(pos, q.toRotationMatrix());
}

EulerPose::EulerPose(const QuaternionPose& qp)
{
	*this = EulerPose(qp.position, qp.gesture);
}

Eigen::Isometry3d EulerPose::IsometryMatrix(void) const
{
	double rz = gesture[0], ry = gesture[1], rx = gesture[2];
	Eigen::Matrix3d m;
	m << cos(ry)*cos(rz), cos(rz)*sin(rx)*sin(ry) - cos(rx)*sin(rz), sin(rx)*sin(rz) + cos(rx)*cos(rz)*sin(ry),
		cos(ry)*sin(rz), cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz), cos(rx)*sin(ry)*sin(rz) - cos(rz)*sin(rx),
		-sin(ry), cos(ry)*sin(rx), cos(rx)*cos(ry);	
	Eigen::Isometry3d isot(m);
	isot.pretranslate(position);
	return isot;
}


Eigen::Vector3d EulerPose::operator*(const Eigen::Vector3d& vec)
{
	return this->IsometryMatrix() * vec;
}

std::ostream& operator<< (std::ostream& ostr, const EulerPose ep)
{
	ostr << " " << ep.position.x() << " " << ep.position.y() << " " << ep.position.z() << ", ";
	ostr << ep.gesture[0] << " " << ep.gesture[1] << " " << ep.gesture[2] << " ";
	return ostr;
}


////////////////////
///
/// referenced basic functions
///
void xyz2lhvd(const double x, const double y, const double z, double& H, double & V, double& D)
{	
	D = sqrt(x *x + y * y + z * z);
	if (D == 0.0) return;
	V = acos(z / D);
	H = atan2(x, y);
}

void lhvd2xyz(const double H, const double V, const double D, double& x, double & y, double& z)
{
	x = D * sin(V) * sin(H);
	y = D * sin(V) * cos(H);
	z = D * cos(V);
}

Eigen::Matrix3d q2r(const Eigen::Quaterniond& q)
{
	Eigen::Matrix3d m;
	double m11 = q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z() * q.z();
	double m12 = 2.0 * (q.x() * q.y() - q.w() * q.z());
	double m13 = 2.0 * (q.x() * q.z() + q.w() * q.y());
	double m21 = 2.0 * (q.x() * q.y() + q.w() * q.z());
	double m22 = q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z() * q.z();
	double m23 = 2.0 * (q.y() * q.z() - q.w() * q.x());
	double m31 = 2.0 * (q.x() * q.z() - q.w() * q.y());
	double m32 = 2.0 * (q.y() * q.z() + q.w() * q.x());
	double m33 = q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z() * q.z();

	m << m11, m12, m13, m21, m22, m23, m31, m32, m33;
	return m;
}

Eigen::Quaterniond r2q(const Eigen::Matrix3d& m,const double epsstd)
{
	// haven't correct the else part
	double w, x, y, z;
	double r11 = m(0, 0), r12 = m(0, 1), r13 = m(0, 2);
	double r21 = m(1, 0), r22 = m(1, 1), r23 = m(1, 2);
	double r31 = m(2, 0), r32 = m(2, 1), r33 = m(2, 2);
	double judge = sqrt(1.0 + r11 + r22 + r33) / 2.0;
	if ( abs(judge) > epsstd) {
		w = judge;
		x = (r32 - r23) / (4.0 * w);
		y = (r13 - r31) / (4.0 * w);
		z = (r21 - r12) / (4.0 * w);
	}
	else
	{
		int _c = 0;
		if (r11 > r22 && r11 > r33) {
			_c = 1;
		}
		else if (r22 > r33)
		{
			_c = 2;
		}
		else
		{
			_c = 3;
		}
		double t;
		switch (_c)
		{
		case 1:
			t = 2.0 * sqrt(1.0 + r11 - r22 - r33);
			w = (r32 - r23) / t;
			x = t / 4;
			y = (r12 + r21) / t;
			z = (r13 + r31) / t;
			break;
		case 2:
			t = 2.0 * sqrt(1.0 - r11 + r22 - r33);
			w = (r13 - r31) / t;
			x = (r12 + r21) / t; 
			y = t / 4;
			z = (r23 + r32) / t;
			break;
		case 3:
			t = 2.0 * sqrt(1.0 - r11 - r22 + r33);
			w = (r21 - r12) / t;
			x = (r13 + r31) / t;
			y = (r23 + r32) / t;
			z = t / 4;
			break;
		default:
			return Eigen::Quaterniond::Identity();
			break;
		}
	} // end of judge
	return Eigen::Quaterniond(w, x, y, z).normalized();
}

Eigen::Matrix3d abc2r(const double rz, const double ry, const double rx)
{
	Eigen::Matrix3d m;
	m << cos(ry)*cos(rz), cos(rz)*sin(rx)*sin(ry) - cos(rx)*sin(rz), sin(rx)*sin(rz) + cos(rx)*cos(rz)*sin(ry),
		cos(ry)*sin(rz), cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz), cos(rx)*sin(ry)*sin(rz) - cos(rz)*sin(rx),
		-sin(ry), cos(ry)*sin(rx), cos(rx)*cos(ry);
	return m;

	//// construct by sequenced angle-axis ratation
	//Eigen::Matrix3d rot_mz(Eigen::AngleAxisd(rz, Eigen::Vector3d(0.0, 0.0, 1.0)));
	//Eigen::Matrix3d rot_my(Eigen::AngleAxisd(ry, Eigen::Vector3d(0.0, 1.0, 0.0)));
	//Eigen::Matrix3d rot_mx(Eigen::AngleAxisd(rx, Eigen::Vector3d(1.0, 0.0, 0.0)));
	//Eigen::Matrix3d rot_m = rot_mz * rot_my * rot_mx;

}

void r2abc(const Eigen::Matrix3d&m, double& a, double& b, double& c)
{
	double m10 = m(1, 0), m00 = m(0, 0);
	double m20 = m(2, 0), m21 = m(2, 1), m22 = m(2, 2);
	double _a1 = atan2(m(1, 0), m(0, 0));
	double _a2 = atan2(-m(1, 0), -m(0, 0));
	double _b1 = atan2(-m(2, 0), sqrt(pow(m(2, 1), 2) + pow(m(2, 2), 2)));
	double _b2 = atan2(-m(2, 0), -sqrt(pow(m(2, 1), 2) + pow(m(2, 2), 2)));
	double _c1 = atan2(m(2, 1), m(2, 2));
	double _c2 = atan2(-m(2, 1), -m(2, 2));
	a = _a1;
	b = _b1;
	c = _c1;
}


void q2abc(const Eigen::Quaterniond& q)
{
	double a = atan2(2.0 * (q.w() * q.x() + q.y() * q.z()), 1 - 2.0 * (q.x() * q.x() + q.y() * q.y()));
	double b = asin(2.0 * (q.w() * q.y() - q.x() * q.z()));
	double c = atan2(2.0 * (q.w() * q.z() + q.x() * q.y()), 1 - 2.0 * (q.y() * q.y() + q.z() * q.z()));
}
Eigen::Quaterniond abc2q(const double a, const double b, const double c)
{
	double w = cos(c / 2.0) * cos(b / 2.0) * cos(a / 2.0) + sin(c / 2.0) * sin(b / 2.0) * sin(a / 2.0);
	double x = sin(c / 2.0) * cos(b / 2.0) * cos(a / 2.0) - cos(c / 2.0) * sin(b / 2.0) * sin(a / 2.0);
	double y = cos(c / 2.0) * sin(b / 2.0) * cos(a / 2.0) + sin(c / 2.0) * cos(b / 2.0) * sin(a / 2.0);
	double z = cos(c / 2.0) * cos(b / 2.0) * sin(a / 2.0) - sin(c / 2.0) * sin(b / 2.0) * cos(a / 2.0);
	return Eigen::Quaterniond(w, x, y, z).normalized();
	
	//// construct by concept
	//Eigen::Quaterniond qrz(Eigen::AngleAxisd(a, Eigen::Vector3d(0.0, 0.0, 1.0)));
	//Eigen::Quaterniond qry(Eigen::AngleAxisd(b, Eigen::Vector3d(0.0, 1.0, 0.0)));
	//Eigen::Quaterniond qrx(Eigen::AngleAxisd(c, Eigen::Vector3d(1.0, 0.0, 0.0)));
	//Eigen::Quaterniond q = qrz.normalized() * qry.normalized() * qrx.normalized();
	//q.normalize();
}



///


