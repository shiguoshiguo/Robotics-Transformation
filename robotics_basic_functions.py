# -*- coding: utf-8 -*-
import sys, os
_path = os.path.abspath(__file__)
_path = os.path.dirname(os.path.dirname(_path))
sys.path.append(_path)

import numpy as np
from math import atan2,sin,cos,asin,acos,pi,sqrt

# ######function list
#
# ### basic rotation matrix
# rotx4
# roty4
# rotz4
# rotx3
# roty3
# rotz3
# trans
#
# ### basic rotation transformation functions
# aa2r
# euler2isom
# euler2q
# euler2r
# q2euler 
# q2r
# r2aa
# r2euler
# r2q
#
# ### senior rotation transformation functions
# euler2isomd
# euler2rv
# isom2eulerangles
# isom2eulerposed
# r2rv
# rv2euler
# rv2r
# xyz2abc
# xyz2abcd

# ### math functions
# compute_common
# make_orth
# mpi2pi
# included_angle
# point2plane
# 
# ### kinematics
# joint_mdh
# joint_sdh
# isom2mdh
# isom2sdh


# ### basic rotation matrix
def rotx4(rx):
    '''Jason Stone, April, 2019'''
    rx = float(rx)
    return [[1, 0, 0, 0],
            [0, cos(rx), -sin(rx), 0],
            [0, sin(rx),  cos(rx), 0],
            [0, 0, 0, 1]]


def roty4(ry):
    '''Jason Stone, April, 2019''' 
    ry = float(ry)
    return [[cos(ry), 0, sin(ry), 0],
            [0, 1, 0, 0],
            [-sin(ry), 0, cos(ry), 0],
            [0, 0, 0, 1]]


def rotz4(rz):
    '''Jason Stone, April, 2019''' 
    rz = float(rz)
    return [[cos(rz), -sin(rz), 0, 0],
            [sin(rz),  cos(rz), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]]

def rotx3(rx):
    '''Jason Stone, May, 2019'''
    return np.asarray(rotx4(rx))[0:3, 0:3].tolist()


def roty3(ry):
    '''Jason Stone, May, 2019'''
    return np.asarray(roty4(ry))[0:3, 0:3].tolist()

def rotz3(rz):
    '''Jason Stone, May, 2019'''
    return np.asarray(rotz4(rz))[0:3, 0:3].tolist()

def trans(pos):
    '''Jason Stone, April, 2019''' 
    return [[1, 0, 0, float(pos[0])],
            [0, 1, 0, float(pos[1])],
            [0, 0, 1, float(pos[2])],
            [0, 0, 0, 1]]

# # attitude_sdh
# def attitude_sdh():
#     '''Jason Stone, April, 2019''' 
#     pass

# # attitude_mdh
# def attitude_mdh():
#     '''Jason Stone, April, 2019''' 
#     pass
# #
# ### end of new arrival
#
# ######

def aa2r(angle, axis):
    '''transform angle-axis to rotation matrix, radians
    Jason Stone, April, 2019'''
    if np.linalg.norm(axis) == 0:
        return 
    t = angle
    axis = np.asarray(axis)
    naxis = axis/np.linalg.norm(axis)
    R = np.zeros((3,3))
    nx = naxis[0]
    ny = naxis[1]
    nz = naxis[2]
    R[0,0] = nx * nx * (1 - cos(t)) + cos(t)
    R[0,1] = nx * ny * (1 - cos(t)) - nz * sin(t)
    R[0,2] = nx * nz * (1 - cos(t)) + ny * sin(t)
    R[1,0] = nx * ny * (1 - cos(t)) + nz * sin(t)
    R[1,1] = ny * ny * (1 - cos(t)) + cos(t)
    R[1,2] = ny * nz * (1 - cos(t)) - nx * sin(t)
    R[2,0] = nx * nz * (1 - cos(t)) - ny * sin(t)
    R[2,1] = ny * nz * (1 - cos(t)) + nx * sin(t)
    R[2,2] = nz * nz * (1 - cos(t)) + cos(t)
    return R.tolist()

def euler2isom(pose):
    ''' transform zyx euler angles pose to isometry, 
    position and attitude or only attitude are both supported.
    angles are represented in radians
    Jason Stone, April, 2019'''
    if len(pose) == 6:
        tx = pose[0]
        ty = pose[1]
        tz = pose[2]
        rz = pose[3]
        ry = pose[4]
        rx = pose[5]
    elif len(pose) == 3:
        tx = 0
        ty = 0
        tz = 0
        rz = pose[0]
        ry = pose[1]
        rx = pose[2]
    else:
        return 

    T = [
    [ cos(ry)*cos(rz), 
      cos(rz)*sin(rx)*sin(ry) - cos(rx)*sin(rz), 
      sin(rx)*sin(rz) + cos(rx)*cos(rz)*sin(ry),
      tx],
    [ cos(ry)*sin(rz), 
      cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz), 
      cos(rx)*sin(ry)*sin(rz) - cos(rz)*sin(rx), 
      ty],
    [ -sin(ry), cos(ry)*sin(rx), cos(rx)*cos(ry), tz],
    [ 0, 0, 0, 1]    ]
    return T

def euler2q( zyx ):
    ''' transform zyx euler angles(radians) to quaternion(w,x,y,z), radians
    Jason Stone, April, 2019'''
    [a,b,c] = zyx
    w = cos(a / 2.0) * cos(b / 2.0) * cos(c / 2.0) + sin(a / 2.0) * sin(b / 2.0) * sin(c / 2.0)
    x = sin(a / 2.0) * cos(b / 2.0) * cos(c / 2.0) - cos(a / 2.0) * sin(b / 2.0) * sin(c / 2.0)
    y = cos(a / 2.0) * sin(b / 2.0) * cos(c / 2.0) + sin(a / 2.0) * cos(b / 2.0) * sin(c / 2.0)
    z = cos(a / 2.0) * cos(b / 2.0) * sin(c / 2.0) - sin(a / 2.0) * sin(b / 2.0) * cos(c / 2.0)
    return [w,x,y,z]

def euler2r(angles):
    ''' transform zyx euler angles pose to rotation matrix, 
    angles are represented in radians
    Jason Stone, April, 2019'''
    rz = angles[0]
    ry = angles[1]
    rx = angles[2]

    R = [
        [ cos(ry)*cos(rz), cos(rz)*sin(rx)*sin(ry) - cos(rx)*sin(rz),
          sin(rx)*sin(rz) + cos(rx)*cos(rz)*sin(ry)],
        [ cos(ry)*sin(rz), cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz),
          cos(rx)*sin(ry)*sin(rz) - cos(rz)*sin(rx)],
        [ -sin(ry), cos(ry)*sin(rx), cos(rx)*cos(ry)]]
    return R

def q2euler(q):
    '''transform a qunatenion to zyx euler angles(radians)
    Jason Stone, April, 2019'''
    q = np.asarray(q).flatten()
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    a = atan2(2.0 * (w * x + y * z), 1 - 2.0 * (x * x + y * y))
    b = asin(2.0 * (w * y - x * z))
    c = atan2(2.0 * (w * z + x * y), 1 - 2.0 * (y * y + z * z))
    return [a,b,c]

def q2r(q):
    '''transform a quaternion to  rotation matrix
    Jason Stone, April, 2019'''
    q = np.asarray(q).flatten()
    if np.shape(q) != (4,):
        return
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    m11 = w * w + x * x - y * y - z * z
    m12 = 2.0 * (x * y - w * z)
    m13 = 2.0 * (x * z + w * y)
    m21 = 2.0 * (x * y + w * z)
    m22 = w * w - x * x + y * y - z * z
    m23 = 2.0 * (y * z - w * x)
    m31 = 2.0 * (x * z - w * y)
    m32 = 2.0 * (y * z + w * x)
    m33 = w * w - x * x - y * y + z * z
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def r2aa(R):
    '''transform rotation matrix to angle-axis, radians
    Jason Stone, April, 2019'''
    R = np.asarray(R)
    if np.shape(R)[-2:] != (3,3):
        return
    angle = acos((np.trace(R) - 1)/2.)
    axis = np.array([R[2,1] - R[1,2],R[0,2] - R[2,0], R[1,0] - R[0,1]]) / sin(angle) / 2
    return [angle, axis]

def r2euler(R):
    ''' transform rotation matrix to zyx euler angles, two solutions are provided, radians
    Jason Stone, April, 2019'''
    cysz = R[1][0]
    cycz = R[0][0]
    msy = R[2][0]
    sxcy = R[2][1]
    cxcy = R[2][2]
    a1 = atan2(cysz, cycz)
    a2 = atan2(-cysz, -cycz)	
    b1 = atan2(-msy, sqrt(sxcy ** 2 + cxcy ** 2))
    b2 = atan2(-msy, -sqrt(sxcy ** 2 + cxcy ** 2))
    c1 = atan2(sxcy, cxcy)
    c2 = atan2(-sxcy, -cxcy)
    return [[a1,b1,c1],[a2,b2,c2]]

def r2q(m, epsstd = 1.0e-3):
    '''transform the rotation matrix to a quaternion
    Jason Stone, April, 2019'''
    r11 = m[0, 0], r12 = m[0, 1], r13 = m[0, 2]
    r21 = m[1, 0], r22 = m[1, 1], r23 = m[1, 2]
    r31 = m[2, 0], r32 = m[2, 1], r33 = m[2, 2]
    judge = sqrt(1.0 + r11 + r22 + r33) / 2.0
    if abs(judge) > epsstd:
        w = judge
        x = (r32 - r23) / (4.0 * w)
        y = (r13 - r31) / (4.0 * w)
        z = (r21 - r12) / (4.0 * w)
    else:
        if r11 > r22 and r11 > r33:
            t = 2.0 * sqrt(1.0 + r11 - r22 - r33)
            w = (r32 - r23) / t
            x = t / 4
            y = (r12 + r21) / t
            z = (r13 + r31) / t
        elif (r22 > r33):
            t = 2.0 * sqrt(1.0 - r11 + r22 - r33)
            w = (r13 - r31) / t
            x = (r12 + r21) / t
            y = t / 4
            z = (r23 + r32) / t
        else:
            t = 2.0 * sqrt(1.0 - r11 - r22 + r33)
            w = (r21 - r12) / t
            x = (r13 + r31) / t
            y = (r23 + r32) / t
            z = t / 4
    q = np.array([w,x,y,z])
    return (q / np.linalg.norm(q)).tolist()

###### senior tansformation functions
def euler2isomd(posed):
    ''' transform zyx euler angles pose to isometry, 
    position and attitude or only attitude are both supported.
    angles are represented in degrees
    Jason Stone, April, 2019'''
    pose = np.asarray(posed).flatten()
    if len(pose) == 6:
        pose[3:6] = np.radians(pose[3:6])
    elif len(pose) == 3:
        pose = np.radians(pose) 
    else:
        return 
    return euler2isom(pose)

def xyzeuler2isomd(posed):
    pose = np.asarray(posed).flatten()
    if len(pose) == 6:
        pose[3:6] = np.radians(pose[3:6])
    elif len(pose) == 3:
        pose = np.radians(pose) 
    else:
        return 
    temp = pose[3]
    pose[3] = pose[5]
    pose[5] = temp
    return euler2isom(pose)

def euler2rv(angles):
    '''transform zyx euler angles to rotation vector, radians
    Jason Stone, April, 2019'''
    r = euler2r(angles)
    ag,ax = r2aa(r)
    return (ag * np.asarray(ax)).tolist()

def isom2eulerangles(isom):
    '''compute only zyx euler angles from isometry matrix, radians
    Jason Stone, April, 2019'''
    isom = np.asarray(isom)
    if np.shape(isom) != (4, 4):
        return
    return r2euler(isom[0:3, 0:3])

def isom2euleranglesd(isom):
    '''compute only zyx euler angles from isometry matrix, degrees
    Jason Stone, April, 2019'''
    isom = np.asarray(isom)
    if np.shape(isom) != (4, 4):
        return
    return np.degrees(r2euler(isom[0:3, 0:3])).tolist()

def isom2eulerposed(isometry):
    '''transform isometry to zyx euler angles, degrees
    Jason Stone, April, 2019'''
    isom = np.asarray(isometry)
    if np.shape(isom) != (4,4):
        return 
    R = isom[0:3,0:3]
    angles = np.degrees(r2euler(R))
    eulerposed = []
    pos = isom[0:3,3]
    for i in range(0,np.shape(angles)[0]):
        eulerposed.append(pos.tolist() + angles[i,...].tolist())
    return eulerposed

def r2rv(R):
    '''rotation matrix to rotation vector, radiansJason Stone, April, 2019'''
    [ag,ax] = r2aa(R)
    return  ag * ax

def rv2euler(rv):
    '''transform rotation vector to zyx euler angles, radians
    Jason Stone, April, 2019'''
    nv = np.linalg.norm(rv)
    if len(rv) !=  3 or nv == 0:
        return
    ax = np.asarray(rv) / nv
    return r2euler(aa2r(nv, ax))

def rv2r(rv):
    '''rotation vector to rotatin matrix
    Jason Stone, April, 2019'''
    rv = np.asarray(rv).flatten()
    t = np.linalg.norm(rv)
    if np.shape(rv) != (3,) or t == 0:
        return     
    ax = rv / t
    return aa2r(t,ax)

def xyz2isom(xyz):
    isom = np.mat(rotx4(xyz[0])) * np.mat(roty4(xyz[1])) * np.mat(rotz4(xyz[2]))
    return isom.tolist()

def abc2isom(abc):
    isom = np.mat(rotz4(abc[0])) * np.mat(roty4(abc[1])) * np.mat(rotx4(abc[2]))
    return isom.tolist()

def xyz2isomd(xyz):
    xyz = np.radians(xyz)
    return xyz2isom(xyz)

def abc2isomd(abc):
    abc = np.radians(abc)
    return abc2isom(abc)



def xyz2abc( xyz ):
    """
Function xyz2abcd transform XYZ Euler angles to ZY'X'' Euler angles
    Jason Stone, July, 2019
    """
    isom = np.asarray(np.mat(rotx4(xyz[0])) * np.mat(roty4(xyz[1])) * np.mat(rotz4(xyz[2])))
    return isom2eulerangles(isom)

def xyz2abcd( xyz ):
    """
Function xyz2abcd transform XYZ Euler angles to ZY'X'' Euler angles
    Jason Stone, July, 2019
    """
    xyz = np.radians(np.asarray(xyz).copy())
    
    return np.degrees(xyz2abc(xyz)).tolist()

def xyzxyz2xyzrpy(xyzxyz):
    xyz = xyzxyz[0:3]
    xyzrpy = xyz + xyz2abc(xyzxyz[3:6])[0]
    return xyzrpy

def xyzxyz2xyzrpyd(xyzxyz):
    xyz = xyzxyz[0:3]
    xyzrpy = xyz + xyz2abcd(xyzxyz[3:6])[0]
    return xyzrpy

### math functions
def compute_common(p1,v1,p2,v2):
    '''compute the common of two lines that are presented in points and vectors,
        returns the two terminal points of the common. 
        The first point will stay still if the lines are parrall.
    Jason Stone, April, 2019'''
    p1 = np.asarray(p1)
    v1 = np.asarray(v1)
    p2 = np.asarray(p2)
    v2 = np.asarray(v2)
    vp = p2 - p1
    m = np.asarray([[ -1 * np.dot(v1,v1), np.dot(v1,v2)],[-1 * np.dot(v1,v2), np.dot(v2,v2)]])
    if abs(np.linalg.det(m)) > 1.0e-3:
        b = np.asarray([-1 * np.dot(vp,v1), -1 * np.dot(vp,v2)])
        t = np.matmul(np.linalg.inv(m),b)
    else:
        t = np.asarray([0,np.dot((p1-p2),v2)/np.dot(v2,v2)])
    print(t)
    return p1 + t[0] * v1, p2 + t[1] * v2

def make_orth(v1,v2):
    ''' expand the intersection angle of two vectors into pi/2
        to make them orth, only for testing.
    Jason Stone, April, 2019'''
    # no check
    in_angle = lambda v1,v2 : np.arccos( np.dot(v1,v2) / np.linalg.norm(v1) / np.linalg.norm(v2) )
    v1 = np.asarray(v1).flatten()
    v2 = np.asarray(v2).flatten()
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    v = [v1,v2]
    vc = (v1 + v2) / np.linalg.norm(v1 + v2) 
    vorth = []
    for vi in v:
        axi = np.cross(vc,vi)
        agi = in_angle(vc,vi)
        rai = pi/4 - agi
        ri = np.asarray(aa2r(rai, axi))
        vorth.append(np.matmul(ri,vi).tolist())
    return vorth

def mpi2pi(angles):
    '''limit all angles into range -pi ~ pi'''
    '''Jason Stone, April, 2019'''
    angles = np.asarray(angles).copy()
    angles = np.mod(angles, 2 * pi)
    def pipi(ag):
        if ag >= pi:
            return ag - 2 * pi
        elif ag < -pi:
            return ag + 2 * pi
        else:
            return ag

    _mark = 0
    if len(np.shape(angles)) == 1:
        _mark = 1
        angles = np.expand_dims(angles,axis=0)
    for i in range(0,np.shape(angles)[0]):
        for j in range(0,np.shape(angles)[1]):
            angles[i, j] = pipi(angles[i,j])
    if _mark:
        return angles.flatten().tolist()
    else:
        return angles.tolist()

def included_angle(v1,v2):
    '''compute the inlcuded angle of two 3d vector,radian
    Jason Stone, April, 2019'''
    v1 = np.asarray(v1).flatten()
    v2 = np.asarray(v2).flatten()
    if np.shape(v1) != (3,) or np.shape(v2) != (3,):
        return
    if np.linalg.norm(v1) == 0 or np.linalg.norm(v2) == 0:
        return
    return acos( np.dot(v1,v2)/ np.linalg.norm(v1) / np.linalg.norm(v2) )

def point2plane(p,nv,p0):
    '''project the point p onto a plane(normal vector and p0)
    Jason Stone, April, 2019'''
    p = np.asarray(p).flatten()
    nv = np.asarray(nv).flatten()
    p0 = np.asarray(p0).flatten()
    if not (np.shape(p) == np.shape(nv) == np.shape(p0) == (3,)):
        return
    dp = p0 - p
    return (p + np.dot(dp, nv) / np.dot(nv, nv) * nv).tolist()

# ### kinematics

# joint_mdh
def joint_mdh(dh, theta = 0):
    '''Forward kinematics of single joint based on modified DH model., radians
    Jason Stone, April, 2019'''
    dh = np.array(dh).copy().flatten()
    if np.shape(dh) != (4,) and np.shape(dh) != (5,):
        return
    elif np.shape(dh) == (4,):
        offset = dh[2]
    else:
        offset = dh[4]

    api = dh[0]
    ai = dh[1]
    thi = theta + offset
    di = dh[3]
    return [[cos(thi), -sin(thi), 0, ai],
            [cos(api) * sin(thi), cos(api) * cos(thi), -sin(api), -di * sin(api)],
            [sin(api) * sin(thi), cos(thi) * sin(api), cos(api), di * cos(api)],
            [0, 0, 0, 1]]

def joint_mdhd(dh, theta = 0.0):
    '''Forward kinematics of single joint based on modified DH model., degrees
    Jason Stone, April, 2019'''
    dh = np.array(dh).copy().flatten()
    if np.shape(dh) != (4,) and np.shape(dh) != (5,):
        return
    elif np.shape(dh) == (4,):
        dh[0] = np.radians(dh[0])
        dh[2] = np.radians(dh[2])
    else:
        dh[0] = np.radians(dh[0])
        dh[2] = np.radians(dh[2])
        dh[4] = np.radians(dh[4])
    return joint_mdh(dh, np.radians(theta))

def joint_sdh(dh, theta = 0):
    '''Forward kinematics of single joint based on standard DH model, radians
    Jason Stone, April, 2019'''
    dh = np.array(dh).copy().flatten()
    if np.shape(dh) != (4,) and np.shape(dh) != (5,):
        return
    elif np.shape(dh) == (4,):
        offset = dh[0]
    else:
        offset = dh[4]
    thi = theta + offset
    di = dh[1]
    ai = dh[2]
    api = dh[3]
    return [[cos(thi), -sin(thi) * cos(api), sin(thi) * sin(api), ai * cos(thi)],
            [sin(thi), cos(thi) * cos(api), -cos(thi) * sin(api), ai * sin(thi)],
            [0, sin(api), cos(api), di],
            [0, 0, 0, 1]]

def joint_sdhd(dh, theta = 0.0):
    '''Forward kinematics of single joint based on modified DH model., degrees
    Jason Stone, April, 2019'''
    dh = np.array(dh).copy().flatten()
    if np.shape(dh) != (4,) and np.shape(dh) != (5,):
        return
    elif np.shape(dh) == (4,):
        dh[0] = np.radians(dh[0])
        dh[3] = np.radians(dh[3])
    else:
        dh[0] = np.radians(dh[0])
        dh[3] = np.radians(dh[3])
        dh[4] = np.radians(dh[4])
    return joint_sdh(dh, np.radians(theta))

def isom2sdh(m):
    ''' calculate Standard D-H parameters from transform isometry,
        paras are arranged in theta - d - a - alpha, radians
        Jason Stone, April, 2019'''
    '''
    A{i}=[ cos(th(i)),   -sin(th(i))*cos(ap(i)),...
               sin(th(i))*sin(ap(i)),   a(i)*cos(th(i));           
           sin(th(i)),  cos(th(i))*cos(ap(i)),...
              -cos(th(i))*sin(ap(i)),  a(i)*sin(th(i));           
           0, sin(ap(i)), cos(ap(i)),  d(i);           
           0,0,0,1];
    '''
    if np.shape(m) != (4,4):
        return
    th = atan2(m[1,0],m[0,0])
    d = m[2,3]
    if abs(m[0,0]) > abs(m[1,0]):
        a = m[0,3] / m[0,0]
    else:
        a = m[1,3] / m[1,0]
    ap = atan2(m[2,1],m[2,2])
    return [th, d, a, ap]

def isom2mdh(m):
    ''' calculate Modified D-H parameters from transform isometry,
        paras are arranged in theta - d - a - alpha, radians
        Jason Stone, April, 2019'''
    '''
    A{i} = [    cos(th(i)), -sin(th(i)),   0,  a(i);
                cos(ap(i))*sin(th(i)), cos(ap(i))*cos(th(i)),...
                    -sin(ap(i)), -d(i)*sin(ap(i));
                sin(ap(i))*sin(th(i)), cos(th(i))*sin(ap(i)), ...
                    cos(ap(i)),  d(i)*cos(ap(i));
                0,   0,   0,   1]; 
    '''
    if np.shape(m) != (4,4):
        return
    ap = atan2(-m[1,3],m[2,2])
    a = m[0,3]
    th = atan2(-m[0,1],m[0,0])
    if abs(m[1,2]) > abs(m[2,2]):
        d = m[1,3] / m[1,2]
    else:
        d = m[2,3] / m[2,2]
    return [ap, a, th, d]

### testing
if __name__ == "__main__":
    # import robotics_basic_functions as rbf
    # help(rbf)
    # angles = [pi/4,pi/3,pi/2]
    # rv = euler2rv(angles)
    # print(rv)
    # angles2 = rv2euler(rv)
    # print(angles2)
    # v1 = [0.707,0.707,0]
    # v2 = [0.707,0,0.707]
    # v_orth = make_orth(v1,v2)
    # print(v_orth[0])
    # print(v_orth[1])

    xyzxyz = [0.371017,	0.184388,	0.26273,	-0.000192,	-0.000024,	0.000685]
    print(xyzxyz2xyzrpyd(xyzxyz))

    xyzxyz = [0.000000,	0.000000,	0.581164,	-0.000350,	0.001075,	0.000000]
    print(xyzxyz2xyzrpyd(xyzxyz))