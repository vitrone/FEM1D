#!/usr/bin/env python
import numpy as np
from ctypes import *
import unittest
import testoob
import math

from legendre import legendre
#------------------------------------------------------------------------------#

def LP10(x):
    lp = (1/256)*(46189*x**10-109395*x**8+90090*x**6-30030*x**4+3465*x**2-63)
    return lp

def LP09(x):
    lp = (1/128)*(12155*x**9-25740*x**7+18018*x**5-4620*x**3+315*x)
    return lp
def gauss(x):
    r = np.exp(-16*x**2)
    return r
#------------------------------------------------------------------------------#

class TestSeq(unittest.TestCase):

    def setUp(self):
        self.libc = cdll.LoadLibrary("../BUILD/libfem1d.so")
        self.lp = legendre()

    def test_LGLpoly(self):
        x = 0.6
        n = 2
        y = self.lp.calc_LGLpoly(x,n)
        z = (1-x*x)*3*x
        tol = 1e-9
        self.assertTrue(abs(z-y)/abs(z)<tol)

    def test_Gauss_points(self):
        n   = 10
        tol = 1e-9
        x   = self.lp.find_Gauss_points(n, tol)
        # 10th order Legendre polynomial
        lp = LP10(x)
        for ix in lp[:,0]:
            self.assertTrue(abs(ix)<tol)

    def test_LGL_points(self):
        n   = 10
        tol = 1e-9
        y, quadW = self.lp.find_LGL_points(n, tol)
        # LGL polynomial is defined as (1-x^2)L'_n(x)
        # (1-x^2)L'_n(x) = nL_{n-1}(x)-nxL_n(x)
        lp09 = LP10(y)
        lp10 = LP09(y)
        lgl  = n*lp09 - n*y*lp10
        for ilgl in lgl:
            self.assertTrue(abs(ilgl)<tol)
        poly = y**(2*n+3)
        poly_quad = np.dot(poly.T, quadW) 
        # print abs(poly_quad)
        self.assertTrue(abs(poly_quad)<tol)
    
    def test_transformation_matrix(self):
        p   = 20
        tol = 1e-9
        x, quadW = self.lp.find_LGL_points(p, tol)
        FM = self.lp.forward_transform_matrix( p, x)
        IM = self.lp.backward_transform_matrix( p, p, x)
        Id = np.eye(p+1, dtype=float)
        err = np.max(np.abs(np.dot(IM, FM)-Id), axis=None)
        # print FM[:,0]
        # print err
        self.assertTrue(err<tol)

    def test_transformation_matrix_colmajor(self):
        p   = 20
        tol = 1e-9
        x, quadW = self.lp.find_LGL_points(p, tol)
        FM = self.lp.forward_transform_matrix_colmajor( p, x)
        IM = self.lp.backward_transform_matrix( p, p, x)
        Id = np.eye(p+1, dtype=float)
        err = np.max(np.abs(np.dot(IM, FM)-Id), axis=None)
        # print FM[:,0]
        # print err
        self.assertTrue(err<tol)
        

    def test_interpolation(self):
        p   = 2
        tol = 1e-9
        x, quadW = self.lp.find_LGL_points(p, tol)
        FM = self.lp.forward_transform_matrix( p, x)
        P = 5
        xi = np.zeros([P+1, 1], dtype=float, order='C')
        xi[:,0] = np.linspace(-1, 1, P+1) 
        IM = self.lp.backward_transform_matrix( P, p, xi)

        v = np.exp(-16*x**2)
        v = np.exp(-16*x**2)
        V = np.dot(FM, v)
        print V
        v_interpol = (np.dot(IM, V))
        v_actual = np.exp(-16*xi**2)
        err = np.linalg.norm(v_interpol - v_actual)/np.linalg.norm(v_actual)
        # print err
        self.assertTrue(err<tol)
    
    def test_anti_aliasing_colmajor(self):
        p   = 100
        d   = 5
        tol = 1e-9

        x, quadW = self.lp.find_LGL_points(p, tol)
        FM = self.lp.forward_transform_matrix( p, x)
        u  = gauss(x)
        U  = np.dot(FM, u) # Forward DLT
        Pxi = int(2e2)
        xi = np.zeros([Pxi+1, 1], dtype=float, order='C')
        xi[:,0] = np.linspace(-1, 1, Pxi+1) 
        IM0 = self.lp.backward_transform_matrix( Pxi, p, xi)
        u2  = (gauss(xi))**d

        # No anti-aliasing
        P = p
        x, quadW = self.lp.find_LGL_points(P, tol)
        FM = self.lp.forward_transform_matrix2_colmajor( p, P, x)
        IM = self.lp.backward_transform_matrix( P, p, x)
        u_interpol  = np.dot(IM, U)
        u2_interpol = u_interpol**d
        U2 = np.dot(FM, u2_interpol) # Forward DLT
        u2_tmp = np.dot(IM0, U2)
        err0 = np.linalg.norm(u2-u2_tmp)/np.linalg.norm(u2)
        
        # With anti-aliasing
        P = (np.ceil(d*p*0.5)+1).astype(np.int)
        x, quadW = self.lp.find_LGL_points(P, tol)
        FM = self.lp.forward_transform_matrix2_colmajor( p, P, x)
        IM = self.lp.backward_transform_matrix( P, p, x)
        u_interpol  = (np.dot(IM, U))
        u2_interpol = u_interpol**d
        U2 = np.dot(FM, u2_interpol)
        u2_tmp = np.dot(IM0, U2)
        err = np.linalg.norm(u2 - u2_tmp)/np.linalg.norm(u2)
        # FM_matlab = np.genfromtxt("FM.dat")         
        # err = np.max(np.abs(tmpM), axis=None)
        # print err/err0
        self.assertTrue(err<err0)

    def test_anti_aliasing(self):
        p   = 100
        d   = 5
        tol = 1e-9

        x, quadW = self.lp.find_LGL_points(p, tol)
        FM = self.lp.forward_transform_matrix( p, x)
        u  = gauss(x)
        U  = np.dot(FM, u) # Forward DLT
        Pxi = int(2e2)
        xi = np.zeros([Pxi+1, 1], dtype=float, order='C')
        xi[:,0] = np.linspace(-1, 1, Pxi+1) 
        IM0 = self.lp.backward_transform_matrix( Pxi, p, xi)
        u2  = (gauss(xi))**d

        # No anti-aliasing
        P = p
        x, quadW = self.lp.find_LGL_points(P, tol)
        FM = self.lp.forward_transform_matrix2( p, P, x)
        IM = self.lp.backward_transform_matrix( P, p, x)
        u_interpol  = np.dot(IM, U)
        u2_interpol = u_interpol**d
        U2 = np.dot(FM, u2_interpol) # Forward DLT
        u2_tmp = np.dot(IM0, U2)
        err0 = np.linalg.norm(u2-u2_tmp)/np.linalg.norm(u2)
        
        # With anti-aliasing
        P = (np.ceil(d*p*0.5)+1).astype(np.int)
        x, quadW = self.lp.find_LGL_points(P, tol)
        FM = self.lp.forward_transform_matrix2( p, P, x)
        IM = self.lp.backward_transform_matrix( P, p, x)
        u_interpol  = (np.dot(IM, U))
        u2_interpol = u_interpol**d
        U2 = np.dot(FM, u2_interpol)
        u2_tmp = np.dot(IM0, U2)
        err = np.linalg.norm(u2 - u2_tmp)/np.linalg.norm(u2)
        # FM_matlab = np.genfromtxt("FM.dat")         
        # err = np.max(np.abs(tmpM), axis=None)
        # print err/err0
        self.assertTrue(err<err0)


if __name__ == '__main__':
    # python <file_unittest> <test_class>.<test_case>

    #unittest.main()
    testoob.main()


