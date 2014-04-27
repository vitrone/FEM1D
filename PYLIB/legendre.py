#!/usr/bin/env python
import numpy as np
from ctypes import *
import math

#------------------------------------------------------------------------------#
class legendre(object):

    def __init__(self):
        self.libc = cdll.LoadLibrary("../BUILD/libfem1d.so")

    def calc_LGLpoly(self, x, n):
        self.libc.Legendre_LGL.argtypes = [ c_double, c_int]
        self.libc.Legendre_LGL.restype  = c_double
        ctx = c_double(x)
        ctn = c_int(n)
        y = self.libc.Legendre_LGL(ctx,ctn)
        return y
    
    def find_Gauss_points(self, n, tol):
        self.libc.find_Gauss_points.argtypes = [ c_int, c_double, POINTER(c_double)]
        self.libc.find_Gauss_points.restype  = None 
        ctn   = c_int(n)
        cttol = c_double(tol)
        # note that ctx here is a numpy.array
        ctx = np.zeros([n, 1], dtype=float, order='C')
        self.libc.find_Gauss_points( ctn, cttol, 
                                     ctx.ctypes.data_as(POINTER(c_double)))
        return ctx

    def find_LGL_points(self, n, tol ):
        gauss_zeros = self.find_Gauss_points(n, tol)
        self.libc.find_LGL_points.argtypes = [ c_int, c_double, 
                POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.libc.find_LGL_points.restype  = None 
        ctn   = c_int(n)
        xi    = np.zeros([n+1, 1], dtype=float, order='C')
        quadW = np.zeros([n+1, 1], dtype=float, order='C')
        self.libc.find_LGL_points( n, tol, 
                                   xi.ctypes.data_as(POINTER(c_double)),
                                   quadW.ctypes.data_as(POINTER(c_double)), 
                                   gauss_zeros.ctypes.data_as(POINTER(c_double)))
        return xi, quadW

    def forward_transform_matrix(self, p, LGL_points ):

        self.libc.forward_transform_matrix.argtypes = [ c_int, 
                                                        POINTER(c_double),
                                                        POINTER(c_double)]
        self.libc.forward_transform_matrix.restypes = None
        ctp = c_int(p)
        FM = np.zeros([p+1, p+1], dtype=float, order='C')
        self.libc.forward_transform_matrix( ctp,
                                            LGL_points.ctypes.data_as(POINTER(c_double)),
                                            FM.ctypes.data_as(POINTER(c_double)))
        return FM
    def forward_transform_matrix2(self, p, P, LGL_points ):

        self.libc.forward_transform_matrix2.argtypes = [ c_int, c_int,
                                                        POINTER(c_double),
                                                        POINTER(c_double)]
        self.libc.forward_transform_matrix2.restypes = None
        ctp = c_int(p)
        ctP = c_int(P)
        FM = np.zeros([p+1, P+1], dtype=float, order='C')
        self.libc.forward_transform_matrix2( ctp, ctP,
                                            LGL_points.ctypes.data_as(POINTER(c_double)),
                                            FM.ctypes.data_as(POINTER(c_double)))
        return FM
    def backward_transform_matrix(self, P, p, x ):

        ctP = c_int(P)
        ctp = c_int(p)

        self.libc.backward_transform_matrix.argtypes = [ c_int, 
                                                         c_int, 
                                                         POINTER(c_double), 
                                                         POINTER(c_double)]
        self.libc.backward_transform_matrix.restypes = None
        IM = np.zeros([P+1, p+1], dtype=float, order='C')
        self.libc.backward_transform_matrix( ctP, ctp, 
                                             x.ctypes.data_as(POINTER(c_double)),
                                             IM.ctypes.data_as(POINTER(c_double)))
        return IM
    def forward_transform_matrix_colmajor(self, p, LGL_points ):

        self.libc.forward_transform_matrix_colmajor.argtypes = [ c_int, 
                                                        POINTER(c_double),
                                                        POINTER(c_double)]
        self.libc.forward_transform_matrix_colmajor.restypes = None
        ctp = c_int(p)
        FM = np.zeros([p+1, p+1], dtype=float, order='F')
        self.libc.forward_transform_matrix_colmajor( ctp,
                                            LGL_points.ctypes.data_as(POINTER(c_double)),
                                            FM.ctypes.data_as(POINTER(c_double)))
        return FM

    def forward_transform_matrix2_colmajor(self, p, P, LGL_points ):

        self.libc.forward_transform_matrix2_colmajor.argtypes = [ c_int, c_int,
                                                        POINTER(c_double),
                                                        POINTER(c_double)]
        self.libc.forward_transform_matrix2_colmajor.restypes = None
        ctp = c_int(p)
        ctP = c_int(P)
        FM = np.zeros([p+1, P+1], dtype=float, order='F')
        self.libc.forward_transform_matrix2_colmajor( ctp, ctP,
                                            LGL_points.ctypes.data_as(POINTER(c_double)),
                                            FM.ctypes.data_as(POINTER(c_double)))
        return FM
