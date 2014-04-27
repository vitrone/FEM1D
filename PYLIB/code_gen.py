#!/usr/bin/env python
import math


#------------------------------------------------------------------------------#

def fem1d_dlp_snorm2_d_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);
        for k in range(0, N):
            num_str = format((2/(2*k+1.0)), num_format )
            mystr = "#define _N%i %s\n" %( k, num_str )
            fobj.write(mystr);
        for p in range(2, N):
            mystr = "double fem1d_dlp_snorm2_d_%i( Index elem_n, double *u);\n" % p
            fobj.write(mystr);

        for p in range( 2, N):
    
            mystr = "/*** Begin Subroutine ***/\n"
            fobj.write(mystr);
            mystr = "double fem1d_dlp_snorm2_d_%i(Index elem_n, double *u)\n{\n" % p
            fobj.write(mystr);
            mystr = "%sdouble snorm = 0, tmp;\n" % ST
            fobj.write(mystr);
            mystr = "%sIndex i;\n" % ST
            fobj.write(mystr);
            mystr = "%sIndex tmp1 = 0;\n" % ST
            fobj.write(mystr);
            mystr = "%sfor (i=0; i<elem_n; i++)\n%s{\n" % (ST, ST)
            fobj.write(mystr);
            mystr = "%stmp = \t  _N0 * *(u+0) * *(u+0)\n" % ST2 
            fobj.write(mystr);
            for k in range(1,p):
                mystr = "%s\t+ _N%i * *(u+%i) * *(u+%i)\n" % (ST2, k, k, k)
                fobj.write(mystr);
            mystr = "%s\t+ _N%i * *(u+%i) * *(u+%i);\n" % (ST2, p, p, p)
            fobj.write(mystr);
            mystr = "%su += %i;\n" % (ST2, p+1) 
            fobj.write(mystr);
            mystr = "%ssnorm += tmp;\n" % ST2 
            fobj.write(mystr);
            mystr = "%s}\n" % ST
            fobj.write(mystr);
            mystr = "%sreturn snorm;\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return
def fem1d_zlp_snorm2_d_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);
        for k in range(0, N):
            num_str = format((2/(2*k+1.0)), num_format )
            mystr = "#define _N%i %s\n" %( k, num_str )
            fobj.write(mystr);
        for p in range(2, N):
            mystr = "double fem1d_zlp_snorm2_d_%i( Index elem_n, MATLIB_COMPLEX *u);\n" % p
            fobj.write(mystr);

        for p in range( 2, N):
    
            mystr = "/*** Begin Subroutine ***/\n"
            fobj.write(mystr);
            mystr = "double fem1d_zlp_snorm2_d_%i(Index elem_n, MATLIB_COMPLEX *u)\n{\n" % p
            fobj.write(mystr);
            mystr = "%sdouble snorm = 0, tmp;\n" % ST
            fobj.write(mystr);
            mystr = "%sIndex i;\n" % ST
            fobj.write(mystr);
            mystr = "%sIndex tmp1 = 0;\n" % ST
            fobj.write(mystr);
            mystr = "%sfor (i=0; i<elem_n; i++)\n%s{\n" % (ST, ST)
            fobj.write(mystr);
            mystr = "%stmp = \t  _N0 * *(u+0) * conj(*(u+0))\n" % ST2 
            fobj.write(mystr);
            for k in range(1,p):
                mystr = "%s\t+ _N%i * *(u+%i) * conj(*(u+%i))\n" % (ST2, k, k, k)
                fobj.write(mystr);
            mystr = "%s\t+ _N%i * *(u+%i) * conj(*(u+%i));\n" % (ST2, p, p, p)
            fobj.write(mystr);
            mystr = "%su += %i;\n" % (ST2, p+1) 
            fobj.write(mystr);
            mystr = "%ssnorm += tmp;\n" % ST2 
            fobj.write(mystr);
            mystr = "%s}\n" % ST
            fobj.write(mystr);
            mystr = "%sreturn snorm;\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return
#------------------------------------------------------------------------------#

def fem_ddot_(N, file_name):

    ST  = "    "
    with open ( file_name, 'w') as fobj:
        for p in range(2, N):
            mystr = "double fem1d_ddot%i( Index l, double *u, double *v);\n" % p
            fobj.write(mystr);
        for p in range( 2, N):
    
            mystr = "/*** Begin Subroutine ***/\n"
            fobj.write(mystr);
            mystr = "double fem1d_ddot%i(Index l, double *u, double *v)\n{\n" % p
            fobj.write(mystr);
            mystr = "%sdouble tmp;\n" % ST
            fobj.write(mystr);
            mystr = "%stmp =\t  u[0] * v[0]\n" % ST 
            fobj.write(mystr);
            for k in range(1,p-1):
                mystr = "%s%s\t+ u[%i] * v[%i]\n" % (ST,ST,k,k)
                fobj.write(mystr);
            mystr = "%s%s\t+ u[%i] * v[%i];\n" % (ST,ST,p-1,p-1)
            fobj.write(mystr);
            mystr = "%sreturn tmp;\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return

def fem_zdot_(N, file_name):

    ST  = "    "
    with open ( file_name, 'w') as fobj:
        for p in range(2, N):
            mystr = "MATLIB_COMPLEX fem1d_zdot%i( Index l, MATLIB_COMPLEX *u, MATLIB_COMPLEX *v);\n" % p
            fobj.write(mystr);
        for p in range( 2, N):
    
            mystr = "/*** Begin Subroutine ***/\n"
            fobj.write(mystr);
            mystr = "MATLIB_COMPLEX fem_zdot%i(Index l, MATLIB_COMPLEX *u, MATLIB_COMPLEX *v)\n{\n" % p
            fobj.write(mystr);
            mystr = "%sMATLIB_COMPLEX tmp;\n" % ST
            fobj.write(mystr);
            mystr = "%stmp =\t  u[0] * v[0]\n" % ST 
            fobj.write(mystr);
            for k in range(1,p-1):
                mystr = "%s%s\t+ u[%i] * v[%i]\n" % (ST,ST,k,k)
                fobj.write(mystr);
            mystr = "%s%s\t+ u[%i] * v[%i];\n" % (ST,ST,p-1,p-1)
            fobj.write(mystr);
            mystr = "%sreturn tmp;\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return
#------------------------------------------------------------------------------#

def fem_shapeFunc2lp_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);

        # Printing the coefficients
        s = 0;
        s = s + 1
        num_str = format(-1.0/math.sqrt(6), num_format )
        coeff   = '_A{0:02d}'.format(s)
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        s = s + 1
        num_str = format(-1.0/math.sqrt(10), num_format )
        coeff   = '_A{0:02d}'.format(s)
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        for k in range(0, N-3):
            s = s + 1
            num_str = format(1.0/math.sqrt(2*(2*k+3)), num_format )
            coeff   = '_A{0:02d}'.format(s)
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
            num_str = format(-1.0/math.sqrt(2*(2*k+7)), num_format )
            coeff   = '_B{0:02d}'.format(s)
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
        s = s + 1
        num_str = format(1.0/math.sqrt(2*(2*N-3)), num_format )
        coeff   = '_A{0:02d}'.format(s)
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        s = s + 1
        num_str = format(1.0/math.sqrt(2*(2*N-1)), num_format )
        coeff   = '_A{0:02d}'.format(s)
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);

        # function enumeration for the header file
        for p in range( 1, N+1 ):
            mystr = "void fem_shapefunc2lp_%i( Index elem_n, double *v, double *b, double *u );\n"% p
            fobj.write(mystr);

        for p in range( 1, N+1 ):
            mystr = "\n/*** Begin Subroutine ***/\n"
            fobj.write(mystr);

            mystr = """
void fem_shapefunc2lp_%i
(
    Index elem_n, 
    double *v,
    double *b,
    double *u
)
{
"""% p
            fobj.write(mystr);
            mystr = """
    Index i;
    for (i=0; i<elem_n; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;\n"""
            fobj.write(mystr);
            s = 2
            for k in range(0,p-3):
                s = s + 1
                coeff = '_A{0:02d}'.format(s)
                mystr = "%s*u = %s**b + " %( ST2, coeff ) 
                coeff = '_B{0:02d}'.format(s)
                mystr = "%s%s**(b+2), u++, b++;\n" %( mystr, coeff )
                fobj.write(mystr);
            coeff1 = '_A{0:02d}'.format(s+1)
            coeff2 = '_A{0:02d}'.format(s+2)
            mystr = """
        *u = %s**b, u++, b++;
        *u = %s**b, u++, b++;\n""" %( coeff1, coeff2)
            fobj.write(mystr);
            mystr = "%s}\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return

#------------------------------------------------------------------------------#

def lp2fem1d_dshapefunc_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);

        # Printing the coefficients
        s = N;
        t = N - 2
        for k in range(0, N-1):
            coeff = '_C{0:02d}'.format(s)
            num_str = format(math.sqrt(2*(2*s-1)), num_format )
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
            s = s - 1
        for k in range(0, N-3):
            coeff = '_D{0:02d}'.format(t)
            num_str = format(math.sqrt(2*(2*t-1))/math.sqrt(2*(2*t+3)), num_format )
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
            t = t - 1
        coeff = '_D{0:02d}'.format(t)
        num_str = format(1/math.sqrt(2*(2*t+3)), num_format )
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        t = t - 1
        coeff = '_D{0:02d}'.format(t)
        num_str = format(1/math.sqrt(2*(2*t+3)), num_format )
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        # function enumeration for the header file
        for p in range( 2, N+1 ):
            mystr = \
"void lp2fem1d_dshapefunc_%i( Index elem_n, double *u, double *v, double *b );\n" %p
            fobj.write(mystr);

        for p in range( 2, N+1 ):
            mystr = "\n/*** Begin Subroutine ***/\n"
            fobj.write(mystr);

            mystr = """
void lp2fem1d_dshapefunc_%i
(
    Index elem_n,
    double *u,   
    double *v,   
    double *b    
)
{
    Index i;
    b  = b + %i*elem_n -1;
    u  = u + %i*elem_n -1;
    v  = v + elem_n;
""" % (p, p-1, p+1 )
            fobj.write(mystr);
            s = p
            t = p - 1
            coeff1 = '_C{0:02d}'.format(s)
            s = s - 1
            coeff2 = '_C{0:02d}'.format(s)
            mystr = """
    *b = %s**u; u--; b--;
    *b = %s**u; u--; b--;
""" %( coeff1, coeff2 )
            fobj.write(mystr);

            for k in range(0,p-3):
                s = s - 1
                t = t - 1
                coeff1 = '_C{0:02d}'.format(s)
                coeff2 = '_D{0:02d}'.format(t)
                mystr = "%s*b = %s**u + %s**(b+2); u--; b--;\n" %(ST, coeff1, coeff2)
                fobj.write(mystr);
            t = t - 1
            coeff1 = '_D{0:02d}'.format(t)
            t = t - 1
            coeff2 = '_D{0:02d}'.format(t)
            mystr = "%s*v = %s**(b+1) + %s**(b+2) + *u + *(u-1), v--;\n" % (ST, coeff2, coeff1) 
            fobj.write(mystr);
            mystr = "%s*v = %s**(b+1) - %s**(b+2) - *u + *(u-1); v--;\n" % (ST, coeff2, coeff1)
            fobj.write(mystr);

            mystr = """
    for (i=1; i<elem_n; i++)
    {
"""
            fobj.write(mystr);
            mystr = """
        u -= 2;
"""
            fobj.write(mystr);
            s = p
            t = p - 1
            coeff1 = '_C{0:02d}'.format(s)
            s = s - 1
            coeff2 = '_C{0:02d}'.format(s)
            mystr = """
        *b = %s**u; u--; b--;
        *b = %s**u; u--; b--;
""" %( coeff1, coeff2 )
            fobj.write(mystr);

            for k in range(0,p-3):
                s = s - 1
                t = t - 1
                coeff1 = '_C{0:02d}'.format(s)
                coeff2 = '_D{0:02d}'.format(t)
                mystr = "%s*b = %s**u + %s**(b+2); u--; b--;\n" %(ST2, coeff1, coeff2)
                fobj.write(mystr);
            t = t - 1
            coeff1 = '_D{0:02d}'.format(t)
            mystr = "%s*v = *(v+1) - 2*%s**(b+2) - 2**u, v--;\n" % (ST2, coeff1) 
            fobj.write(mystr);
            mystr = "%s}\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return

def lp2fem1d_zshapefunc_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);

        # Printing the coefficients
        s = N;
        t = N - 2
        for k in range(0, N-1):
            coeff = '_C{0:02d}'.format(s)
            num_str = format(math.sqrt(2*(2*s-1)), num_format )
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
            s = s - 1
        for k in range(0, N-3):
            coeff = '_D{0:02d}'.format(t)
            num_str = format(math.sqrt(2*(2*t-1))/math.sqrt(2*(2*t+3)), num_format )
            mystr   = "#define %s %s\n" %( coeff, num_str )
            fobj.write(mystr);
            t = t - 1
        coeff = '_D{0:02d}'.format(t)
        num_str = format(1/math.sqrt(2*(2*t+3)), num_format )
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        t = t - 1
        coeff = '_D{0:02d}'.format(t)
        num_str = format(1/math.sqrt(2*(2*t+3)), num_format )
        mystr   = "#define %s %s\n" %( coeff, num_str )
        fobj.write(mystr);
        # function enumeration for the header file
        for p in range( 2, N+1 ):
            mystr = \
"void lp2fem1d_zshapefunc_%i( Index elem_n, MATLIB_COMPLEX *u, MATLIB_COMPLEX *v, MATLIB_COMPLEX *b );\n" %p
            fobj.write(mystr);

        for p in range( 2, N+1 ):
            mystr = "\n/*** Begin Subroutine ***/\n"
            fobj.write(mystr);

            mystr = """
void lp2fem1d_zshapefunc_%i
(
    Index           elem_n,
    MATLIB_COMPLEX *u,   
    MATLIB_COMPLEX *v,   
    MATLIB_COMPLEX *b    
)
{
    Index i;
    b  = b + %i*elem_n -1;
    u  = u + %i*elem_n -1;
    v  = v + elem_n;
""" % (p, p-1, p+1 )
            fobj.write(mystr);
            s = p
            t = p - 1
            coeff1 = '_C{0:02d}'.format(s)
            s = s - 1
            coeff2 = '_C{0:02d}'.format(s)
            mystr = """
    *b = %s**u; u--; b--;
    *b = %s**u; u--; b--;
""" %( coeff1, coeff2 )
            fobj.write(mystr);

            for k in range(0,p-3):
                s = s - 1
                t = t - 1
                coeff1 = '_C{0:02d}'.format(s)
                coeff2 = '_D{0:02d}'.format(t)
                mystr = "%s*b = %s**u + %s**(b+2); u--; b--;\n" %(ST, coeff1, coeff2)
                fobj.write(mystr);
            t = t - 1
            coeff1 = '_D{0:02d}'.format(t)
            t = t - 1
            coeff2 = '_D{0:02d}'.format(t)
            mystr = "%s*v = %s**(b+1) + %s**(b+2) + *u + *(u-1), v--;\n" % (ST, coeff2, coeff1) 
            fobj.write(mystr);
            mystr = "%s*v = %s**(b+1) - %s**(b+2) - *u + *(u-1); v--;\n" % (ST, coeff2, coeff1)
            fobj.write(mystr);

            mystr = """
    for (i=1; i<elem_n; i++)
    {
"""
            fobj.write(mystr);
            mystr = """
        u -= 2;
"""
            fobj.write(mystr);
            s = p
            t = p - 1
            coeff1 = '_C{0:02d}'.format(s)
            s = s - 1
            coeff2 = '_C{0:02d}'.format(s)
            mystr = """
        *b = %s**u; u--; b--;
        *b = %s**u; u--; b--;
""" %( coeff1, coeff2 )
            fobj.write(mystr);

            for k in range(0,p-3):
                s = s - 1
                t = t - 1
                coeff1 = '_C{0:02d}'.format(s)
                coeff2 = '_D{0:02d}'.format(t)
                mystr = "%s*b = %s**u + %s**(b+2); u--; b--;\n" %(ST2, coeff1, coeff2)
                fobj.write(mystr);
            t = t - 1
            coeff1 = '_D{0:02d}'.format(t)
            mystr = "%s*v = *(v+1) - 2*%s**(b+2) - 2**u, v--;\n" % (ST2, coeff1) 
            fobj.write(mystr);
            mystr = "%s}\n" % ST
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return
#------------------------------------------------------------------------------#

def prjLP2FEM_ShapeFunc_(N, file_name, PRECISION):
    
    num_format = '0.%sf' % PRECISION
    ST  = "    "
    ST2 = "%s%s" % (ST, ST)
    with open ( file_name, 'w') as fobj:
        mystr = "/* Printed with precision 0.%sf.*/\n" % PRECISION
        fobj.write(mystr);

        for k in range(0, N-1):
            coeff1 = '_E{0:02d}'.format(k)
            coeff2 = '_F{0:02d}'.format(k)
            tmp =  1.0/math.sqrt(4*k+6);
            num_str1 = format(tmp/(k+2.5), num_format )
            num_str2 = format(-tmp/(k+0.5), num_format )
            mystr = "#define %s %s\n" %( coeff1, num_str1 )
            fobj.write(mystr);
            mystr = "#define %s %s\n" %( coeff2, num_str2 )
            fobj.write(mystr);

        for p in range(2, N+1):
            mystr = \
"void fem1d_zprjLP2FEM_ShapeFunc_%i( Index elem_n, MATLIB_COMPLEX *u, MATLIB_COMPLEX *Pv, MATLIB_COMPLEX *Pb);\n" %p
            fobj.write(mystr);

        for p in range(2, N+1):
            mystr = "\n/*** Begin Subroutine ***/\n"
            fobj.write(mystr);
            mystr = """
void fem1d_zprjLP2FEM_ShapeFunc_%i
(
    Index           elem_n,
    MATLIB_COMPLEX* u, 
    MATLIB_COMPLEX* Pv,   
    MATLIB_COMPLEX* Pb    
)
{
    Index i;
    double tmp = 0;
    for (i=0; i<elem_n; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
""" % p
            fobj.write(mystr);
            for k in range(0, p-1):
                coeff1 = '_E{0:02d}'.format(k)
                coeff2 = '_F{0:02d}'.format(k)
                mystr = "%s*Pb = %s**(u+2) + %s**u; Pb++; u++;\n" %(ST2, coeff1, coeff2)
                fobj.write(mystr);
            mystr = "%su += 2;\n" %ST2
            fobj.write(mystr);

            mystr = """
    }
    *Pv = tmp;
"""
            fobj.write(mystr);
            mystr = "}\n"
            fobj.write(mystr);
            mystr = "/*** End of subroutine ***/\n"
            fobj.write(mystr);
    return

            



#------------------------------------------------------------------------------#

if __name__ == '__main__':
    
    file_name = 'test.c'
    PRECISION = '20'
    # fem1d_zlp_snorm2_d_(11, file_name, PRECISION)
    # fem1d_dlp_snorm2_d_(11, file_name, PRECISION)
    # fem_ddot_(17, file_name)
    # fem_zdot_(17, file_name)
    # fem_shapeFunc2lp_(10, file_name, '20')
    # lp2fem1d_dshapefunc_(10, file_name, PRECISION)
    lp2fem1d_zshapefunc_(10, file_name, PRECISION)
    # prjLP2FEM_ShapeFunc_(10, file_name, PRECISION)
    




