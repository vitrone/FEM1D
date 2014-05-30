function mxLGLdataIM
% mxLGLdataIM Inverse Transform Matrix for LGL points
%     IM = LGLdataIM( p, xi) finds the (P+1)-by-(p+1) matrix for computing the 
%     inverse Legendre transformation of a vector containing the sampled values of 
%     a function over P+1 LGL points (xi) in the interval [-1, 1]. The highest 
%     degree of the Legendre polynomials inviolved is p. For arbitrary intervals, 
%     create a linear map from [-1,1]-->[a,b]: 
%                   x = J*xi+xm, where J = (b-a)/2, xm = (a+b)/2.
%     
%     Example: 
%     Let v be the Legendre transform over p+1 LGL points.
%     xi = (-1:0.001:1);
%     P = length(xi)-1;
%     IM = mxLGLdataIM( p, xi);
%     u = IM*v; 
%     Here, 'u' is the vector containing values of the function at points xi.

%------------x------------------------------------------------------------------%
