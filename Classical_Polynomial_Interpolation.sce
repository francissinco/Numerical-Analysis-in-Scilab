//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 271.1
//University of the Philippines - Diliman
//29 July 2010


//This program implements the Classical Polynomial Interpolation to the function f(x)=1/(1+25x^2) on [-1,1] using equidistant points and using Chebyshev points


//THE QUANTITIES USED
//n = degree of the interpolating polynomials
//n+1 = number of data points; degree of the Chebyshev polynomial
//a = left endpoint of the interval
//b = right endpoint of the interval
//B = vector containing the 0th to the nth powers of the variable x
//X = vector containing the equidistant points
//Y = vector containing the function values at the equidistant points
//V = Vandermonde matrix of the equidistant points
//A = vector containing the coefficients of the interpolating polynomial obtained using equidistant points
//U = vector containing the Chebyshev points
//W = vector containing the function values at the Chebyshev points
//Z = Vandermonde matrix of the Chebyshev points
//D = vector containing the coefficients of the interpolating polynomial obtained using Chebyshev points


//FUNCTION DEFINITIONS
//function to be interpolated
function y = f(x)
  y = 1 / (1 + 25*x^2);
endfunction;

//////////////////////for equidistant points//////////////////////////
//interpolating polynomial using equidistant points
function p = classic1(x)
  B=zeros(n+1,1);
  for i=1:n+1
    B(i)=x^(i-1);
  end;
  p=(A')*B;
endfunction;

//interpolation error using equidistant points
function q = err1(x)
  q = abs(f(x) - classic1(x));
endfunction;

////////////////////for Chebyshev data points/////////////////////////
//Chebyshev polynomial
function T = chebyshev(x)
  C=zeros(n+2,1);
  C(1)=1;
  C(2)=x;
  for i=2:n+1
    C(i+1) = 2*x*C(i) - C(i-1);
  end;
  T=C(n+2);
endfunction;

//interpolating polynomial using Chebyshev points
function r = classic2(x)
  B=zeros(n+1,1);
  for i=1:n+1
    B(i)=x^(i-1);
  end;
  r=(D')*B;
endfunction;

//interpolating polynomial using Chebyshev points
function s = err2(x);
  s = abs(f(x) - classic2(x));
endfunction;

///////////////////////////MAIN////////////////////////////////////
a=-1;
b=1;
n=10;
x=poly(0,'x');


//for equidistant points
X=linspace(a,b,n+1)';
Y=zeros(n+1,1);
for i=1:n+1
  Y(i)=f(X(i));
end;

V=zeros(n+1,n+1);
for j=1:n+1
  for i=1:n+1
    V(i,j)=X(i)^(j-1);
  end;
end;

A=zeros(n+1,1);
A=V\Y;
////////////////////////////////


//for Chebyshev data points

U=zeros(n+1,1);
W=zeros(n+1,1);
for i=1:n+1
  U(i)=cos((2*(n+1-i+1)-1)*%pi / (2*(n+1)));
  W(i)=f(U(i)); 
end;

Z=zeros(n+1,n+1);
for j=1:n+1
  for i=1:n+1
    Z(i,j)=U(i)^(j-1);
  end;
end;

D=zeros(n+1,1);
D=Z\Y;

[X Y]
V
A
classic1(x)

[U W]
chebyshev(x)
Z
D
classic2(x)


