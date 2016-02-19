//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 271.1
//University of the Philippines - Diliman
//16 September 2010



//This program implements the Jacobi, Gauss-Seidel, and Conjugate Gradient methods for finding the solution of a large system Ax=b, where A is a symmetric positive-definite matrix

n=30; //size of the matrix
epsilon=.00001;
p=4;
m=10^p; //number of iterations

x0=ones(n,1) //initial iterate
b=ones(n,1) 

Z=rand(n,n); //random n by n matrix
A=Z'*Z //symmetric positive definite (spd) matrix
spec(A) //eigenvalues of A; all eigenvalues should be positive so that A is spd

D=diag(diag(A)) //diagonal part of A
L=tril(A,-1) //lower-triangular part of A
R=triu(A,1) //upper-triangular part of A

//JACOBI//
function J=jacobi(x)
for k=1:m
  x1=-inv(D)*(L+R)*x + inv(D)*b;
  if(norm(x1-x)<epsilon | norm(A*x1-b)<epsilon) then
    break;
  end;
  x=x1;  
end;
disp(k);
J=x;
endfunction;

//GAUSS-SEIDEL
function GS=gaussseidel(x)
for k=1:m
  x1=-inv(D+L)*R*x + inv(D+L)*b;
  if(norm(x1-x)<epsilon | norm(A*x1-b)<epsilon) then
    break;
  end;
  x=x1;  
end;
disp(k);
GS=x;
endfunction;

//CONJUGATE GRADIENT
function CG=conjgrad(x)
r=b-A*x0
rho=(r')*r   
  for k=1:m
    if (k==1) then 
       p=r;
       w=A*p;
       alpha=rho/(p'*w);
       x=x+alpha*p;
       r=r-alpha*w;
       Beta=(r'*r)/rho;
       rho=r'*r;
    else
       p=r+Beta*p;
       w=A*p; 
       alpha=rho/(p'*w);
       x=x+alpha*p;
       r=r-alpha*w;
       Beta=(r'*r)/rho;
       rho=r'*r;
     end;  
     
     if(sqrt(rho) < epsilon*sqrt(b'*b)) then
       break;
     end;
  end;
disp(k);
disp(sqrt(rho));  
CG=x;
endfunction;      
    `
    
