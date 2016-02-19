//Francis Brylle G. Sinco 2005-06538
//MS Applied Mathematics
//Math 271.1
//University of the Philippines - Diliman
//02 October 2010

//This program implements the Lanczos Method for finding the eigenvalues of a given symmetric matrix


n=10; //size of the matrix
b=0; //initial guess for beta
q0=zeros(n,1); //initial guess for the orthogonal vector
q=[1; zeros(n-1,1)]; 
dummy=b*q0;

Z=rand(n,n);
A=Z'*Z; //symmetric positive definite (spd) matrix

Q=zeros(n,n); //matrix whose column vectors are orthogonal each
F=zeros(n,1); // vector of alphas
B=zeros(n,1); //vector of betas

for i=1:n //iteration for the Lanczos Method
  Q(:,i)=q;
  p=A*q;
  alpha=q'*p;
  w=p-(alpha)*(q)-dummy;
  b=sqrt(w'*w);
  dummy=b*q;
  q=w/b;
  F(i)=alpha;
  B(i)=b;
end;

F=diag(F);
B=B(1:n-1);
B=diag(B);
B1=[zeros(1,n); B zeros(n-1,1)]; //lower-diagonal matrix
B2=[zeros(n-1,1) B; zeros(1,n)]; //upper-diagonal matrix
B=B1+B2; 
S=B+F; //tridiagonal matrix
disp(Q);
disp(S);

A0=Q*S*Q'; //approximation of A using the Lanczos method
disp(A);
disp(A0);

//note that since A=QSQ' and Q is orthogonal, then A and S have the same characteristic polynomial, and hence, eigenvalues.

spec(A) //eigenvalues of A
spec(S) //eigenvalues of S
