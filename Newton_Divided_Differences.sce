//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 271.1
//University of the Philippines - Diliman
//15 July 2010



//This program implements the Newton Divided Differences Interpolation to a function, defined and continuous on a closed interval [a,b]

//FUNCTION DEFINITIONS

//function to be interpolated
function y = f(x)
  y = cos(x/2) + sin(x);
endfunction;



//Newton Divided Differences polynomial
function P = newton(x) 
  X=zeros(m,1);  //vector containing the n+1 nodes
  X(1)=a;
  for i=1:n
    X(i+1) = a + i*h;
  end;
  
  Y=zeros(m,1); //vector containing the n+1 function values
  for i=1:m
    Y(i)=f(X(i));
  end;
  
  for j=2:m
    for i=1:m-j+1
      Y(i,j)=(Y(i+1,j-1)-Y(i,j-1))/(X(i+j-1)-X(i));
    end;
  end;
    
  P=Y(1,m);
  for i=2:m
    P=P*(x-X(i))+Y(i,m-i+1); 
  end;
endfunction;



//interpolation error
function e = err(x)
  e = abs(f(x) - newton(x));
endfunction;


//MAIN 
a=-4;          //left endpoint(x_0)
b=4;           //right endpoint(x_n)
n=4;           //number of sub-intervals; degree of Newton polynomial
h=(b-a)/n;     //equal distance between nodes
m=n+1;         //number of nodes(or sample x values)
x=poly(0,'x');

X=zeros(m,1);  //vector containing the n+1 nodes
X(1)=a;
for i=1:n
  X(i+1) = a + i*h;
end;
  
Y=zeros(m,1); //vector containing the n+1 function values
for i=1:m
  Y(i)=f(X(i));
end;

disp([X Y]);   //displays the nodes (first column) and function values (second column)

for j=2:m      //binomial tree for the Newton's Divided Differences
  for i=1:m-j+1
    Y(i,j)=(Y(i+1,j-1)-Y(i,j-1))/(X(i+j-1)-X(i));
  end;
end;

disp([X Y]);  //displays the binomial tree for the Newton's Divided Differences
  
newton(x)     //displays the Newton Divided Differences Interpolating Polynomial

//////////////////////////////END////////////////////////////////////////
