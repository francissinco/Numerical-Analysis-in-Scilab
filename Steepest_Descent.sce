//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 288
//University of the Philippines - Diliman
//15 January 2011


// Steepest Descent Algorithm with Backtracking (Armijo) Linesearch

//INITIALIZATION OF QUANTITIES
t=1;//initial guess for the steplength
rho=0.5;
alpha=0.5;
k=1;
tau=10e-8;//relative error tolerance
x0=[0;0];//initial guess for the minimizer


//FUNCTION DEFINITIONS
function z = f(x)
  z = 2*x(1)^2 + 2*x(1)*x(2) + x(2)^2 + x(1) - x(2);
endfunction;

function g = gradient(x)
  g = [4*x(1)+2*x(2)+1; 2*x(1)+2*x(2)-1];
endfunction;


//MAIN
x=x0;
d=-gradient(x);

while norm(gradient(x)) > tau*norm(gradient(x0))//relative error criterion loop
  d = -gradient(x);
  while (f(x + t*d) > f(x)+ alpha*t*gradient(x)'*d) //Armijo condition loop
    t = rho*t;
  end;  
  x = x + t*d;
  k=k+1;
end;

x//final answer (computed minimizer)
k//number of iterations performed


///////////////////////////////END////////////////////////////////////
