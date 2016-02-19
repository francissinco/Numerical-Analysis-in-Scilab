//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 288
//University of the Philippines - Diliman
//29 January 2011



// Steepest Descent Algorithm with Backtracking (Armijo) Linesearch Applied to Minimize the Rosenbrock Function

//INITIALIZATION OF QUANTITIES

rho=0.9;
a=0.5;
k=1;
tau=10e-6;//relative error tolerance
x0=[1.2;1.2];//initial guess for the minimizer


//FUNCTION DEFINITIONS
function z = f(x)//Rosenbrock Function
  z = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction;

function g = gradient(x)
  g = [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
endfunction;


//MAIN
xstar = [1;1]; //actual minimizer
x=x0;
while (norm(gradient(x)) > tau)//relative error criterion loop
  //disp(k)//prints the iteration number
  //disp(x)//prints the iterate
  d = -gradient(x);
  t=1;//initial guess for the steplength
  while (f(x+t*d) > f(x)+a*t*gradient(x)'*d) //Armijo condition loop
    t = rho*t;  
  end;  
  //disp(t)//prints the steplength used for each iteration
  x = x + t*d;
  k=k+1;
end;

x//final answer (computed minimizer)
k//the iteration number when the Relative Error Criterion is satisfied
err = norm(x-xstar) //prints the absolute difference between the actual and the computed minimizer


///////////////////////////////END////////////////////////////////////
