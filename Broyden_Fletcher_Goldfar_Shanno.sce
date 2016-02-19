//Francis Brylle G. Sinco
//MS Applied Mathematics
//Math 288
//University of the Philippines - Diliman
//12 March 2011


//Broyden-Fletcher-Goldfar-Shanno(BFGS) Algorithm with Line Search Algorithm using the Strong Wolfe Conditions applied to minimize the Rosenbrock Function

//INITIALIZATION OF QUANTITIES

k=1;
M=1000;
alpha=10e-4;//alpha
b=0.8;//beta
rho=0.5;
epsilon=10e-6;
x0=[1.2;1.2];//initial guess for the minimizer
B0=eye(2,2);//initial guess for the approximate inverse of the Hessian
I=eye(2,2);


//FUNCTION DEFINITIONS
function z = f(x)//Rosenbrock Function
  z = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
  //z = 2*x(1)^2 + 2*x(1)*x(2) + x(2)^2 + x(1) - x(2);//dummy function
endfunction;

function g = gradient(x)
  g = [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
  //g = [4*x(1)+2*x(2)+1; 2*x(1)+2*x(2)-1];//dummy gradient
endfunction;


//MAIN
xstar = [1;1]; //actual minimizer
x=x0;
B=B0;
while (norm(gradient(x)) > epsilon)//main loop
  d = -B*gradient(x);
  t=1;//initial guess for the steplength
  while ((f(x+t*d) > f(x)+a*t*gradient(x)'*d) | (abs(gradient(x+t*d)'*d) > b*abs(gradient(x)'*d))) //Strong Wolfe conditions loop
    t = rho*t;  
  end;  
  if (abs(t)>epsilon) then
    temp = gradient(x);
    x = x + t*d;
    s = t*d;
    y = gradient(x)-temp;
    P = y'*s;
    Q = s*y';
    R = y*s';
    S = s*s';
    B = [I-Q/P]*B*[I-R/P] + S/P;
    k=k+1;
  else
    break;
  end;
end;

t//prints final step length used before termination
k//prints the latest iteration number before termination
x//prints the final answer (computed minimizer)
err = norm(x-xstar) //prints the absolute difference between the actual and the computed minimizer

