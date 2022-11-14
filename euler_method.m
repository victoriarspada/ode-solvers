%% The following function is an implememntation of Euler's method for a differential equation, with the follwing inputs:
%% f: the input differential equation defined as an inline function. It is defined as f @ (t,y), where t is a parameter and y is the equation we are solving for.
%% t0: the start time of the differential equation.
%% tN: the end time of the differential equation.
%% y0: the value of the function at the start time t0.
%% h: the number of timesteps chosen for the approximation.

%% The solution outputs a solution in two arrays, deltax and y, and graphs them against a solution using the ode45 function.

function euler_method(f, t0, tN, y0, h)
    deltax = linspace(t0, tN, h) %% X values in the given interval.
    y = zeros(t0, tN, h) %% The template of our final solution to y.
    j = (tN-t0)/h %% The number of timesteps.
    
    y(1,1) = y0 %Input initial value of y.
    
    for (i=2:j)
        y(1,i) = y(1,i-1) + h*f( (i-1)*t, y(1,i-1)); % Solve for every point over j points.
    end
    disp(deltax)
    disp(y)
   %% Graph the Euler method with ode45 soln
soln0 = ode45(f, [t0, tN], y0)
figure
title('Euler approximation of function y')
ylabel('Approximation y')
xlabel('t')
plot(deltax, y, soln0.x, soln0.y)
legend('euler method', 'ode45')
end


