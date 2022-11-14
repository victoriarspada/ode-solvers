%% The following function is an implememntation of Improved Euler's method for a differential equation, with the follwing inputs:
%% f: the input differential equation defined as an inline function. It is defined as f @ (t,y), where t is a parameter and y is the equation we are solving for.
%% t0: the start time of the differential equation.
%% tN: the end time of the differential equation.
%% y0: the value of the function at the start time t0.
%% h: the number of timesteps chosen for the approximation.

%% The solution outputs a solution in two arrays, deltaxi and yi, and graphs them against a solution using the ode45 function.

function [yi, deltaxi] = euler_method_improved(f, t0, tN, y0, h) 
    j = (tN-t0)/h; %% Number of timesteps.
    deltaxi = linspace(t0, tN, j); %% X values in given interval.
    yi = zeros(1,j); %% The template of our final solution to y.
    tol = 1e-8; % Error tolerance [Choice of user].
    
    yi(1,1) = y0 % Initial value of y.
    
    for (i=2:j) % For j points in our solution...
        
        Y = y(1,i-1) + h*f( (i-1)*t, y(1,i-1));
        Z = y(1,i-1) + 0.5*h*f( (i-1)*t, y(1,i-1));
        D = Z-Y;
        
        if abs(D)<tol
            yi(1,i) = Z+D; % If the point has tolerable error, O(h^3)...
        end
        if abs(D)>=tol % If the error is too high...
           while abs(D)>=tol  
               h = 0.9*h*min(max(tol/abs(D),0.3),2);
                Y = y(1,i-1) + h*f( (i-1)*t, y(1,i-1));
                Z = y(1,i-1) + 0.5*h*f( (i-1)*t, y(1,i-1));
                D = Z-Y;
                yi(1,i) = Z+D; % The point now has tolerable error, O(h^3)
            end
        end
        
    end
    
end

% the formula for this stepsize is attempting to achieve O(h^3) error by
% setting the tolerance to e^(-1/8), or (e^(-1/2))^3. This makes the
% improved Euler method have 3 times less error than O(h).