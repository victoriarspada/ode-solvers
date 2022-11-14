%% Objective: An ODE system solver using the Heun/Improved Euler Method.
%% Details: Consider the system of 2 ODEs:
% |x1'=f(t,x1,x2), x2'=g(t,x1,x2)|
%
% This m-file is a function which accepts as variables 
% (t0,tN,x0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, h is the stepsize, and x0 is a vector 
% for the initial condition of the system of ODEs |x(t0)=x0|.
%
% The m-file returns a row vector of times and a matrix of 
% approximate solution values (the first row has the approximation for |x1|
% and the second row has the approximation for |x2|). 

function [t,y] = solvesystem(x1, x2, t0,tN,x0,h)
 	j = (tN-t0)/h; % Number of timesteps.
 	t = zeros(1,j); % Initialize the output t matrix.
 	y = zeros(2,j); % Create a square matrix for x1,x2.
 
 	t = linspace(t0, tN, j); % Set values into output t matrix.
 	
 	y(1,1) = x0(1,1); % Insert initial values into output y matrix.
 	y(2,1) = x0(1,2);
 
 	n=2;
 	for (n=1:j-1)
 		k11 = x1(t(1,n),y(1,n); % First slope estimates.
 		k12 = x2(t(1,n),y(2,n));
 
 		un1 = y(1,n) + h*k11; % First y estimates.
 		un2 = y(2,n) + h*k12;
 
 		k21 = x1(t(1,n+1), un1); % Second slope estimates.
 		k22 = x2(t(1,n+1), un2);
 
 		y(1,n+1) = y(1,n) + h*0.5*(k11+k21); % Final y estimates.
 		y(2,n+1) = y(2,n) + h*0.5*(k12+k22);
 	end
 
end