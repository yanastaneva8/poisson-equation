%% 2019862s

%% Main file to run the PolarPoisson function,
%% produce a plot of the rectangular grid, convert to
%% polar coordinates, and plot the solution of the PDE.

[r,theta,u,niter]=PolarPoisson(30,33,10^(-5));

%% This figure plots the rectangular grid from r=1 to r=3, and
%% theta=0 to theta=2*pi in polar coordinates. 
figure
surf(r,theta,u);
    title('Rectangular grid in polar coordinates')
    xlabel('Radius r')
ylabel('Angle theta')

%% This figure converts the rectangular grid in polar coordinates 
%% back to Cartesian and plots the solutions, obtained by
%% the Gauss-Jacobi iteration for solving the PDE.
[r,theta]=meshgrid(r,theta);
[x,y]=pol2cart(theta,r);
figure
surf(x,y,u);
title('Numerical approximations of the solutions to the PDE')
xlabel('x')
ylabel('y')