%% 2019862s

%% This function obtains the numerical approximations
%% to the solutions of the given Poission equation
%% using the Gauss-Jacobi iterative method in polar
%% coordinates.

function [r,t,u,niter]=PolarPoisson(n,m,tol)

%% Parameters and Grid
% Step size in radial direction
deltaR=2/n;
% Step size in azimuthal direction
deltaTheta=(2*pi)/m;
% Radius goes from 1 to 3
r=1:deltaR:3;
% Theta goes from 0 to 2*pi
t=0:deltaTheta:(2*pi);
% Source term
f=-2;
% For loop, which populates the radial partition
ri=zeros(n-1,m-1);
for matDimR = 1:m-1
    ri(:,matDimR)=r(2:n)';
end;
%% Initial guess for Gauss-Jacobi
u=zeros(n+1,m+1);
uold=ones(n+1,m+1);
% Counter for the number of iterations
niter=0;

%% Gauss-Jacobi iteration based on FD equations
i=2:n;
j=2:m;
while norm(u-uold)>tol
    niter=niter+1;
    uold=u;
%%   Boundary conditions
    % Left vertical BC
    u(1,:)=sin(t);
    % Right vertical BC
    u(n+1,:)=cos(3*t);
    % Periodic BCs
    u(i,1)=((-2*f.*ri(:,1).^2.*deltaR.^2.*deltaTheta.^2)...
        +(2.*deltaTheta.^2.*ri(:,1).^2.*(u(i+1,1)+u(i-1,1)))...
        +(deltaR.*deltaTheta.^2.*ri(:,1).*(u(i+1,1)-u(i-1,1)))...
        +(2.*deltaR.^2.*(u(i,2)+u(i,m))))./...
        (4.*ri(:,1).^2.*deltaTheta.^2+4.*deltaR.^2);
    u(i,m+1)=u(i,1);
%%   Internal Nodes
    u(i,j)=((-2*f.*ri.^2*deltaR.^2.*deltaTheta.^2)...
        +(2.*deltaTheta.^2.*ri.^2.*(u(i+1,j)+u(i-1,j)))...
        +(deltaR.*deltaTheta.^2.*ri.*(u(i+1,j)-u(i-1,j)))...
        +(2.*deltaR.^2.*(u(i,j+1)+u(i,j-1))))./...
        (4.*ri.^2.*deltaTheta.^2+4.*deltaR.^2);
end

%% Return solution
r
t
u=u(1:n+1,1:m+1)'
niter

end