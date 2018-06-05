function [x,y,u,niter]=CartesianPoisson(n,tol)

%% Parameters and Grid
h=1/n;
x=0:h:1;
y=0:h:1;

%% Exact Solution - useful for validation
% [X,Y]=meshgrid(x,y);
% U=-X.^2+X;

%% Initial guess for Gauss-Jacobi
u=zeros(n+1,n+3);
uold=ones(n+1,n+3);
niter=0;

%% Source term
f=-2*ones(n-1,n+1);

%% Gauss-Jacobi iteration based on FD equations
i=2:n;
j=2:n+2;
while norm(u-uold)>tol
    niter=niter+1;
    uold=u;
    
    % Boundary conditions
    u(1,:)=0;
    u(n+1,:)=0;
    u(:,1)=u(:,3);
    u(:,n+3)=u(:,n+1);
    
    % Internal nodes 
    u(i,j)=(f*h^2-u(i-1,j)-u(i+1,j)-u(i,j-1)-u(i,j+1))/(-4);
end

%% Validation if needed
%norm(u(1:n+1,2:n+2)'-U)

%% Return solution
u=u(1:n+1,2:n+2)';
end