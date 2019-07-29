function [Z Pi] = tauchen(mu,sigma,rho,lambda,znum)

m = zeros(1,znum);
Pi= zeros(znum);
sigma_z = sigma / sqrt(1-rho^2);

% Define the bounds on the grid%
zmin=mu-lambda*sigma_z;
zmax=mu+lambda*sigma_z;

%Assign the quadrature points
Z = linspace(zmin,zmax,znum)';

m(2:znum)=(Z(1:znum-1)+Z(2:end))/2;
m(1)=-inf; %truncate the distribution below
m(znum+1)=inf; %truncate the distribution above 

%size(m)
%size(Z)



upper = repmat(m(2:end),znum,1)-(1-rho)*mu-rho*repmat(Z,1,znum);
lower = repmat(m(1:znum),znum,1)-(1-rho)*mu-rho*repmat(Z,1,znum);
%Pi=normcdf(upper/sigma)-normcdf(lower/sigma);
Pi=normcdf(upper,0,sigma)-normcdf(lower,0,sigma);

%Pi=eye(znum);