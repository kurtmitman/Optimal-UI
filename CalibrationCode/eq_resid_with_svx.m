function resid=eq_resid_with_svx(tsx,xi,h,chi,A,psi,tau,eps,b,sig,p,Pi,znum,z)

del = 0.0081;         % Exogenous separation rate
bet = .99^(1/12);     % discount factor
k=.58;

%del=0.0082*exp(-2.4772*log(z));
%Error variance: 0.0026


t=tsx(1:znum);
%s=exp(tsx(znum+1:2*znum))./(1+exp(tsx(znum+1:2*znum)));
%x=exp(tsx(2*znum+1:3*znum))./(1+exp(tsx(2*znum+1:3*znum)));
%s=min(0.99,tsx(znum+1:2*znum));
%x=min(0.999999,tsx(2*znum+1:3*znum));
s=tsx(znum+1:2*znum);
x=tsx(2*znum+1:3*znum);

f=t./((1+t.^chi).^(1/chi));
q=1./((1+t.^chi).^(1/chi));
%f=t.^(1-chi);
%q=f./t;

w=upinv((1-xi)*cp(s,A,psi)./(k*xi*t),sig);


JJ=k./q;
DD=cp(s,A,psi)./f;
XX=cp(x,A,psi)./f;


resid(1:znum)=1-(w+tau+JJ-bet*Pi*((1-del).*JJ))./z;

resid(znum+1:2*znum)=DD-(u(w,sig)-u(h+b,sig)+bet*(1-eps).*(Pi*(c(s,A,psi)+(1-del-s.*f).*DD))+...
    bet*eps.*(Pi*(c(x,A,psi)+(1-x.*f).*XX-del.*DD)));

resid(2*znum+1:3*znum)=XX-(u(w,sig)-u(h+p,sig)+bet*Pi*(c(x,A,psi)+(1-x.*f).*XX-del.*DD));



%resid(znum+1:2*znum)=1-(exp(DD-bet*Pi*(A*log(1./(1-s))+(1-del-alpha).*DD)).*(BB+ch)./x);    
% resid(znum+1:2*znum)=1-(exp(DD-bet*Pi*(log(1./(1-A*s))+(1-del-alpha).*DD)).*(BB+ch)./x);    
 



%resid(2*znum+1:3*znum)=1-(1/A-(1-xi)*x./(xi*k*q))./s;