% Optimal UI Over the Cycle
% Kurt and Stan
% Matlab Dynare Code

close all

% Initialize endogenous and exogenous variables, and params

var L z b epsilon dur mu nu gam lam thet s phi w f q fp qp c cp cpp uu uup ue uep uepp D alpha x cx cxp cxpp Pi ux profit fhat thethat Dfrac expend welfare unemp reprate eb vac output;

varexo eps;

parameters h chi xi del k bet A rho sigma_z psi_c sigma_c eta tau;


% Assign calibrated values to the parameters


%%%%%%%%% Baseline %%%%%%%%% 
rho=0.9895;
sigma_z=0.0034;
bet=0.999162822640879;
k=0.58;
del=0.0081;
sigma_c=1;
tau_adj=0;
chi=0.492;

chi=0.441994028491430;
xi=0.077765818398739;
h=0.583738591997052;
psi_c=1.911283133102631;
A=0.005265199133822;

tau=     		 -3.8457;
eta=     		 1.04704;


model; 


c = A*((1/(1-s))^(1+psi_c)-1)/(1+psi_c)-A*s;
cp = A*(1/(1-s))^(2+psi_c)-A;
cpp = (2+psi_c)*A/((1-s)^(3+psi_c));
cx = A*((1/(1-x))^(1+psi_c)-1)/(1+psi_c)-A*x;
cxp = A*(1/(1-x))^(2+psi_c)-A;
cxpp = (2+psi_c)*A/((1-x)^(3+psi_c));
ux=log(h);
uu=log(h+b);
uup=1/(h+b);
ue=log(w);
uep=1/w;
uepp=-1/w^2;
f = thet/(1.0+thet^chi)^(1/chi);
q = 1.0/((1.0+thet^chi)^(1/chi));
qp = -thet^(chi-1.0)/((1+thet^chi)^(1/chi+1));
fp = q+thet*qp;
z = rho*z(-1)+sigma_z*eps;
fhat=(s*D(-1)+x*(1-L(-1)-D(-1)))*f/(1-L(-1));
thethat=(s*D(-1)+x*(1-L(-1)-D(-1)))*thet/(1-L(-1));
profit=k/q;
Dfrac=D/(1-L);
expend=b*D/(1-epsilon);
welfare=L*ue+(D/(1-epsilon))*uu+(1-L-(D/(1-epsilon)))*ux-D(-1)*c-(1-L(-1)-D(-1))*cx;
unemp=1-L;
vac=thethat*unemp;

%%%%%
%L LOM
%%%%%
L = (1-del)*L(-1)+f*(s*D(-1)+x*((1-L(-1)-D(-1))));

%%%%%
%D LOM
%%%%%
D = (1-epsilon)*(del*L(-1)+((1-s*f)*D(-1)));

%%%%%
%FE
%%%%%

k/q=exp(z)-w-exp(tau)+bet*(1-del)*k/q(+1);

%%%%%
%NB
%%%%%

xi*uep*k*thet=(1-xi)*cp;


% FOC b
(D-mu*(1-epsilon))*uup=eta*(D);

%%%%%
% FOC w 
%%%%%
gam=(L+mu+nu)*uep-phi*xi*k*thet*uepp;


%%%%%
% FOC L
%%%%%
lam=ue-ux+eta*exp(tau)+bet*(cx(+1)+lam(+1)*(1-del-x(+1)*f(+1))+alpha(+1)*del);


%%%%%
% FOC D
%%%%%
eta*b=uu-ux-alpha+(1-epsilon)*bet*(cx(+1)-c(+1)+lam(+1)*f(+1)*(s(+1)-x(+1))+alpha(+1)*(1-s(+1)*f(+1)));

%%%%%
% FOC s
%%%%%
phi*(xi-1)*cpp=D(-1)*(lam*f-alpha*f-cp)+(cpp/f)*(mu(-1)*((1-epsilon(-1))*(1-s*f)-del)-mu-del*nu(-1));

%%%%%
% FOC x
%%%%%
(1-L(-1)-D(-1))*(cxp-lam*f)=(cxpp/f)*(mu(-1)*epsilon(-1)*(1-x*f)+nu(-1)*(1-x*f)-nu);

%%%%%
% PK E
%%%%%
cp=f*(ue-uu+bet*((1-epsilon)*c(+1)+epsilon*cx(+1)+((1-epsilon)*(1-s(+1)*f(+1))-del)*cp(+1)/f(+1)+epsilon*(1-x(+1)*f(+1))*cxp(+1)/f(+1)));

%%%%%
% PK X
%%%%%
cxp=f*(ue-ux+bet*(cx(+1)+(1-x(+1)*f(+1))*cxp(+1)/f(+1)-del*cp(+1)/f(+1)));


%%%%%
% FOC epsilon
%%%%%

D*(uu-ux-eta*(b)-alpha)+(mu*(1-epsilon))*((cxp-cp)/f-uu+ux)+Pi=0;

%%%%%
% FOC thet
%%%%%
phi*xi*k*uep=fp*(lam*(s*D(-1)+x*(1-L(-1)-D(-1)))-alpha*s*D(-1))+(k*qp/q^2)*(gam-(1-del)*gam(-1))+(cp*fp/f^2)*(mu-mu(-1)*(1-epsilon(-1)-del)+nu(-1)*del)+(cxp*fp/f^2)*(nu-nu(-1)-mu(-1)*epsilon(-1));

output=exp(z)*L;

%Budget
%(1-epsilon)*exp(tau)*L=b*D;

%FOC tau
%L*eta=gam;
Pi*epsilon=0;

eb=0; %b/(1e-6+epsilon);
reprate=(b+h)/w;
dur=1/epsilon;
%FOC H

end;

% Calculate the Analytic Steady State Values
% Use them to initialize dynare
  Lss=0.942412; %0.95349; 
  bss=0.382; %0.25; %0.47;
  epsilonss=0.022;

%  Lss=0.9;

  Gss=0.9855; %0.98907;
  Mss=-0.0032;
  lamss=2.84;
  thetss=1.425;
  phiss=-0.02121303;
  wss=0.957;
  Delta=0.921477985;
  
  fss = thetss/(1+thetss^chi)^(1/chi);
  qss = 1/((1+thetss^chi)^(1/chi));

  qpss = -thetss^(chi-1.0)/((1+thetss^chi)^(1/chi+1));
  fpss = qss+thetss*qpss;
    pss=0.1;

  sss=0.50;
  xss=0.72;

  uxss=log(h);
  uxpss=1/(h);
  uuss=log(h+bss);
  uupss=1/(h+bss);
  uess=log(wss);
  uepss=1/wss;
  ueppss=-1/wss^2;
  Dss=.0321916;
  nuss=-0.00119666;
  alphass=lamss;
  

css = A*((1/(1-sss))^(1+psi_c)-1)/(1+psi_c)-A*sss;
cpss = A*(1/(1-sss))^(2+psi_c)-A;
cppss = (2+psi_c)*A/((1-sss)^(3+psi_c));
cxss = A*((1/(1-xss))^(1+psi_c)-1)/(1+psi_c)-A*xss;
cxpss = A*(1/(1-xss))^(2+psi_c)-A;
cxppss = (2+psi_c)*A/((1-xss)^(3+psi_c));



  etass=0.581;

initval;
  
 % eta=1.04;
 % tau=-4;
  gam = Gss;
  L = Lss;
  mu = Mss;
  b = bss;
  epsilon=epsilonss;
  lam = lamss;
  thet = thetss;
  s=sss;
  phi=phiss;
  w=wss;
  reprate=bss/wss; 
  f=fss;
  q=qss;
  fp=fpss;
  qp=qpss;
  c=css;
  cp=cpss;
  cpp=cppss;
  ux=uxss;
  uu=uuss;
  uup=uupss;
  ue=uess;
  uep=uepss;
  uepp=ueppss;
  D=Dss;
  nu=nuss;
  alpha=alphass;
  x=xss;
  cx=cxss;
  cxp=cxpss;
  cxpp=cxppss;
  Pi=0.0;
  z = 0.1; 
  eps = 0;
 %   zeta=0;
end;
    
shocks;
  var eps = 1;
end;


steady(solve_algo = 2);

stoch_simul(irf=240, order=3) z b epsilon s x L unemp vac Dfrac w fhat thethat profit expend welfare mu nu eb reprate gam dur;

%stoch_simul(periods=52000, order=3,nofunctions) z unemp vac thethat b epsilon s x L Dfrac w mu nu eb  fhat profit expend welfare;