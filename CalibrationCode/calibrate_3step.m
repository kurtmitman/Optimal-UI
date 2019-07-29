function resid=calibrate_3step(X,yseq,Pi,znum,z)

resid=zeros(4,1);
%chi=0.5839+0.05*X(1);
%xi=0.1083+0.005*X(2);
%h=0.55+0.01*X(3);
%A=0.03+0.005*X(4);
%psi=X(5);
%psi=-1;
%sig=2+0.05*X(4);

% chi=0.492254+0.1*X(1);
% xi=0.103394+.01*X(2);
% h=0.5736415+0.01*X(3);
% A=0.0006860+0.00005*X(4);
% psi=4.613825+0.01*X(5);

% chi=0.4446+0.1*X(1);
% xi=0.0275+.01*X(2);
% h=0.599044+0.01*X(3);
% A=0.0000003+0.000000005*X(4);
% psi=10.81281+0.1*X(5);

% A=0.001485132203114;
% psi=3.785780171146302;
% h=0.580464837730579;
% chi=0.492046398530831;
% xi=0.114130426526315;

%A=0.001485132203114+0.005*X(3);

%to do 0.9 calib

%A=0.0001*exp(X(3));
%psi=3.785780171146302+X(4);
%h=0.58464837730579+0.01*X(2);
%xi=0.05+0.005*X(1);
%chi=0.492+0.1*X(5);

A=0.002*exp(X(3));
psi=exp(X(4));
h=0.6+0.01*X(2);
xi=0.05+0.005*X(1);
chi=0.414+0.1*X(5);
%b=0.38+0.1*X(6);

b=0.4;


%chi=0.492046398530831+0.1*X(5);

%xi=0.14130426526315+0.005*X(1);


%chi=0.3995;
% A=0.0001485132203114+0.00005*X(4);
% psi=8+0.25*X(5);
% h=0.580464837730579+0.025*X(3);
% chi=0.492046398530831+0.05*X(1);
% xi=0.05+0.01*X(2);

p=0.0;
% psi=real(psi);
% h=real(h);
% xi=real(xi);
% chi=real(chi);
% A=real(A);

[chi xi h psi A]
sig=1;
k=.58;
bet=.99^(1/12);
del=0.0081;

Lss=0.945;
fss=0.139;
ssf=fss;
eps=1/26;
%b=0.4;

thetss=0.634;

q=fss/thetss;
tau=0.023;

w=1-(k/q)*(1-bet*(1-del))-tau;


fap0=[0.3; 0.7; 0.8];
options=optimset('MaxFunEvals',1000000,'MaxIter',10000,'Display','off');

[FAX val]=fsolve(@(x) find_fsx(x,fss,Lss,A,psi,w,h,b,bet,eps,del,sig,p),fap0,options);

f=FAX(1);
s=FAX(2);
x=FAX(3);
%A=FAX(4);
%psi=FAX(5);
%[f s x A psi]
xf=x*f;
sf=s*f;

Esb=-0.9;
%Esb=-0.5;

Ese=-0.1;
denom=(bet*((1-eps)*(1-s*f-del)-1)*(1-bet*(1-xf))-bet*bet*del*eps*(1-xf));
GG=c(x,A,psi)+(1-x*f)*cp(x,A,psi)/f-c(s,A,psi)-(1-s*f)*cp(s,A,psi)/f;

elas_sb=(f/(s*cpp(s,A,psi)))*(b*up(h+b,sig)*(1-bet*(1-xf)))/denom;
resid(4)=(1-elas_sb/Esb);

elas_se=(eps*f/(s*cpp(s,A,psi)))*bet*GG*(1-bet*(1-x*f))/denom;

%resid(5)=(1-elas_se/Ese);
keyboard
eps=1/26*ones(znum,1);
eps(1:8)=1/59; %1/59; %1/52;
eps(9:14)=1/46; %1/46; %1/39; 


Dyn=0;

if(Dyn==1)
   eps(1:34)=1/59; 
   eps(34:42)=1/46; 
    
end
const=0;
if(const==1)
    b=0.417697;
    eps(:)=0.0373989;
    tau=0.0213742;
end
%eps=1./(1.01*1./eps);

%[h xi chi A psi f s x]
X0=[0.634*linspace(0.9,1.1,znum)';max(min(s*linspace(0.9,1.1,znum),0.999),0)'; min(x*linspace(0.98,1.02,znum),0.999)'];
%X0=[0.634*linspace(0.9,1.1,znum)';linspace(0,1,znum)'; 1.1*linspace(0,1,znum)'];
options=optimset('MaxFunEvals',100000,'MaxIter',500,'Display','off');
[tsx val]=fsolve(@(x) eq_resid_with_svx(x,xi,h,chi,A,psi,tau,eps,b,sig,p,Pi,znum,z),X0);

t=tsx(1:znum);
s=tsx(znum+1:2*znum);
x=tsx(2*znum+1:3*znum);
%eq_resid_with_svx(tsx,xi,h,chi,A,psi,tau,eps,b,sig,p,Pi,znum,z)
%eq_resid_with_svx([t;s;s],xi,h,chi,A,psi,tau,eps,b,sig,p,Pi,znum,z)


f=t./((1+t.^chi).^(1/chi));
q=1./((1+t.^chi).^(1/chi));

%f=t.^(1-chi);
%q=f./q;


w=upinv((1-xi)*cp(s,A,psi)./(k*xi*t),sig);

%[s-x, s, x];

L0=.942;
D0=(1-L0)*.8;
Lpr=L0;
Dpr=D0;
ynum=length(yseq)/10;
Lvec=zeros(ynum,1)';
tvec=zeros(ynum,1)';
xvec=zeros(ynum,1)';
zvec=zeros(ynum,1)';
Dvec=zeros(ynum,1)';
fvec=zeros(ynum,1)';
evec=zeros(ynum,1)';
svec=zeros(ynum,1)';
wvec=zeros(ynum,1)';
frawvec=zeros(ynum,1)';
trawvec=zeros(ynum,1)';
welfare=zeros(ynum,1)';
budget=0;
%welfare=0;
for N=1:ynum
    if(N==5200)
        budget=0;
%        welfare=0;
    end
    L0=Lpr;
    D0=Dpr;
    zz=yseq(N);
    zvec(N)=z(zz);
    xx=x(zz);
    tt=t(zz);
    ss=s(zz);
    ff=f(zz);
    qq=q(zz);
    ww=w(zz);
    Lpr=(1-del)*L0+ff*(ss*D0+xx*(1-L0-D0));
    Dpr=(1-eps(zz))*(del*L0+(1-ss*ff)*D0);
    Lvec(N)=Lpr;
    Dvec(N)=Dpr;
    tvec(N)=(ss*D0+xx*(1-L0-D0))*tt/(1-L0);
    frawvec(N)=ff;
    fvec(N)=(ss*D0+xx*(1-L0-D0))*ff/(1-L0);
    evec(N)=eps(zz);
    trawvec(N)=tt;
    svec(N)=ss;
    xvec(N)=xx;
    wvec(N)=ww;
    budget=budget+(tau*Lpr-b*(Dpr/(1-eps(zz))));
%    welfare=welfare+bet^(N-5200)*Lpr*u(ww,sig)+(Dpr/(1-eps(zz)))*u(b+h,sig)+(1-Lpr-(Dpr/(1-eps(zz))))*u(b+p,sig)-D0*c(ss,A,psi)-(1-L0-D0)*c(xx,A,psi);
    welfare(N)=Lpr*u(ww,sig)+(Dpr/(1-eps(zz)))*u(b+h,sig)+(1-Lpr-(Dpr/(1-eps(zz))))*u(h+p,sig)-D0*c(ss,A,psi)-(1-L0-D0)*c(xx,A,psi);
end

vvec=tvec.*(1-Lvec);

tq=weektoquarter12(tvec(5200:end));
uq=weektoquarter12(1-Lvec(5200:end));
fq=weektoquarter12(fvec(5200:end));
wq=weektoquarter12(wvec(5200:end));
sq=weektoquarter12(svec(5200:end));
xq=weektoquarter12(xvec(5200:end));
zq=weektoquarter12(zvec(5200:end));
eq=weektoquarter12(1./evec(5200:end));
vq=weektoquarter12(vvec(5200:end));
rq=weektoquarter12(b./wvec(5200:end));
dq=weektoquarter12(Dvec(5200:end));

%piq=weektoquarter(profit(5200:end));
%exq=weektoquarter(expend(5200:end));
%bq=weektoquarter(b(5200:end));

smoothing=1600;
%smoothing=10000;
hpzq=hpfilter(zq,smoothing);
hptq=hpfilter(tq,smoothing);
hpuq=hpfilter(uq,smoothing);
hpwq=hpfilter(wq,smoothing);
hpsq=hpfilter(sq,smoothing);
hpxq=hpfilter(xq,smoothing);
hpeq=hpfilter(eq,smoothing);
hpfq=hpfilter(fq,smoothing);
hpvq=hpfilter(vq,smoothing);

%hppiq=hpfilter(piq,smoothing);
%hpexq=hpfilter(exq,smoothing);
%hpbq=hpfilter(bq,smoothing);


logdevzq=log(zq)-log(hpzq);
logdevwq=log(wq)-log(hpwq);
logdevtq=log(tq)-log(hptq);
logdevuq=log(uq)-log(hpuq);
logdevsq=log(sq)-log(hpsq);
logdevxq=log(xq)-log(hpxq);
logdeveq=log(eq)-log(hpeq);
logdevfq=log(fq)-log(hpfq);
logdevvq=log(vq)-log(hpvq);

%logdevpiq=log(piq)-log(hppiq);
%logdevexq=log(exq)-log(hpexq);
%logdevbq=log(bq)-log(hpbq);


 statsvec=[logdevzq logdevuq logdevvq logdevtq  logdevfq logdevwq  logdevsq logdevxq logdeveq];% logdevpiq logdevexq logdevbq];
 meanstatsvec=[hpzq uq hpvq hptq hpfq hpwq  hpsq hpxq hpeq dq];% hppiq hpexq hpbq]; 



bhat=regress(logdevwq,[logdevzq ones(length(logdevzq),1)]);
std(logdevzq);

resid(2)=(1-std(logdevtq)/.259);
resid(5)=(1-std(logdevuq)/.125);
resid(3)=(1-mean(fq)/.139);
%resid(6)=(1-mean(rq)/.4);
resid(1)=(1-bhat(1)/.449);
%resid(2)=(1-mean(tq)/.634);
%resid(3)=(1-mean(fq)/ssf);


%resid(5)=(1-bhat(1)/.449);
avef1=mean(fvec(5200:end));
avefraw=mean(frawvec(5200:end));
aves1=mean(svec(5200:end));
avex1=mean(xvec(5200:end));
avee1=mean(evec(5200:end));
meandur1=(1+(1-aves1*avef1)*avee1/(avex1*avef1))/(1-(1-aves1*avef1)*(1-avee1));
meandur2=(1+(1-aves1*avefraw)*avee1/(avex1*avefraw))/(1-(1-aves1*avefraw)*(1-avee1));
meandur3=mean((1+(1-svec(5200:end).*frawvec(5200:end)).*evec(5200:end)./(xvec(5200:end).*frawvec(5200:end)))./(1-(1-svec(5200:end).*frawvec(5200:end)).*(1-evec(5200:end))));
%resid(5)=(1-meandur1/13.2);
mean(welfare(5200:end))
budget


% durtime=zeros(ynum-5304+1,1);
% durweight=zeros(ynum-5304+1,1);
% for utime=5200:ynum-104
%     durtime=sum(
%     
% end

 mean(meanstatsvec)'
 std(statsvec)'
 corr(statsvec)
[chi xi h psi A b]
[std(logdevtq) mean(tq) elas_sb bhat(1) std(logdevuq) mean(1./fvec(5200:end))]


%[mean(uq) std(logdevuq) mean(vq) std(logdevfq) bhat(1) elas_se]
mean(Dvec./(1-Lvec));

%sprintf('%15.10f',mean(uq))
%sprintf('%15.10f',mean(tq))
%sprintf('%15.10f',mean(vq))


% meanL=mean(Lvec);
% meanD=mean(Dvec);
% budget1=budget;
% 
