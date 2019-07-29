function resid=find_fsx(X,fss,L,A,psi,w,h,b,bet,eps,del,sig,p)

resid=zeros(3,1);
f=X(1);
s=X(2);
x=X(3);
%A=X(4);
%psi=X(5);

sf=s*f;
xf=x*f;

resid(1)=u(w,sig)-u(b+h,sig)+bet*( (1-eps)*c(s,A,psi)+eps*c(x,A,psi)+...
    ((1-eps)*(1-sf)-del)*cp(s,A,psi)/f+eps*(1-xf)*cp(x,A,psi)/f)-cp(s,A,psi)/f;

resid(2)=u(w,sig)-u(h+p,sig)+bet*( c(x,A,psi)+(1-xf)*cp(x,A,psi)/f-del*cp(s,A,psi)/f)-cp(x,A,psi)/f;

resid(3)=fss*(1-L)-s*f*del*L*(1-eps)/(1-(1-s*f)*(1-eps))-x*f*(1-L-del*L*(1-eps)/(1-(1-s*f)*(1-eps)));

% Esb=-0.9;
% Ese=-0.16;
% denom=(bet*((1-eps)*(1-s*f-del)-1)*(1-bet*(1-xf))-bet*bet*del*eps*(1-xf));
% GG=c(x,A,psi)+(1-x*f)*cp(x,A,psi)/f-c(s,A,psi)-(1-s*f)*cp(s,A,psi)/f;
% 
% elas_sb=(f/(s*cpp(s,A,psi)))*(b*up(h+b,sig)*(1-bet*(1-xf)))/denom;
% if(elas_sb>-0.3 || elas_sb<-0.9)
%     resid(4)=1-elas_sb/Esb;
% else
%     resid(4)=0;
% end
% 
% elas_se=(eps*f/(s*cpp(s,A,psi)))*bet*GG*(1-bet*(1-x*f))/denom;
% 
% if(elas_se>-0.1 || elas_se<-0.6)
%     resid(5)=1-elas_se/Ese;
% else
%     resid(5)=0;
% end

%resid(4)=(1-elas_sb/Esb)^20;
%resid(5)=(1-elas_se/Ese)^20;

%[elas_sb elas_se resid(4) resid(5)]