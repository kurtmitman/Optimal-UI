function cp=cp(s,A,psi)

cp=A*(1./(1-s)).^(2+psi)-A;

%cp=A*(1./(1-s)).^(2+psi);


%cp=A*s.^(psi);
%cp=A*(1+psi)*s.^(psi);
%cp=A./(1-s);
%cp=(A^(1+psi))*(1+psi)*s.^(psi);
