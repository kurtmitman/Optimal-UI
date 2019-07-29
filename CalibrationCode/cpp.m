function cpp=cpp(s,A,psi)

cpp=(2+psi)*A*(1./(1-s)).^(3+psi);


%cpp=A*psi*s.^(psi-1);
%cpp=A*psi*(1+psi)*s.^(psi-1);
%cpp=A./((1-s).^2);
%cpp=(A^(1+psi))*psi*(1+psi)*s.^(psi-1);
