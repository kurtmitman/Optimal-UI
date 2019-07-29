function c=c(s,A,psi)

if(psi==-1)
   c=A*log(1./(1-s))-A*s;
else
   c=(A/(1+psi))*((1./(1-s)).^(1+psi)-1)-A*s;
end

% if(psi==-1)
%    c=A*log(1./(1-s));
% else
%    c=(A/(1+psi))*((1./(1-s)).^(1+psi)-1);
% end





%c=(A/(1+psi))*s.^(1+psi);
%c=A*s.^(1+psi);
%c=A*log(1./(1-s));
%c=(A*s).^(1+psi);
