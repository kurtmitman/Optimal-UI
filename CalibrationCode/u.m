function u=u(x,s)

if(s==1)
    u=log(x);
else
    u=(x.^(1-s)-1)/(1-s);
end