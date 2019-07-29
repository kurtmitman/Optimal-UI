function qt=weektoquarter12(x)
    T=floor(length(x)/12);
    qt=zeros(T,1);
    for t=1:T
        qt(t)=mean(x(12*(t-1)+1:12*t));
    end
    
