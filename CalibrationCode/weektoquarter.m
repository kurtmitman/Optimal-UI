function qt=weektoquarter(x)
    T=floor(length(x)/13);
    qt=zeros(T,1);
    for t=1:T
        qt(t)=mean(x(13*(t-1)+1:13*t));
    end
    
