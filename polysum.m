function [res] = polysum(a,b)
    %вычисляет сумму двух полиномов(по убывающим степеням)
    la=length(a);
    lb=length(b);
    if(la>lb)
        c=a;
        for i_=0:lb-1 
           c(la-i_)=c(la-i_)+b(lb-i_);
        end 
    else
        c=b;
        for i_= 0:la-1 
            c(lb-i_)=c(lb-i_)+a(la-i_);
        end
    end
    res=c;
    end