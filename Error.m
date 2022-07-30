function [ error ] = Error(U_half,n_N )
error=zeros(1,1);
%%%%Calculate the relative error value
    for i=1:4
        for j=1:n_N-1
            error(i,j)=U_half(i,j+1)-U_half(i,j);
        end
    end
end