function [ R_e_delta_y ] = Error_Slope( Er,H ,n_N)
R_e=zeros(1,1);
R_delta_e=zeros(1,1);
R_e_delta_y=zeros(1,1);

for jj=1:n_N-1
     R_delta_e(jj)=(H(jj+1)/H(jj));
end
for i=1:4
    for j=1:n_N-2
        R_e(i,j)=(Er(i,j+1)/Er(i,j));

    end

    
    for j=1:n_N-2
          R_e_delta_y(i,j)=abs(log10(R_e(i,j)/R_delta_e(j)));
    end

end

end