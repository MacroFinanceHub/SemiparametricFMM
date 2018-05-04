% This function takes the matrix D0, and reshapes it into a vector
% The vector contains blocks of coefficients according to 
% wavespecs and wavespecs2

function [D0_reshaped]=reshape_3d_hybrid(D0_vectorized_comb,wavespecsD0)
Kjs_col=wavespecsD0.Kjs2D;
Kjs_row=wavespecsD0.Kjs1D;
temprow=[0,cumsum(Kjs_row)];
tempcol=[0,cumsum(Kjs_col)];

D0_reshaped=zeros(sum(Kjs_row),sum(Kjs_col));
blocksum=[0,cumsum(wavespecsD0.Kj)];
c=1;
for i=1:length(Kjs_col)
    for j=1:length(Kjs_row)
      D0_reshaped(temprow(j)+1:temprow(j+1),tempcol(i)+1:tempcol(i+1))=reshape(D0_vectorized_comb(blocksum(c)+1:blocksum(c+1),1),Kjs_row(j),Kjs_col(i));
      c=c+1;
    end
end

