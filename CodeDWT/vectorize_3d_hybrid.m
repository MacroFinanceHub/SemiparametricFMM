% This function takes the matrix D0, and reshapes it into a vector
% The vector contains blocks of coefficients according to 
% wavespecs and wavespecs2

function [D0_vectorized_comb,wavespecsD0]=vectorize_3d_hybrid(M,Kjs_row, Kjs_col)
wavespecsD0.Kjs2D=Kjs_col;
wavespecsD0.Kjs1D=Kjs_row;
temprow=[0,cumsum(Kjs_row)];
tempcol=[0,cumsum(Kjs_col)];

block=1;
D0_vectorized=cell(length(Kjs_col)*length(Kjs_row),1);
for i=1:length(Kjs_col)
    for j=1:length(Kjs_row)
     
        D0_vectorized{block}=reshape(M(temprow(j)+1:temprow(j+1),tempcol(i)+1:tempcol(i+1)),Kjs_row(j)*Kjs_col(i),1);
        block=block+1;
    end
end

wavespecsD0.Kj=kron(Kjs_col,Kjs_row);
aux=zeros(size(D0_vectorized,1),1);
D0_vectorized_comb=[];
for block=1:size(D0_vectorized,1)
    aux(block,1)=size(D0_vectorized{block,1},1);
    D0_vectorized_comb=[D0_vectorized_comb;D0_vectorized{block}];
end
