
%Compute Phi_2^C, the wavelet transform



clearvars -except wavespecs wavespecs3 wavespecs3_2

X=zeros(1,91*109*91);
X(1,1)=1;
index1=find(wavespecs3_2.keep==1);
[C]=wavedec3(reshape(X(1,:),wavespecs.t(1),wavespecs.t(2),wavespecs.t(3)),wavespecs3.nlevels,wavespecs.wavelet,wavespecs3.boundary);
temp1=vectorize_3d_wavelet(C);
D(1,:)=temp1(index1);
    for i=4:size(X,2)
        tic
        X=zeros(1,91*109*91);
        X(1,i)=1;
        [C]=wavedec3(reshape(X,wavespecs.t(1),wavespecs.t(2),wavespecs.t(3)),wavespecs3.nlevels,wavespecs.wavelet,wavespecs3.boundary);
        temp1=vectorize_3d_wavelet(C);
        clear C
        [D(i,:)]=temp1(index1);
        toc
    end
    
    save('/Users/miranda/Documents/MATLAB/HCP_Analysis/3D_aftercompression/WavProjection.mat','D','-v7.3')
    clear D