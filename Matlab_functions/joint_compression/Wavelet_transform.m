function [D]=Wavelet_transform(X_raw,wavespecs,Image_matrix)

tic
    
if nargin<3
    Image_matrix=1;
end;

if nargin<2
    wavespecs.wavelet='db4';
    wavespecs.nlevels=floor(log(min(size(X_raw{1})))/log(2));
    wavespecs.ndim=2;
end;

if wavespecs.ndim==2
    if Image_matrix==1
        X=reshape(X_raw(1,:),wavespecs.t1,wavespecs.t2);
        n=size(X_raw,1);
    else
        n=length(X_raw);
        X=X_raw{1}; 
    end;
    if wavespecs.rectangular==1  %%% if rectangular, then wavespecs.nlevels
        %%%% and wavespecs.boundary refer to rows, and then .nlevels2 and
        %%%% .boundary2 refer to columns;
        dwtmode(wavespecs.boundary);
        [C11,S1]=wavedec(X(1,:),wavespecs.nlevels,wavespecs.wavelet);
        C1=zeros(size(X,1),size(C11,2));
        C1(1,:)=C11;
        for j=2:size(X,1);
            [C1(j,:)]=wavedec(X(j,:),wavespecs.nlevels,wavespecs.wavelet);
        end;
        dwtmode(wavespecs.boundary2);
        [C21,S2]=wavedec(C1(:,1),wavespecs.nlevels2,wavespecs.wavelet2);
        C2=zeros(size(C21,1),size(C1,2));
        C2(:,1)=C21;
        for j=2:size(C1,2);
            [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels2,wavespecs.wavelet2);
        end;
        %%%% Should be reordered in blocks of K1 x K2 
        [reorder,wavespecs.Kj]=reorder_rectangular(S1(1:(end-1)),S2(1:(end-1)));
        wavespecs.reorder=reorder;
        wavespecs.Kj1=S1(1:(end-1));
        wavespecs.Kj2=S2(1:(end-1));
        wavespecs.J1=length(wavespecs.Kj1);
        wavespecs.J2=length(wavespecs.Kj2);
        wavespecs.T1=S1(end);
        wavespecs.T2=S2(end);
        wavespecs.T=S1(end)*S2(end);
        C=reshape(C2,1,size(C2,1)*size(C2,2));
        D=zeros(n,length(C));  
    else
        [C]=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
        D=zeros(n,length(C));    
    end
    for i=1:n
        if Image_matrix==1
            X=reshape(X_raw(i,:),wavespecs.t1,wavespecs.t2);
        else
            X=X_raw{i};
        end;
        if wavespecs.rectangular==1  
            dwtmode(wavespecs.boundary);
            C1=zeros(size(X,1),size(C2,2));
            for j=1:size(X,1);
                [C1(j,:)]=wavedec(X(j,:),wavespecs.nlevels,wavespecs.wavelet);
            end;
            C2=zeros(size(C2,1),size(C2,2));
            dwtmode(wavespecs.boundary2);
            for j=1:size(C1,2);
                [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels2,wavespecs.wavelet2);
            end;
            %%%% Should be reordered in blocks of K1 x K2 but leave alone for
            %%%% now
            C=reshape(C2,1,size(C2,1)*size(C2,2));
            D(i,:)=C(reorder);
        else
            D(i,:)=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
        end;
    end;
elseif wavespecs.dim==1
        n=size(X_raw,1);
        X=X_raw(1,:);
        [C]=wavedec(X,wavespecs.nlevels,wavespecs.wavelet);
        D=zeros(n,length(C));    
    for i=1:n
        D(i,:)=wavedec(X_raw(i,:),wavespecs.nlevels,wavespecs.wavelet);
    end;
end;
'Done with Wavelet Transforms ',toc


