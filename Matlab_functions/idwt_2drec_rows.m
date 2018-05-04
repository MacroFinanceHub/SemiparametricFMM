function Y=idwt_2drec_rows(infile,wavespecs)

% help file needs to be updated
% dwt_rows(infile,nlevels,wavelet): Compute IDWT on rows of infile, using wavelet and 
%                                   number of levels specified.
%
% Input: infile : matrix, with each row containing wavelet coefficients for function on which to perform IDWT (n x k)
%        nlevels: Number of levels of wavelet transform
%       wavespecs: structure with information on wavelet settings,
%                   including:
%           wavelet: Name of wavelet basis to use
%
% Output: D : matrix with each row containing a function (n x T)
%         K : vector of number of coefficients per level (1 x nlevels)
%         T : number of sampled points per function 
%
% Once can change the extension mode  to periodic by using dwtmode('per').
%tic
if wavespecs.compress==1
    temp=zeros(size(infile,1),length(wavespecs.keep2));
    temp(:,wavespecs.keep2==1)=infile;
    infile=temp;
end;

Y=NaN(size(infile,1),wavespecs.t1*wavespecs.t2);
K1=[wavespecs.Kj1,wavespecs.t1];
K2=[wavespecs.Kj2;wavespecs.t2];

infile(:,wavespecs.reorder)=infile;
for i=1:size(infile,1)
    d=reshape(infile(i,:),sum(wavespecs.Kj2),sum(wavespecs.Kj1));
    y2=NaN(wavespecs.t2,size(d,2));
    dwtmode(wavespecs.boundary2,'nodisp');
    for j=1:size(d,2);
        [y2(:,j)]=waverec(d(:,j),K2,wavespecs.wavelet2);
    end;
    y1=NaN(wavespecs.t2,wavespecs.t1);
    dwtmode(wavespecs.boundary,'nodisp');
    for j=1:size(y2,1);
        y1(j,:)=waverec(y2(j,:),K1,wavespecs.wavelet);
    end;
    Y(i,:)=reshape(y1,1,wavespecs.t1*wavespecs.t2);
    if (i/100==floor(i/100))
        i,toc
    end;
end



