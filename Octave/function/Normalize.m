function A_norm = Normalize(A)
%NORMALIZE input the three-dimensional tensor A, normalize it, and output the normalized A_norm
if numel(size(A))==2
    istensor=0;
    ss=A;
elseif  numel(size(A))==3
    istensor=1;
    [N,M,L]=size(A);
    ss=reshape(A,[],L);
else
    error('It must be matrix or tensor!\n');
end
[m,n]=size(ss);
sss = zeros(m,n);
for i=1:n
    col=ss(:,i);
    minc = min(col);
    maxc = max(col);
    sss(:,i) = (col-minc)/(maxc-minc);
end

if istensor
A_norm=reshape(sss,N,M,L);
else
   A_norm=sss; 
end
end

