function  [A ,DK]  =  sparse_tucker( SpatialDictionary1,SpatialDictionary2,D3,X,par )
%input：2 spatial dictionaries and 1 elemental dictionary：SpatialDictionary1,SpatialDictionary2,D3
% element dictionary is a matrix and the spatial dictionaries are tensor.
K=size(SpatialDictionary1,3);
%solve for dictioanry DK.
for k=1:K
    D21=kron(SpatialDictionary2(:,:,k),SpatialDictionary1(:,:,k));   %SIZE=256*16
    z(:,:,k)=kron(D3(k,:),D21);
end
DK=reshape(permute(z,[1,3,2]),[],size(z,2));  %9728*32
%solve for DTZ
ZZ=reshape(X,[],size(X,4));  %9728*97
DKTZ=(DK')*ZZ;
%inverse
num=size(DK,2);
DTD_I=(DK')*DK+0.01*eye(num);
inv_DTDI=inv(DTD_I);

A= zeros(size(SpatialDictionary1,2),size(SpatialDictionary2,2),size(D3,2),size(X,4));
V = zeros( size(A));
T =100;
mu=0.01;
bbb=DTD_I\DKTZ;
%bbb=inv_DTDI*DKTZ;
for  i  =  1:T
    vec=mu*(A-V/(2*mu));
    vec=vec(:);
    ccc=inv_DTDI*reshape(vec,[],size(X,4));  
    S  =bbb+ccc;
    S=reshape(S,size(A));

    A         =   soft(S+V/(2*mu), par/(2*mu));
   V         =   V + 2*mu*( S - A );
end



