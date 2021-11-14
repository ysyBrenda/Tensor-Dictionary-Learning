function  [A ,DK]  =  My_sparse_tucker( SpatialDictionary1,SpatialDictionary2,D3,X,par )
%两个空间字典+1个元素字典：SpatialDictionary1,SpatialDictionary2,D3
%其中元素字典为1个矩阵，空间字典都是张量。
K=size(SpatialDictionary1,3);
%求解大字典DK.
for k=1:K
    D21=kron(SpatialDictionary2(:,:,k),SpatialDictionary1(:,:,k));   %SIZE=256*16
    z(:,:,k)=kron(D3(k,:),D21);
end
DK=reshape(permute(z,[1,3,2]),[],size(z,2));  %9728*32
%求DTZ
ZZ=reshape(X,[],size(X,4));  %9728*97
DKTZ=(DK')*ZZ;
%求逆
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
    vec=vec(:);%拉直成向量
    ccc=inv_DTDI*reshape(vec,[],size(X,4));  
    S  =bbb+ccc;%将括号拆开计算
    S=reshape(S,size(A));

    A         =   soft(S+V/(2*mu), par/(2*mu));
   V         =   V + 2*mu*( S - A );
end
% for  i  =  1:T
%     vec=2*mu*(A-V/(2*mu));
%     vec=vec(:);%拉直成向量
%     ccc=inv_DTDI*reshape(vec,[],97);  
%     S  =bbb+ccc;%将括号拆开计算
%     S=reshape(S,size(A));
% 
%     A         =   soft(S+V/(2*mu), par/(2*mu));
%    V         =   V + mu*( S - A );
% end



