function [blocks,paraCube]= ExtractCubes(Y,patchsize,overlap)
%EXTRACTCUBES Y(original data)，patchsize(block size)，overlap(block overlap size)；
%blocks are patches of Y
sz=size(Y);    
% obtain patches 
img=Y;
block_num = floor((sz(1:2) - overlap)./(patchsize - overlap));
Is_addrow=0;Is_addcol=0;
if mod((sz(1)-patchsize),(patchsize-overlap))~=0
    block_num(1)= block_num(1)+1;
    Is_addrow=1;
end
if mod((sz(2)-patchsize),(patchsize-overlap))~=0
    block_num(2)= block_num(2)+1;
    Is_addcol=1;
end
range=find(~img(:,:,1));
number=reshape(1:(sz(1)*sz(2)),sz(1),sz(2));
blocks=[];
%blocks=zeros([patchsize,patchsize, sz(3),1]);
idx=0;
for j=1:block_num(2)  
    for i= 1:block_num(1)
        ii = 1 + (i - 1)*(patchsize- overlap);
        jj = 1 + (j - 1)*(patchsize- overlap);
        
        if Is_addrow==1 && i==block_num(1)     
            ii= sz(1)-patchsize+1;
        end
        if Is_addcol==1 && j==block_num(2)
            jj=sz(2)-patchsize+1;
        end
        
        is_out=ismember(number(ii:ii+patchsize-1, jj:jj+patchsize-1),range);
        
        if ~sum(sum(is_out))
            %idx = (j-1)*block_num(1) + i;
            idx=idx+1;
            blocks(:,:,:,idx)=img(ii:ii+patchsize-1, jj:jj+patchsize-1, :);  
        end
    end
end
paraCube.block_num=block_num;
paraCube.Is_addrow=Is_addrow;
paraCube.Is_addcol=Is_addcol;
paraCube.patchsize=patchsize;
paraCube.overlap=overlap;

end

