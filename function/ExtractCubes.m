function [blocks,paraCube]= ExtractCubes(Y,patchsize,overlap)
%EXTRACTCUBES Y是原始张量，patchsize是块的大小，overlap是重叠的个数；
%blocks是Y分块后的结果，paraCube中包含着block_num块的行列个数、是否单独为行/列采样的参数。
sz=size(Y);    %96 100 39
% 换步长采取cubes：
img=Y;
block_num = floor((sz(1:2) - overlap)./(patchsize - overlap));
Is_addrow=0;Is_addcol=0;%表示是否要单独采cube，如果为1则需要单独采cube
if mod((sz(1)-patchsize),(patchsize-overlap))~=0
    block_num(1)= block_num(1)+1;
    Is_addrow=1;
end
if mod((sz(2)-patchsize),(patchsize-overlap))~=0
    block_num(2)= block_num(2)+1;
    Is_addcol=1;
end
range=find(~img(:,:,1));%边界 0 的位置index
number=reshape(1:(sz(1)*sz(2)),sz(1),sz(2));
blocks=[];
%blocks=zeros([patchsize,patchsize, sz(3),1]);%blocks数量先设置小一点。
idx=0;
for j=1:block_num(2)  %一行有num（2）个cube
    for i= 1:block_num(1)
        ii = 1 + (i - 1)*(patchsize- overlap);%ii和jj是cube左上角的位置
        jj = 1 + (j - 1)*(patchsize- overlap);
        
        if Is_addrow==1 && i==block_num(1)     %如果要单独采cube并且到了边缘位置，则更新ii，jj的位置
            ii= sz(1)-patchsize+1;
        end
        if Is_addcol==1 && j==block_num(2)
            jj=sz(2)-patchsize+1;
        end
        
        is_out=ismember(number(ii:ii+patchsize-1, jj:jj+patchsize-1),range);%超出边界，isout返回1
        %is_out全为0才执行,取cube
        if ~sum(sum(is_out))
            %idx = (j-1)*block_num(1) + i;%第几个
            idx=idx+1;
            blocks(:,:,:,idx)=img(ii:ii+patchsize-1, jj:jj+patchsize-1, :);  %取cube的操作
        end
    end
end
paraCube.block_num=block_num;
paraCube.Is_addrow=Is_addrow;
paraCube.Is_addcol=Is_addcol;
paraCube.patchsize=patchsize;
paraCube.overlap=overlap;

end

