function [EX,mult] = JointCubes( X_blocks,Y,paraCube,Lambda,mu )
%函数作用：输入块blocks，Y原始数据，参数λ，关于块的参数params==>输出（与原始数据Y一起）平均后的EX
%其中，λ控制Y的影响大小。即实现Image denoising via learned dictionaries and sparse
%representation中的公式6.2
%Lambda表示参数λ；X是blocks，所有的块，(8*8*39*552)
%params是块的参数，包括block_num 块的个数;overlap_sz 块重叠大小;block_sz 块大小
%mu 表示拉普拉斯算子的系数  如果没有梯度项，mu=0；

block_num= paraCube.block_num;
Is_addrow= paraCube.Is_addrow;
Is_addcol= paraCube.Is_addcol;
patchsize= paraCube.patchsize;
overlap= paraCube.overlap;

sz=size(Y); 
number=reshape(1:(sz(1)*sz(2)),sz(1),sz(2));
range=find(~Y(:,:,1));%边界 0 的位置index
mult = zeros(size(Y));%mult记录元素叠加次数
EX0  = zeros(size(Y));
idx=0;
RTR=zeros(sz(1)*sz(2));
for j=1:block_num(2)
    for i= 1:block_num(1)
        ii = 1 + (i - 1)*(patchsize-overlap);%ii和jj是cube左上角的位置
        jj = 1 + (j - 1)*(patchsize-overlap);
        
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
            % blocks(:,:,:,idx)=img(ii:ii+patchsize-1, jj:jj+patchsize-1, :);
            mult(ii:ii+patchsize-1, jj:jj+patchsize-1, :)...
                =mult(ii:ii+patchsize-1, jj:jj+patchsize-1, :) + 1;
            EX0(ii:ii+patchsize-1, jj:jj+patchsize-1, :)...
                =EX0(ii:ii+patchsize-1, jj:jj+patchsize-1,:) +X_blocks(:,:,:,idx);
%             Mask=zeros(sz(1),sz(2));
%             Mask(ii:ii+patchsize-1, jj:jj+patchsize-1)=1;
%             R=diag(Mask(:));% 位置（i，j）patch的提取矩阵R
%             RTR=RTR+(R')*R;
        end
    end
end

if mu == 0    
    EX= (EX0+Lambda*Y)./(mult+Lambda);
else    
    %-----------add laplas-----------
  if patchsize==16
      load('D:\01Code\matlab\RTR.mat');
  elseif patchsize==10
    load('D:\01Code\matlab\RTR_10.mat');
  else
          disp("patchsize Error")
  end
  
    load('D:\01Code\matlab\Laplace matrix.mat');
    filter = [0 -1 0;  -1 4 -1;  0 -1 0];
    for k=1:size(Y,3)
        y=Y(:,:,k);
        Ly(:,:,k)=conv2(y,filter,'same');
    end
    
    temp_A=RTR+Lambda*eye(size(RTR,2))+mu*LapL;
    temp_B=(EX0+Lambda*Y+mu*Ly);
    for bb=1:size(temp_B,3)
        EXX=temp_A\reshape(temp_B(:,:,bb),[],1);  % inv(A)*b=A\b
        EX(:,:,bb)=reshape(EXX,sz(1),sz(2));
    end
end



