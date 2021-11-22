function [EX,mult] = JointCubes( X_blocks,Y,paraCube,Lambda,mu )
% input X(ALL blocks)，Y(original data)，Lambda (parameter  λ),params(parametersabout blocks)
%==>output (together with the original data Y) averaged EX
%where λ controls the effect of Y.
%params:including block_num (number of blocks);overlap_sz (block overlap size);block_sz(block size)
%mu indicates the coefficients of the Laplace operator. If there is no gradient term, mu=0.

block_num= paraCube.block_num;
Is_addrow= paraCube.Is_addrow;
Is_addcol= paraCube.Is_addcol;
patchsize= paraCube.patchsize;
overlap= paraCube.overlap;

sz=size(Y); 
number=reshape(1:(sz(1)*sz(2)),sz(1),sz(2));
range=find(~Y(:,:,1));
mult = zeros(size(Y));
EX0  = zeros(size(Y));
idx=0;
RTR=zeros(sz(1)*sz(2));
for j=1:block_num(2)
    for i= 1:block_num(1)
        ii = 1 + (i - 1)*(patchsize-overlap);
        jj = 1 + (j - 1)*(patchsize-overlap);
        
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
            % blocks(:,:,:,idx)=img(ii:ii+patchsize-1, jj:jj+patchsize-1, :);
            mult(ii:ii+patchsize-1, jj:jj+patchsize-1, :)...
                =mult(ii:ii+patchsize-1, jj:jj+patchsize-1, :) + 1;
            EX0(ii:ii+patchsize-1, jj:jj+patchsize-1, :)...
                =EX0(ii:ii+patchsize-1, jj:jj+patchsize-1,:) +X_blocks(:,:,:,idx);
%             Mask=zeros(sz(1),sz(2));
%             Mask(ii:ii+patchsize-1, jj:jj+patchsize-1)=1;
%             R=diag(Mask(:));
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



