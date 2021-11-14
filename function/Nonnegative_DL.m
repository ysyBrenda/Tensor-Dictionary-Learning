% =========================================================================
% NSSR for Hyperspectral image super-resolution, Version 1.0
% Copyright(c) 2016 Weisheng Dong
% All Rights Reserved.
% 
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
% ----------------------------------------------------------------------
% 
% This is an implementation of the algorithm for Hyperspectral image super-
% resolution from a pair of low-resolution hyperspectral image and a high-
% resolution RGB image.
%  
% Please cite the following paper if you use this code:
% 
% Weisheng Dong, Fazuo Fu, et. al.,"Hyperspectral image super-resolution via 
% non-negative structured sparse representation", IEEE Trans. On Image Processing, 
% vol. 25, no. 5, pp. 2337-2352, May 2015.

%--------------------------------------------------------------------------
function    [D,B,fun]    =  Nonnegative_DL( X, par )
rng('default');%Add:用默认的随机数环境，每次生成的随机数相同。
Q            =    randperm( size(X,2) );%随机数列 以X的列数为参数
D            =    X(:, Q(1:par.K));%随机选取X的前K列
D            =    D./repmat((sqrt(sum(D.^2))+eps), size(D,1), 1);%标准化
Iter         =   10;  
X_s          =    X;
for  t   =  1 : Iter
    % Nonnegative Sparse coding         %非负稀疏编码
    [B,mu, fun(t)]    =    Nonnegative_SC( D, X_s, par );        
    b     =    sum(B.^2, 2);
    R     =    X - D*B;   
    % update dictionary
    for k = 1 : par.K       
        d_k_pre   =   D(:,k);
        d_k       =   max( D(:,k) + R*B(k,:)'/b(k), 0 );
        d_k       =   d_k./max(norm(d_k), 1);
        D(:,k)    =   d_k;
      
        R0        =   R;
        R         =   R0 - (d_k - d_k_pre)*B(k,:);
    end
% disp("mu:"+mu)   %ADD
% disp("fun:"+fun) %ADD
end

