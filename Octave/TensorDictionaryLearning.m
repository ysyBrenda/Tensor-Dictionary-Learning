clearvars
%%------------load data--------------
addpath(genpath('function'))
load('.\data\geoc_data.mat')

Y=d_value*100;    %d_value is transformed geochemical data after ilr 
%set atoms of dictionary
para.W=10;
para.H=10;
para.S=2;
para.K=1; 
para.patchsize=16;
para.overlap=para.patchsize/2;
para.J_Iter=1;
para.Lambda=1; 
paraD.lambda=0.1; 
paraD.lambda1= 1e-2; paraD.lambda2=1e-2;paraD.lambda3=1e-2

%%------------partition to get patches--------------
[blocks,paraCube]= ExtractCubes(Y,para.patchsize,para.overlap);

if para.K==1  
    [cluster_kmeans]=ones(1,size(blocks,4));
else
    Y2=Unfold(blocks,size(blocks),4); 
    [cluster_kmeans]=fkmeans(Y2,para.K);  
end

blocks_upd=blocks;
% for J=1:para.J_Iter
X_blocks=[];
%elemental dictionary
BB=blocks_upd;
cc3=Unfold(BB,size(BB),3);
%              cc3=unique(cc3','rows');
%              cc3=cc3';
D = cc3;
[~ ,nn]=size(cc3);
par.lambda=paraD.lambda3;%=1e-5
par.K= para.S;  %atoms numebr
[D3, SparseCode3 ,fun3]= Nonnegative_DL( D, par );   

for mn=1:max(cluster_kmeans) 
    gg=find(cluster_kmeans==mn);
    BB=blocks_upd(:,:,:,gg);
    
    for J=1:size(Y,3)
        Bb=BB(:,:,J,:);
        
        %y-dictionary learning
        cc1=Unfold(Bb,size(Bb),1);
        cc1=unique(cc1','rows'); 
        cc1=cc1';
        D = cc1;
        [~ ,nn]=size(cc1);
        par.lambda=paraD.lambda1;%=1e-5
        par.K=min(para.W,nn); 
        [D1 B1 ,fun1(:,J)] =Nonnegative_DL( D, par );   
        
        %x-dictionary learning
        cc2=Unfold(Bb,size(Bb),2);
        cc2=unique(cc2','rows');
        cc2=cc2';
        D = cc2;
        [~,nn]=size(cc2);
        par.lambda=paraD.lambda2
        par.K=min(para.H,nn);  
        [D2 B2 ,fun2(:,J)] = Nonnegative_DL( D, par );
        SpatialDictionary1(:,:,J)=D1;
        SpatialDictionary2(:,:,J)=D2;
        SparseCode1(:,:,J)=B1;
        SparseCode2(:,:,J)=B2;
    end
    [CoreTensor,Dk]=sparse_tucker(SpatialDictionary1,SpatialDictionary2,D3,BB,paraD.lambda);
    x=Dk*reshape(CoreTensor,[],size(BB,4));  %vector form
    X_blocks(:,:,:,gg)=reshape(x,size(BB,1),size(BB,2),[],size(BB,4));
end


%%-------------reconstruction--------------
[EX,mult] = JointCubes(X_blocks,Y,paraCube,1,5); % 0.1 & 0


%%-------------calculate anomaly score(Err)--------------
E =Y-EX;
Err=sqrt(sqrt(sum(E.*E,3)));

%%--------------Evaluation--------------
Au_Err=griddata(X_Long,Y_Lat,Err,Au_x,Au_y); %Positive Sample


Auc=[];
N_Au=numel(Au_Err);  

%%--- AUC----
%NIndex=logical(Err);
yi=Err(NIndex); %Negative sample
y_random=[];
for i=1:100
 ##   rng(i);
    y_random(:,i)=yi(randperm(numel(yi),N_Au))';
    [Auc(i) Zauc(i)]=computeAUC(Au_Err,y_random(:,i));
end
mean(Auc)



%% draw anomaly score
% clf
% x1=[119.75 120.67];
% y1=[36.76 37.66];
% imagesc(x1,y1,Err);
% colormap jet
% hold on
% plot(Au_x,Au_y,'r^')
% axis xy
% colorbar
% T=title('Err');
% T.FontSize = 14;


%% draw
% clf
% figure(2)
% i=1;
% clims=[0,80];%10,60
% subplot('Position',[0.1 0.4 0.25 0.45]);
% imagesc(EX(:,:,i),clims);colorbar;colormap jet;axis xy;
% T=title('reconstructed data');
% T.FontSize = 14;
% subplot('Position',[0.4 0.4 0.25 0.45]);
% imagesc(Y(:,:,i),clims);colorbar;axis xy
% T=title('raw data');
% T.FontSize = 14;
% subplot('Position',[0.7 0.4 0.25 0.45]);
% imagesc(Y(:,:,i)-EX(:,:,i));colorbar;axis xy
% T=title('differences');
% T.FontSize = 14;

