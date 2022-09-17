clear all
clear memory;
clc
addpath('./tools');
addpath('./data_rec');
max_iter = 50;          % Maximum number of iterations
rho = 1.2;
miu = 1.0000e-02;
d = 100;                % dimension of latent space
max_miu = 1e6; 
load('scene3.mat');  	
c = length(unique(Y)); 
gt =Y;

for i = 1:length(X)
    X{i} = NormalizeData(X{i}');
end

LM_ini = cell(1,length(X));
Z_ini = cell(1,length(X));
for i = 1:length(X)   
   % ---------- initilization for Z  -------- %     
    options.NeighborMode = 'KNN';
    options.k = 15;
    options.WeightMode = 'Binary';      % Binary  HeatKernel
    Z= constructW(X{i}',options);
    Z = full(Z);
    Z1 = Z-diag(diag(Z));             
    Z = (Z1+Z1')/2;
    DZ= diag(sum(Z));
    LM_ini{i} = DZ - Z;      
    Z_ini{i} = Z;   
end
clear Z DZ


tic
         lambda1 = 100;
         lambda2 = 1;
         lambda3 = 0.01;                      
         lambda4 = 0.001;                         
         [Z] = ConsensusAffinity(X,Z_ini,LM_ini,d,lambda1,lambda2,lambda3,lambda4,max_iter,rho,miu, max_miu);   

         Zt = sparse(length(Y));
         for j=1:length(Z)
               Zt = Zt  + abs(Z{j}) + abs(Z{j}');
         end
         SS = Zt/length(Z);      

         result = clustering(SS, c, gt) 
        