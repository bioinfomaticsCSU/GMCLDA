% GMCLDA: Prediction of lncRNA-disease associations based on inductive matrix completion
%  
%%%%%%%%%%
%   lncdis_demo.mat: an m*n association matrix between lncRNAs and diseases, m is 
%the number of lncRNAs, and n is the number of diseases
%   lncseqSim_demo.mat: an m*m sequence similarity matrix of lncRNAs
%   disSim_demo.mat: an n*n similarity matrix of disease
%% load data
LD=importdata('lncdis_demo.mat');
lncSim=importdata('lncseqSim_demo.mat');
dissim=importdata('disSim_demo.mat');
%% initialization
[nl,nd]=size(LD);
posIndices= find(LD==1);    % indices of positive samples
disp(['Total: ' num2str(length(posIndices)) ' interactions']);
%% set parameters
k=70;
weight=1; %
prob_params.gamma_n = 0.005;    % alpha_n
prob_params.gamma_l = 0.9;  % alpha_l
prob_params.gamma_d = 1;    % alpha_d
prob_params.size_X = size(LD);
%% GMCLDA geometric matrix completion
    
% complete interaction information for a new lncRNA
for i=1:nl
    if length(find(LD(i,:)))==0       
        rowVec=lncSim(i,:);                 
        rowVec(i)=0;        
        [simVal,simIndex]=sort(rowVec,'descend');        
        new_row=zeros(1,nd);
        for j=1:k
             new_row=new_row+weight^(j-1)*LD(simIndex(j),:);
        end
        new_row=new_row/k;
        LD(i,:)=new_row;
    end
end
    
% computing Gaussian interaction profile kernel of lncRNAs
[lncKernel,~]=gKernel(nl,nd,LD);

% computing Laplacian matrix
Wl=lncKernel;
Ll=Laplacian(Wl);

Wd=dissim;
Ld=Laplacian(Wd);

% mark observed indices
mask_train=false(numel(LD),1);
mask_train(posIndices)=true;   % indices of training samples
y_train=LD(mask_train); % training samples

mask_val=false(numel(LD),1);
y_val=LD(mask_val);

y_test=[];
mask_test=[];

% set problem parameters
prob_params.Ll=Ll;
prob_params.Ld=Ld;
prob_params.mask_val = mask_val;
prob_params.mask_test = mask_test;

prob_params.At_op = @(x) sample_sparse_t(x, mask_train);
prob_params.AtA_op = @(x) sample_sparse_AtA(x, mask_train);
    
% optimization using ADMM
[X_MC_graphs, stat_MC_graphs] = ADMM_PCG(y_train, y_val, y_test, prob_params);     
