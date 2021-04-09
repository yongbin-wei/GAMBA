clc, clear, close all;

%% Example 1
% I have an imaging map (.nii file) and a gene-set. I want to test if the 
% imaging pattern correlates to the pattern of gene expression





%% Example 2
% I have an imaging data matrix, a gene expression data matrix, and a gene-
% set. I want to test if the pattern of imaging data correlates to the 
% pattern of gene expression.

load('data/example_HSE.mat');
imgDescriptions = {'NOS','FA','SD','Degree','FC'};

% null-coexpression model
res_nullcoexp = permutation_null_coexp(img_data, geneset);

% null-brain model
res_nullbrain = permutation_null_brain(img_data, geneset);

% null-spin model

%% Example 3
% I have a gene-set. I want to test in which brain regions the gene-set is 
% over-expressed.




%% Example 4
% I have a gene expression matrix and a gene-set. I want to test in which 
% brain regions the gene-set is over-expressed.




%% Example 5
% I have an imaging map (.nii file), and I want to look for the most
% correlated genes









