clc, clear, close all;

% Example of human-supragranular-enriched genes
load('data/example_HSE.mat');
imgDescriptions = {'NOS','FA','SD','Degree','FC'};

% null-coexpression model
res_nullcoexp = permutation_null_coexp(img_data, geneset);

% null-brain model
res_nullbrain = permutation_null_brain(img_data, geneset);

% null-spin model



