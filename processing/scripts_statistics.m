clear, clc, close all

disp('Preprocessing data for GAMBA');

% settings: data path
filePath = fileparts(mfilename('fullpath'));
disp('## The current path is:')
disp(filePath);

dataPath = fullfile(fileparts(filePath), 'data');
disp('## The data path is:')
disp(dataPath)

imgIn = fullfile(dataPath, 'IMG_DATA_ALL_20200803.mat');
disp('## The input file is: ');
disp(imgIn)

outputPath = fullfile(filePath, 'output');
disp('## The output file path is');
disp(outputPath)

NN = 1000; % set the max GOI size
disp('##The maximal gene set size is:')
disp(NN)

N_par = 2; % parallel

% load data
load(imgIn);
load(fullfile(dataPath, 'gene_expression.mat'));
save('./input/input.mat', 'mDataGEctx', 'mDataGEall', 'staIMG', 'BRAINgene_idx');

% gene set sizes
gs_all_size = [1:499, 500:10:NN];
m = numel(gs_all_size);


%% Null-random model
disp('## Start permutation analysis');

inputDir = './input/input.mat';

disp('## Start null-random model ...')
poolobj = parpool(N_par);
parfor ii = 1: m
    gsSize = gs_all_size(ii);
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation(inputDir, outputPath, gsSize, 'nullrandom');
    y_func_permutation_ge(inputDir, outputPath, gsSize, 'nullrandom');
end
delete(poolobj)


%% Null-brain model
inputDir = './input/input.mat';

disp('## Start null-brain model ...')
poolobj = parpool(N_par);
parfor ii = 1: m
    gsSize = gs_all_size(ii);
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation(inputDir, outputPath, gsSize, 'nullbrain');
    y_func_permutation_ge(inputDir, outputPath, gsSize, 'nullbrain');
end
delete(poolobj)


%% Null-coexpression model
inputDir = './input/input.mat';

disp('## Start null-coexpression');

disp('## ## Searching random gene sets');
% obtain/generate the pool of gene index of random gene sets
outDir = './output/idx_null_coexp/';
if ~exist(outDir, 'dir')
    mkdir(outDir);
    % this script generate all needed variables
    y_func_generate_null_coexp_index('./input/input.mat', outDir);
end

disp('## ## Association with imaging');

load('./output/idx_null_coexp/idx_rand_10.mat', 'gs_all_size', 'targ_coexp');
m = numel(gs_all_size);

for ii = 1: m
    gsSize = gs_all_size(ii);
    coexpDir = ['./output/idx_null_coexp/idx_rand_', num2str(gsSize), '.mat'];
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation_coexp(inputDir, outputPath, gsSize, coexpDir);
    y_func_permutation_ge_coexp(inputDir, outputPath, gsSize, coexpDir);
end
