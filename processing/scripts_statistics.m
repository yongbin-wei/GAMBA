% Performing statistical analyses of (1) null-random, (2) null-brain, and
% (3) null-coexpression models

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
outputDir = fullfile(outputPath, 'null_random');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

disp('## Start null-random model ...')
poolobj = parpool(N_par);
parfor ii = 1: m
    gsSize = gs_all_size(ii);
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation(inputDir, outputDir, gsSize, 'nullrandom');
    y_func_permutation_ge(inputDir, outputDir, gsSize, 'nullrandom');
end
delete(poolobj)


%% Null-brain model
inputDir = './input/input.mat';
outputDir = fullfile(outputPath, 'null_brain');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

disp('## Start null-brain model ...')
poolobj = parpool(N_par);
parfor ii = 1: m
    gsSize = gs_all_size(ii);
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation(inputDir, outputDir, gsSize, 'nullbrain');
    y_func_permutation_ge(inputDir, outputDir, gsSize, 'nullbrain');
end
delete(poolobj)


%% Null-coexpression model
inputDir = './input/input.mat';
outputDir = fullfile(outputPath, 'null_random');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

disp('## Start null-coexpression');

disp('## ## Searching random gene sets');
% obtain/generate the pool of gene index of random gene sets
outDir = './output/idx_null_coexp/';
if ~exist(outDir, 'dir')
    mkdir(outDir);
    % this script generate all needed variables
    gs_all_size = [5:5:495, 500:10:1000];
    m = numel(gs_all_size);
    for ii = 1:m
      y_func_generate_null_coexp_index('./input/input.mat', outDir, num2str(gs_all_size(ii)));
    end
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
    y_func_permutation_coexp(inputDir, outputDir, gsSize, coexpDir);
    y_func_permutation_ge_coexp(inputDir, outputDir, gsSize, coexpDir);
end

