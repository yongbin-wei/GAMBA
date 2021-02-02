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

imgIn = fullfile(dataPath, 'img_data.mat');
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
save('./input/input.mat', 'mDataGEctx', 'mDataGEall', 'staIMG',...
    'DATA_IMG_IDX', 'DATA_IMG_DESCRIPTIONS', 'BRAINgene_idx', ...
    'gene_symbols', 'regionDescriptionAll', 'regionDescriptionCtx');

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
outputDir = fullfile(outputPath, 'null_coexpression');

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
      y_func_generate_null_coexp_index('./input/input.mat', outDir, ...
          num2str(gs_all_size(ii)));
    end
end

disp('## ## Association with imaging');

load('./output/idx_null_coexp/idx_rand_10.mat', 'gs_all_size', 'targ_coexp');
m = numel(gs_all_size);

for ii = 1: m
    gsSize = gs_all_size(ii);
    coexpFile = ['./output/idx_null_coexp/idx_rand_', num2str(gsSize), '.mat'];
    if mod(gsSize, 10) == 0
        disp(['#### ', num2str(gsSize)]);
    end
    y_func_permutation_coexp(inputDir, outputDir, gsSize, coexpFile);
    y_func_permutation_ge_coexp(inputDir, outputDir, gsSize, coexpFile);
end


%% Null-spin model
disp('## Start null-spin model ...')

inputDir = './input/input.mat';
outputDir = fullfile(outputPath, 'null_spin');
spinDir = '../data/gene_expression_spinned';

poolobj = parpool(N_par);
parfor ii = 1: nGenes
    gene = gene_symbols{ii};
    disp(['#### ', num2str(ii), ': ', gene]);    

    y_func_permutation_spin(inputDir, spinDir, gene, outputDir);
end
delete(poolobj)

