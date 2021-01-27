clear, clc, close all

disp('Preprocessing data for GAMBA');

% data path
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
m = 2%numel(gs_all_size);


%% Null-random model
% association with imaging phenotypes
disp('## Start permutation analysis');

inputDir = './input/input.mat';

disp('## Start null-random model ...')
disp('## ## Association with imaging');

poolobj = parpool(N_par);
parfor ii = 1: m
    gs_size = gs_all_size(ii);
    if mod(gs_size, 10) == 0
        disp(['####', num2str(gs_size)]);
    end
    y_func_permutation(inputDir, outputPath, gs_size, 'nullrandom');
end
delete(poolobj)

% over-expression
disp('## ## Gene expression values');
for ii = 1: m
    gs_size = gs_all_size(ii);    
    if mod(gs_size, 10) == 0
        disp(['####', num2str(gs_size)]);
    end
    y_func_permutation_ge(inputDir, outputPath, gs_size, 'nullrandom');
end


%% Null-brain model
disp('## Start null-brain model ...')

inputDir = './input/input.mat';

disp('## ## Association with imaging');
poolobj = parpool(N_par);
parfor ii = 1: m
    gs_size = gs_all_size(ii);
    if mod(gs_size, 10) == 0
        disp(['####', num2str(gs_size)]);
    end
    y_func_permutation(inputDir, outputPath, gs_size, 'nullbrain');
end
delete(poolobj)

% over-expression
disp('## ## Gene expression values');
for ii = 1: m
    gs_size = gs_all_size(ii);    
    if mod(gs_size, 10) == 0
        disp(['####', num2str(gs_size)]);
    end
    y_func_permutation_ge(inputDir, outputPath, gs_size, 'nullbrain');
end

return



%% Null-coexpression model
disp('## Start null-coexpression');

% obtain/generate the pool of gene index of random gene sets
outDir = './output/idx_null_coexp/';
if ~exist(outDir, 'dir')
    mkdir(outDir);
    % this script generate all needed variables
    y_func_generate_null_coexp_index('./input/input.mat', outDir);
end




return





%% Random genes with the similar co-expression level
disp('## Loading files');
load(imgIn);
load(fullfile(matrixPath, 'gene_expression.mat'));

C = ones(57, 1); % set const

disp('## ## Association with imaging');

% loop target coexperession level
for ii = 1:numel(targ_coexp)

    % loop gene set size
    for jj = 5:5:1000
        tmp = num2str(round(targ_coexp(ii) * 100));
        fname = fullfile(outputPath, 'random_gene_coexp', ...
            ['coexp_gs_',tmp,'_',num2str(jj),'.json']);
        disp(fname);
        BETA = nan(nPerm, size(staIMG, 2));

        % 1000 permutations
        for n = 1:nPerm
            % select random genes
            tmp_status = false;
            while tmp_status ~= true
                [rid, coexp, tmp_status] = y_func_random_coexp_samesize(...
                    mDataGEctx, targ_coexp(ii), jj);
            end     
            disp(num2str(n));
            
            % standardize
            GG = nanmean(mDataGEctx(:, rid), 2);
            GG = (GG - nanmean(GG)) ./ nanstd(GG);

            % for each imaging phenotype
            for pheno_id = 1: size(staIMG, 2)            
                stat = regress(staIMG(:, pheno_id), ...
                    [C, GG, stacov]);            
                BETA(n, pheno_id) = stat(2);            
            end
        end

        mBeta = round(nanmean(BETA, 1), 4);
        stdBeta = round(nanstd(BETA, '', 1), 4);

        % save results to json
        datastruct = struct();
        val  = mBeta;
        val = fillmissing(val, 'constant', 0);
        datastruct.mbeta = val;
        val  = stdBeta;
        val = fillmissing(val, 'constant', 0);
        datastruct.stdbeta = val;

        str = jsonencode(datastruct);
        fid = fopen(fname, 'Wb');
        fprintf(fid,'%s',str);
        fclose(fid);
    end
end