function y_func_permutation(inputDir, outputDir, gsSize, nullType)
% Input
%   inputDir -- path to the input .mat file
%               This file should include following fields: mDataGEctx, staIMG
%   outputDir -- path to the output file
%   gsSize -- size of random gene set
%   nullType -- type of null models, "null-random" or "null-brain"
% by Yongbin Wei 2020, VU University Amsterdam


% load data
load(inputDir, 'mDataGEctx', 'staIMG', 'BRAINgene_idx');
nGenes = size(mDataGEctx, 2);

if isequal(nullType, 'nullrandom')
    idx_background = 1:nGenes;
    % output file name
    fname = fullfile(outputDir, ['gs_', num2str(gsSize),'.json']);
elseif isequal(nullType, 'nullbrain')
    idx_background = BRAINgene_idx;
    % output file name
    fname = fullfile(outputDir, ['brain_gs_', num2str(gsSize),'.json']);    
end

if ~isnumeric(gsSize)
    gsSize = str2num(gsSize);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% number of permutations
N = 10000;

% set const
C = ones(size(staIMG, 1), 1); 

% number of background index
n_bg = numel(idx_background);

% initialize beta
BETA = nan(N, size(staIMG, 2));

% permutation analysis
for n = 1: N
    rid = idx_background(randperm(n_bg, gsSize));    
    GG = nanmean(mDataGEctx(:, rid), 2);
    GG = (GG - nanmean(GG)) ./ nanstd(GG);
    
    % for each imaging phenotype
    for pheno_id = 1: size(staIMG, 2)            
        stat = regress(staIMG(:, pheno_id), ...
            [C, GG]);            
        BETA(n, pheno_id) = stat(2);            
    end
end
 
% mean and std
mBeta = nanmean(BETA, 1);
stdBeta = nanstd(BETA, '', 1);

% save results to json
datastruct = struct();

val  = round(mBeta, 5);
val = fillmissing(val, 'constant', 0);
datastruct.mbeta = val;
    
val  = round(stdBeta, 5);
val = fillmissing(val, 'constant', 0);
datastruct.stdbeta = val;

str = jsonencode(datastruct);        
fid = fopen(fname,'Wb');
fprintf(fid,'%s',str);
fclose(fid);

end