function y_func_permutation_ge(inputFile, outputDir, gs_size, nullType)
% =========================================================================
% Input
%   inputFile -- path to the input .mat file
%               This file should include following fields: mDataGEctx, staIMG
%   outputDir -- path to the output folder
%   gsSize -- size of random gene set
%   nullType -- type of null models, "null-random" or "null-brain"
%
% by Yongbin Wei 2020, VU University Amsterdam
% =========================================================================

% load data
load(inputFile, 'mDataGEall', 'BRAINgene_idx');
nGenes = size(mDataGEall, 2);
nRegions = size(mDataGEall, 1);

if isequal(nullType, 'nullrandom')
    idx_background = 1:nGenes;
    % output file name
    fname = fullfile(outputDir, ['ge_gs_', num2str(gs_size),'.json']);
elseif isequal(nullType, 'nullbrain')
    idx_background = BRAINgene_idx;
    % output file name
    fname = fullfile(outputDir, ['brain_ge_gs_', num2str(gs_size),'.json']);
end

if ~isnumeric(gs_size)
    gs_size = str2num(gs_size);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% number of permutations
N = 10000;

% number of background index
n_bg = numel(idx_background);

% initialize
tmpGE = nan(nRegions, N);

for n = 1: N
    rid = idx_background(randperm(n_bg, gs_size));
    tmpGE(:, n) = nanmean(mDataGEall(:, rid), 2);
end

mGE = nanmean(tmpGE, 2);
stdGE = nanstd(tmpGE, '', 2);

% save results to json
datastruct = struct();

val  = round(mGE, 5);
val = fillmissing(val, 'constant', 0);
datastruct.mGE = val;

val  = round(stdGE, 5);
val = fillmissing(val, 'constant', 0);
datastruct.stdGE = val;

str = jsonencode(datastruct);
fid = fopen(fname,'Wb');
fprintf(fid,'%s',str);
fclose(fid);

end