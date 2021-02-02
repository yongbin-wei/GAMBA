function y_func_permutation_ge_coexp(inputFile, outputDir, gsSize, coexpDir)
% =========================================================================
% Input
%   inputFile -- path to the input .mat file
%               This file should include following fields: mDataGEctx, staIMG
%   outputDir -- path to the output folder
%   gsSize -- size of random gene set
%   coexpDir -- path to the folder containing processed random gene index
%
% by Yongbin Wei 2020, VU University Amsterdam
% =========================================================================


% load data
load(inputFile, 'mDataGEall');
nGenes = size(mDataGEall, 2);
nRegions = size(mDataGEall, 1);

load(coexpDir, 'targ_coexp', 'gene_idx');
nPerm = size(gene_idx, 2);

% each coexp level
for ii = 1:numel(targ_coexp)
    tmp = num2str(round(targ_coexp(ii) * 100));
    fname = fullfile(outputDir, ['coexp_ge_gs_', tmp, '_', num2str(gsSize),'.json']);

    % initialize
    tmpGE = nan(nRegions, nPerm);

    for n = 1: nPerm
        rid = squeeze(gene_idx(ii, n, :))';
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