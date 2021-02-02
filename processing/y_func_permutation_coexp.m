function y_func_permutation_coexp(inputDir, outputDir, gsSize, coexpDir)
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


if ~isnumeric(gsSize)
    gsSize = str2num(gsSize);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% load data
load(inputDir, 'mDataGEctx', 'staIMG');
load(coexpDir, 'targ_coexp', 'gene_idx');
nPerm = size(gene_idx, 2);

% set const
C = ones(size(staIMG, 1), 1); 

% each coexp level
for ii = 1:numel(targ_coexp)
    tmp = num2str(round(targ_coexp(ii) * 100));
    fname = fullfile(outputDir, ['coexp_gs_', tmp, '_', num2str(gsSize),'.json']);
    
    BETA = nan(nPerm, size(staIMG, 2));
    
    for n = 1:nPerm
        rid = squeeze(gene_idx(ii, n, :))';
            
        % standardize
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

end