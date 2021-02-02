function y_func_permutation_spin(inputFile, spinDir, gene, outputDir)
% =========================================================================
% Input
%   inputFile -- path to the input .mat file
%               This file should include following fields: staIMG
%   spinDir -- path to the processed spinned gene expression data (i.e.,
%   GE_spin_SYMBOL.txt)
%   outputDir -- path to the output folder
%   gene -- gene symbol
%
% by Yongbin Wei 2020, VU University Amsterdam
% =========================================================================

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% load data
load(inputFile, 'mDataGEctx', 'mDataGEall', 'staIMG');

spinFile = fullfile(spinDir, ['GE_spin_', gene, '.txt']);

if exist(spinFile, 'file')
    geSpin = dlmread(spinFile);

    if size(geSpin, 1) == 1000
        geSpin = geSpin';
    end
    
    % set const
    C = ones(size(staIMG, 1), 1); 
    N = size(geSpin, 2);

    % linear regression
    BETA = nan(N, size(staIMG, 2)); % 1000 by 387
    
    for n = 1: N
        GG = geSpin(:, n);
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
    
    outputData = [mBeta', stdBeta'];
    dlmwrite(fullfile(outputDir, [gene, '_beta.txt']), outputData);
else
    warning([spinFile, ' not exist!!'])
end





