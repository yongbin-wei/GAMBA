function scripts_all_singlegene_all_pheno(expressionData, imgData, braingenesData, spinDir, outputDir)
% compute associations between single gene expression profile and all
% phenotypes involved in GAMBA

% load processed gene expression data
ge = load(expressionData);
II_ctx = contains(ge.regionDescriptions, 'ctx-lh-'); % only lh
regionDescriptions = ge.regionDescriptions(II_ctx);
dataGE = ge.mean_gene_expression(II_ctx, :);

% load phenotypic imaging data
IMG = load(imgData);
dataIMG = IMG.DATA_IMG;

% initialize
pval = nan(size(dataGE, 2), size(dataIMG, 2));
beta = nan(size(dataGE, 2), size(dataIMG, 2));
pval_nullspin = nan(size(dataGE, 2), size(dataIMG, 2));
pval_nullrandom = nan(size(dataGE, 2), size(dataIMG, 2));
pval_nullbrain = nan(size(dataGE, 2), size(dataIMG, 2));


% loop all genes
for ii = 1: size(dataGE, 2)
    disp(ii);
    x = ( dataGE(:, ii) - nanmean(dataGE(:, ii)) ) ./ nanstd(dataGE(:, ii));
    if nnz(isnan(x)) <= 5
        % linear regression
        for jj = 1:size(dataIMG, 2)
            y = dataIMG(:, jj);
            reg = regstats(y, x, 'linear', 'tstat');
            beta(ii, jj) = reg.tstat.beta(2);
            pval(ii, jj) = reg.tstat.pval(2);
        end 
    end
end


% null-random
mbeta_nullrandom = nanmean(beta, 1);
stdbeta_nullrandom = nanstd(beta, '', 1);
for ii = 1:size(beta, 1)
    zval_nullrandom = (beta(ii, :) - mbeta_nullrandom) ./ stdbeta_nullrandom;
    pval_nullrandom(ii, :) = 2 * (1 - normcdf(abs(zval_nullrandom)));
end


% null-brain
braingenes = load(braingenesData);
[~, J] = ismember(braingenes.gene_brain, ge.gene_symbol);
J(J==0) = [];
mbeta_nullbrain = nanmean(beta(J, :), 1);
stdbeta_nullbrain = nanstd(beta(J, :), '', 1);
for ii = 1:size(beta, 1)
    zval_nullbrain = (beta(ii, :) - mbeta_nullbrain) ./ stdbeta_nullbrain;
    pval_nullbrain(ii, :) = 2 * (1 - normcdf(abs(zval_nullbrain)));
end

save([outputDir, '/LR_all_singlegene_all_pheno.mat'], 'beta', 'pval', ...
    'pval_nullspin', 'pval_nullrandom', 'pval_nullbrain', 'regionDescriptions');

end


