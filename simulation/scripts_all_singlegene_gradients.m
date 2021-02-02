function scripts_all_singlegene_gradients(expressionData, gradientsData, spinDir, outputDir)
% compute associations between single gene expression profile and three
% gradients and their combinations

% load processed gene expression data
ge = load(expressionData);
II_ctx = contains(ge.regionDescriptionCtx, 'ctx-lh-'); % only lh
regionDescriptions = ge.regionDescriptionCtx(II_ctx);
dataGE = ge.mDataGEctx(II_ctx, :);

% load gradients
grad = load(gradientsData);
grad.regionProperties(:, 4) = (grad.regionProperties(:, 1) + ...
    grad.regionProperties(:, 2))./2;
grad.regionProperties(:, 5) = (grad.regionProperties(:, 3) + ...
    grad.regionProperties(:, 2))./2;
grad.regionProperties(:, 6) = (grad.regionProperties(:, 3) + ...
    grad.regionProperties(:, 1))./2;
grad.regionProperties(:, 7) = (grad.regionProperties(:, 2) + ...
    grad.regionProperties(:, 3) + grad.regionProperties(:, 1))./3;
dataIMG = grad.regionProperties;

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

        % null-spin
        spinfile = fullfile(spinDir, ...
            ['GE_spin_', ge.gene_symbol{ii}, '.txt']);

        if exist(spinfile, 'file')
            ge_spin = dlmread(spinfile);

            if nnz(isnan(ge_spin(1,:))) <= 5
                for jj = 1:size(dataIMG, 2)
                    for kk = 1:size(ge_spin, 1)
                        y = dataIMG(:, jj);
                        x = ge_spin(kk, :)';
                        x = (x - nanmean(x)) ./ nanstd(x);
                        reg = regstats(y, x, 'linear', 'tstat');
                        beta_spin(kk, jj) = reg.tstat.beta(2);
                        pval_spin(kk, jj) = reg.tstat.pval(2);
                    end
                end
            end

            mbeta = nanmean(beta_spin, 1);
            stdbeta = nanstd(beta_spin, '', 1);

            zval_nullspin = (beta(ii, :) - mbeta) ./ stdbeta;
            pval_nullspin(ii, :) = 2 * (1 - normcdf(abs(zval_nullspin)));
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
J = ge.BRAINgene_idx;
mbeta_nullbrain = nanmean(beta(J, :), 1);
stdbeta_nullbrain = nanstd(beta(J, :), '', 1);
for ii = 1:size(beta, 1)
    zval_nullbrain = (beta(ii, :) - mbeta_nullbrain) ./ stdbeta_nullbrain;
    pval_nullbrain(ii, :) = 2 * (1 - normcdf(abs(zval_nullbrain)));
end
save([outputDir, '/LR_NULLSPIN_all_singlegene_gradients.mat'], 'beta', 'pval', ...
    'pval_nullspin', 'pval_nullrandom', 'pval_nullbrain', 'regionDescriptions');
end

