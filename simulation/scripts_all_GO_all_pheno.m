function scripts_all_GO_all_pheno(expressionData, imgData, goData, spinDir, outputDir)
% compute associations between gene expression profiles of GO gene set and all phenotypes involved in GAMBA

% load processed gene expression data
ge = load(expressionData);
II_ctx = contains(ge.regionDescriptionCtx, 'ctx-lh-'); % only lh
regionDescriptions = ge.regionDescriptionCtx(II_ctx);
dataGE = ge.mDataGEctx(II_ctx, :);

% load phenotypic imaging data
IMG = load(imgData);
dataIMG = IMG.DATA_IMG;

% load GO terms
go = load(goData);
GOgs = cell(numel(go.GOnames), 1);
for ii = 1: numel(GOgs)
    gs = go.GOgenes(ii, :);
    II = cellfun(@(X) isempty(X), gs);
    gs = gs(~II);
    GOgs{ii, 1} = gs;
end

% initialize
nGO = numel(GOgs);
nImg = size(dataIMG, 2);
pval = nan(nGO, nImg);
beta = nan(nGO, nImg);
pval_nullspin = nan(nGO, nImg);
pval_nullrandom = nan(nGO, nImg);
pval_nullbrain = nan(nGO, nImg);

% loop all genes
for ii = 1: nGO
    disp(ii);
    
    gs = GOgs{ii,1};
    [~, idx_gs] = ismember(gs, ge.gene_symbol);
    gs(idx_gs == 0) = '';
    idx_gs(idx_gs == 0) = [];
    ngs = numel(idx_gs); % the actual number of genes

    % mean gene expression within the gs
    G1 = nanmean(dataGE(:, idx_gs), 2);
    x = ( G1 - nanmean(G1) ) ./ nanstd(G1);

    if nnz(isnan(x)) <= 5
        % linear regression
        for jj = 1:size(dataIMG, 2)
            y = dataIMG(:, jj);
            reg = regstats(y, x, 'linear', 'tstat');
            beta(ii, jj) = reg.tstat.beta(2);
            pval(ii, jj) = reg.tstat.pval(2);
        end 

        % null-spin
        for jj = 1:ngs % for each randomization
            spinfile = fullfile(spinDir,...
                ['GE_spin_', gs{jj}, '.txt']);
            ge_spin = [];
            if exist(spinfile, 'file')
                ge_spin(:,:,jj) = dlmread(spinfile);
            end
        end
        ge_spin = nanmean(ge_spin, 3); % 1000 randomizations by 57 regions
        
        for jj = 1:size(dataIMG, 2) % 1:384
            for kk = 1:size(ge_spin, 1) % 1:1000
                y = dataIMG(:, jj);
                x = ge_spin(kk, :)';
                x = (x - nanmean(x)) ./ nanstd(x);
                reg = regstats(y, x, 'linear', 'tstat');
                beta_spin(kk, jj) = reg.tstat.beta(2);
                pval_spin(kk, jj) = reg.tstat.pval(2);
            end
        end
        mbeta = nanmean(beta_spin, 1);
        stdbeta = nanstd(beta_spin, '', 1);

        zval_nullspin = (beta(ii, :) - mbeta) ./ stdbeta;
        pval_nullspin(ii, :) = 2 * (1 - normcdf(abs(zval_nullspin)));
    end
end

save([outputDir, '/LR_NULLSPIN_all_GO_all_pheno.mat'], 'beta', 'pval', ...
    'pval_nullspin', 'regionDescriptions');


% null-random
disp('## null-random');
idx_pool = 1:numel(ge.gene_symbol);

for ii = 1:nGO
    disp(ii);
    
    gs = GOgs{ii,1};
    [~, idx_gs] = ismember(gs, ge.gene_symbol);
    gs(idx_gs == 0) = '';
    idx_gs(idx_gs == 0) = [];
    ngs = numel(idx_gs); % the actual number of genes

    beta_tmp = nan(1000, size(dataIMG, 2));
    pval_tmp = nan(1000, size(dataIMG, 2));
    for jj = 1:1000
        ridx = idx_pool(randperm(numel(idx_pool), ngs));
        
        % mean gene expression within the gs
        G1tmp = nanmean(dataGE(:, ridx), 2);
        x = ( G1tmp - nanmean(G1tmp) ) ./ nanstd(G1tmp);
        if nnz(isnan(x)) <= 5
            % linear regression
            for kk = 1:size(dataIMG, 2)
                y = dataIMG(:, kk);
                reg = regstats(y, x, 'linear', 'tstat');
                beta_tmp(jj, kk) = reg.tstat.beta(2);
                pval_tmp(jj, kk) = reg.tstat.pval(2);
            end 
        end
    end
    mbeta_nullrandom = nanmean(beta_tmp, 1);
    stdbeta_nullrandom = nanstd(beta_tmp, '', 1);
    zval_nullrandom = (beta(ii, :) - mbeta_nullrandom) ./ stdbeta_nullrandom;
    pval_nullrandom(ii, :) = 2 * (1 - normcdf(abs(zval_nullrandom)));
end
save([outputDir, '/LR_NULLSPIN_all_GO_all_pheno.mat'], 'pval_nullrandom', '-append');


% null-brain
disp('## null-brain');
idx_pool = ge.BRAINgene_idx;

for ii = 1:nGO
    disp(ii);
    
    gs = GOgs{ii,1};
    [~, idx_gs] = ismember(gs, ge.gene_symbol);
    gs(idx_gs == 0) = '';
    idx_gs(idx_gs == 0) = [];
    ngs = numel(idx_gs); % the actual number of genes

    beta_tmp = nan(1000, size(dataIMG, 2));
    pval_tmp = nan(1000, size(dataIMG, 2));
    for jj = 1:1000
        ridx = idx_pool(randperm(numel(idx_pool), ngs));
        
        % mean gene expression within the gs
        G1tmp = nanmean(dataGE(:, ridx), 2);
        x = ( G1tmp - nanmean(G1tmp) ) ./ nanstd(G1tmp);
        if nnz(isnan(x)) <= 5
            % linear regression
            for kk = 1:size(dataIMG, 2)
                y = dataIMG(:, kk);
                reg = regstats(y, x, 'linear', 'tstat');
                beta_tmp(jj, kk) = reg.tstat.beta(2);
                pval_tmp(jj, kk) = reg.tstat.pval(2);
            end 
        end
    end
    mbeta_nullbrain = nanmean(beta_tmp, 1);
    stdbeta_nullbrain = nanstd(beta_tmp, '', 1);
    zval_nullbrain = (beta(ii, :) - mbeta_nullbrain) ./ stdbeta_nullbrain;
    pval_nullbrain(ii, :) = 2 * (1 - normcdf(abs(zval_nullbrain)));
end
save([outputDir, '/LR_NULLSPIN_all_GO_all_pheno.mat'], 'pval_nullbrain', '-append');

end


