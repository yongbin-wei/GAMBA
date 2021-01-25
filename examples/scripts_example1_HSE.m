clc, clear, close

projectPath = fileparts(fileparts(mfilename('fullpath')));
figurePath = fullfile(projectPath, 'figures');

addpath(genpath(fullfile(projectPath, 'misc')));

% load supragranular genes
tbl = readtable(fullfile(projectPath, 'examples', 'gene_supragranular.txt'), ...
    'Delimiter', '\t', 'ReadVariableNames', false);
gene_set_hse = unique(tbl.Var1);

% load brain genes
data = load(fullfile(projectPath, 'data', 'GTEx_brain_genes_0.05_updated.mat'));
gene_set_brain = data.gene_brain;

% load gene expression
ge = load(fullfile(projectPath, 'data', 'gene_expression.mat'));
regionDescription = ge.regionDescriptionCtx;

% load cov
IMG = load(fullfile(projectPath, 'data', 'IMG_DATA_ALL_SINGLE_GENE_20200803.mat'));
dataIMG = IMG.staIMG;

disp(['# HSE genes: ', num2str(numel(gene_set_hse))]);
disp(['# BRAIN genes: ', num2str(numel(gene_set_brain))]);

geneset = intersect(gene_set_hse, ge.gene_symbols);
[~, idx_gs] = ismember(geneset, ge.gene_symbols);
ngenes = numel(idx_gs);


% correlation with connectome metrics (LR)
GG = nanmean(ge.mDataGEctx(:, idx_gs), 2);
x = GG./ nanstd(GG);

% get connectome items
II = IMG.DATA_IMG_IDX == 5;
imgDATA = IMG.staIMG(:, II);

% linear regression
beta = nan(size(imgDATA,2), 1);
pval = nan(size(imgDATA,2), 1);
for ii = 1:size(imgDATA,2)
    y = imgDATA(:, ii);
    reg = regstats(y, x, 'linear', 'tstat');
    beta(ii,1) = reg.tstat.beta(2);
    pval(ii,1) = reg.tstat.pval(2);
end

% FDR
pval_adj = mafdr(pval, 'BHFDR', true);

disp('## Table 1. Correlation with connectome metrics:')
disp(table(IMG.DATA_IMG_DESCRIPTIONS(II), beta, pval, pval_adj))


% bar plot
sig_index = find(pval_adj < 0.05);
y_label = {'Nodal strength (NOS)'; 'Nodal strength (FA)'; 'Nodal strength (SD)';...
    'Nodal degree'; 'Nodal strength (FC)'};
y_barplot(beta, sig_index, ...
    'standardized \beta', y_label, ...
    fullfile(figurePath, 'HSE_connectome_chart.svg'), 1, 7, 4)


% scatter plot
y = imgDATA(:,1);
y_scatter(x, y, 'Gene expression', 'NOS', ...
    '../figures/HSE_nos.svg')

y = imgDATA(:,3);
y_scatter(x, y, 'Gene expression', 'SD', ...
    '../figures/HSE_svd.svg')

y = imgDATA(:,5);
y_scatter(x, y, 'Gene expression', 'FC', ...
    '../figures/HSE_fc.svg')
