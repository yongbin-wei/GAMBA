clc, clear, close

projectPath = fileparts(fileparts(mfilename('fullpath')));
figurePath = fullfile(projectPath, 'figures');

if ~exist(figurePath, 'dir')
    mkdir(figurePath)
end

addpath(genpath(fullfile(projectPath, 'misc')));

% load gene expression
ge = load(fullfile(projectPath, 'data', 'gene_expression.mat'));
regionDescription = ge.regionDescriptionCtx;

% load img
IMG = load(fullfile(projectPath, 'data', 'img_data.mat'));
dataIMG = IMG.staIMG;

% load gene sets
T = readtable('./asdTop25.txt', 'ReadVariableNames', false);
geneset = intersect(T.Var1, ge.gene_symbols);
[~, idx_gs] = ismember(geneset, ge.gene_symbols);


%% linear regression
GG = nanmean(ge.mDataGEctx(:, idx_gs), 2);
x = (GG - mean(GG))./ nanstd(GG);
data_bmfmri = IMG.staIMG(:, IMG.DATA_IMG_IDX == 6);
data_bmfmri_descriptions = ...
    IMG.DATA_IMG_DESCRIPTIONS(IMG.DATA_IMG_IDX == 6);

% linear regression
for ii = 1:size(data_bmfmri, 2)
    y = data_bmfmri(:, ii);
    reg = regstats(y, x, 'linear', 'tstat');
    beta(ii,1) = reg.tstat.beta(2);
    pval(ii,1) = reg.tstat.pval(2);
end

% FDR
pval_adj = mafdr(pval,'BHFDR', true);

table(beta(pval<0.05), pval(pval<0.05), pval_adj(pval<0.05),...
    data_bmfmri_descriptions(pval<0.05), ...
    'VariableNames', {'beta', 'pval', 'pval_adj', 'items'})
% plot bar
[tmp, idx] = sort(beta, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = data_bmfmri_descriptions(idx(1:10));
y_barplot(vals, sig_index, 'standardized \beta', y_label, ...
    fullfile(figurePath, 'ASD_fmri_chart.svg'), 1, 4.5, 4)

[~, j] = ismember('Asperger', data_bmfmri_descriptions);
y = data_bmfmri(:, j);
y_scatter(x, y, 'Gene expression', 'mALE', ...
    fullfile(figurePath,'ASD_asperger.svg'));


%% null-spin
spinPath = '../processing/output/genes/';
II = IMG.DATA_IMG_IDX == 6;
for ii = 1: numel(geneset)
    jsonfile = fullfile(spinPath, [geneset{ii},'.json']);
    fid = fopen(jsonfile); 
    raw = fread(fid, inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    mtmp(ii, :) = val.spin_mbeta(II); 
    stmp(ii, :) = val.spin_sbeta(II);   
end
mbeta = nanmean(mtmp, 1)';
stdbeta = nanmean(stmp, 1)';

zval = (beta - mbeta) ./ stdbeta;
pval = 2* (1 - normcdf(abs(zval)));
pval_adj = mafdr(pval,'BHFDR', true);

table(zval, pval, pval_adj, data_bmfmri_descriptions)

[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = data_bmfmri_descriptions(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'ASD_fmri_chart_nspin.svg'), 1, 4, 4)


%% null-random
nrandomPath = '../processing/output/null_random/';
jsonfile = fullfile(nrandomPath, ['gs_', num2str(numel(idx_gs)),'.json']);

fid = fopen(jsonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

% extract mean and std
mbeta = val.mbeta(IMG.DATA_IMG_IDX == 6);
stdbeta = val.stdbeta(IMG.DATA_IMG_IDX == 6);

% compute z score
zval = (beta - mbeta) ./ stdbeta;
pval = 2 * (1 - normcdf(abs(zval)));
pval_adj = mafdr(pval,'BHFDR', true);
table(zval, pval, pval_adj, data_bmfmri_descriptions)

% plot bar
[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = data_bmfmri_descriptions(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'ASD_fmri_chart_nrand.svg'), 1, 4, 4)


%% null-brain
nrandomPath = '../processing/output/null_brain/';
jsonfile = fullfile(nrandomPath, ['brain_gs_', num2str(numel(idx_gs)),'.json']);

fid = fopen(jsonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

mbeta = val.mbeta(IMG.DATA_IMG_IDX == 6);
stdbeta = val.stdbeta(IMG.DATA_IMG_IDX == 6);

zval = (beta - mbeta) ./ stdbeta;
pval = 2 * (1 - normcdf(abs(zval)));
pval_adj = mafdr(pval,'BHFDR', true);
table(zval, pval, pval_adj, data_bmfmri_descriptions)

[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = data_bmfmri_descriptions(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'ASD_fmri_chart_nbrain.svg'), 1, 4, 4)


%% null-coexpression
G = ge.mDataGEctx(:, idx_gs);
Rtmp = corr(G);
Rtmp = Rtmp - diag(diag(Rtmp));
Rtmp = triu(Rtmp);
orgCoexp = nanmean(nonzeros(Rtmp));

disp(['## Mean Coexpression level of HSE genes:', num2str(orgCoexp)]);

MM = round(numel(idx_gs)/5)*5;

NN = round(orgCoexp * 100);
NN = floor(NN./5).*5;

nCoexpPath = '../processing/output/null_coexpression/';
jsonfile = fullfile(nCoexpPath, ...
    ['coexp_gs_', num2str(NN), '_', num2str(MM), '.json'])

fid = fopen(jsonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

mbeta = val.mbeta(II);
stdbeta = val.stdbeta(II);
zval = (beta - mbeta) ./ stdbeta;
pval = 2 * (1 - normcdf(abs(zval)));
pval_adj = mafdr(pval,'BHFDR', true);
table(zval, pval, pval_adj, data_bmfmri_descriptions)

[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = data_bmfmri_descriptions(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'ASD_fmri_chart_ncoexp.svg'), 1, 4, 4)

