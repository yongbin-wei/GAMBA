clc, clear, close

projectPath = fileparts(fileparts(mfilename('fullpath')));
figurePath = fullfile(projectPath, 'figures');

if ~exist(figurePath, 'dir')
    mkdir(figurePath)
end

addpath(genpath(fullfile(projectPath, 'misc')));

% load brain genes
data = load(fullfile(projectPath, 'data', 'GTEx_brain_genes_0.05_updated.mat'));
gene_set_brain = data.gene_brain;

% load gene expression
ge = load(fullfile(projectPath, 'data', 'gene_expression.mat'));
regionDescription = ge.regionDescriptionCtx;

% load img
IMG = load(fullfile(projectPath, 'data', 'IMG_DATA_ALL_20200803.mat'));
dataIMG = IMG.staIMG;

% apoe
[~, idx_gs] = ismember('APOE', ge.gene_symbols);
G = ge.mDataGEall(:, idx_gs);


%% linear regression
GG = IMG.staGE(:, idx_gs);
x = GG;

II = IMG.DATA_IMG_IDX == 7;
imgDATA = IMG.staIMG(:, II);

for ii = 1:size(imgDATA,2)
    y = imgDATA(:, ii);
    reg = regstats(y, x, 'linear', 'tstat');
    beta(ii,1) = reg.tstat.beta(2);
    pval(ii,1) = reg.tstat.pval(2);
end

pval_adj = mafdr(pval,'BHFDR', true);

mDataBMVBM_description = IMG.DATA_IMG_DESCRIPTIONS(II);
mDataBMVBM_description{21} =  'Semantic dementia';
mDataBMVBM_description{9} =  'Alzheimer''s';

table(mDataBMVBM_description(pval_adj<0.05), ...
    beta(pval_adj<0.05),pval(pval_adj<0.05),pval_adj(pval_adj<0.05))

for i = 1:numel(mDataBMVBM_description)
    mDataBMVBM_description{i}(1) = upper(mDataBMVBM_description{i}(1));
end

% plot bar
[tmp, idx] = sort(beta, 'descend');
vals = tmp(1:10);
pval_adj_tmp = pval_adj(idx(1:10));
sig_index = find(pval_adj_tmp < 0.05);
y_label = mDataBMVBM_description(idx(1:10));
y_label{8} = 'Semantic dementia';

y_barplot(vals, sig_index, 'standardized \beta', y_label, ...
    fullfile(figurePath, 'APOE_vbm_chart.svg'), 1, 7, 4)

% scatter plot
y = imgDATA(:,1);
y_scatter(x, y, 'Gene expression', 'mALE',  fullfile(figurePath, 'apoe_adhd.svg'));

y = imgDATA(:,4);
y_scatter(x, y, 'Gene expression', 'mALE',  fullfile(figurePath, 'apoe_ftd.svg'));

y = imgDATA(:,9);
y_scatter(x, y, 'Gene expression', 'mALE',  fullfile(figurePath, 'apoe_ad.svg'));

y = imgDATA(:,13);
y_scatter(x, y, 'Gene expression', 'mALE',  fullfile(figurePath, 'apoe_dementia.svg'));

y = imgDATA(:,21);
y_scatter(x, y, 'Gene expression', 'mALE',  fullfile(figurePath, 'apoe_sdementia.svg'));


%% null-spin
datapath = '../processing/output/genes/';
jasonfile = fullfile(datapath, 'APOE.json');
fid = fopen(jasonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

II = IMG.DATA_IMG_IDX == 7;
mbeta = val.spin_mbeta(II,1);
stdbeta = val.spin_sbeta(II,1);
zval = (beta - mbeta) ./ stdbeta;
pval = 2* (1 - normcdf(abs(zval)));

pval_adj = mafdr(pval,'BHFDR', true);

disp(table(mDataBMVBM_description(find(pval<0.05)), ...
    zval(pval< 0.05), pval(pval< 0.05), pval_adj(pval< 0.05)));

% plot bar
[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = mDataBMVBM_description(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'APOE_vbm_chart_nspin.svg'), 1);


%% null random
jasonfile = fullfile('../processing/output/null_random/gs_1.json');
fid = fopen(jasonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

mbeta = val.mbeta(II);
stdbeta = val.stdbeta(II);

zval = (beta - mbeta) ./ stdbeta;
pval = 2* (1 - normcdf(abs(zval)));
pval_adj = mafdr(pval,'BHFDR', true);
id_sig = find(pval < 0.05);

disp(table(mDataBMVBM_description(find(pval<0.05)), ...
    zval(pval< 0.05), pval(pval< 0.05), pval_adj(pval< 0.05)));

% plot bar
[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = mDataBMVBM_description(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'APOE_vbm_chart_nrand.svg'), 1)


%% null brain
jasonfile = fullfile('../processing/output/null_brain/brain_gs_1.json');
fid = fopen(jasonfile); 
raw = fread(fid, inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

mbeta = val.mbeta(II);
stdbeta = val.stdbeta(II);

zval = (beta - mbeta) ./ stdbeta;
pval = 2* (1 - normcdf(abs(zval)));

pval_adj = mafdr(pval,'BHFDR', true);

disp(table(mDataBMVBM_description(find(pval<0.2)), ...
    zval(pval< 0.2), pval(pval< 0.2), pval_adj(pval< 0.2)));

id_sig = find(pval<0.5);

% plot bar
[tmp, idx] = sort(zval, 'descend');
vals = tmp(1:10);
pval_tmp = pval(idx(1:10));
sig_index = find(pval_tmp < 0.05);
y_label = mDataBMVBM_description(idx(1:10));
y_barplot(vals, sig_index, 'z score', y_label, ...
    fullfile(figurePath, 'APOE_vbm_chart_nbrain.svg'), 1)

