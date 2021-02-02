% Make json file for GAMBA
clear, clc, close all;

filePath = fileparts(mfilename('fullpath'));
disp('## The current path is:')
disp(filePath);

% path to the web-application source 
outputPath = fullfile(filePath, 'output');
disp('## The output folder is: ');
disp(outputPath)

load('./input/input.mat')

%% write profiles of phenotypes
datastruct = struct();

val  = round(staIMG', 4);
% N x 57, so pheno[1] will be an array of the profile used in GLM

val = fillmissing(val, 'constant', 0);
datastruct.img = val;
    
% % % val  = round(stacov, 4);
% % % val = fillmissing(val, 'constant', 0);
% % % datastruct.cov = val;

str = jsonencode(datastruct);
fname = fullfile(outputPath, 'phenotypes.json');
fid = fopen(fname,'Wb');
fprintf(fid,'%s',str);
fclose(fid);


%% write phenotype index
datastruct = struct();

val  = DATA_IMG_IDX;
datastruct.idx = val;
str = jsonencode(datastruct);
   
fname = fullfile(outputPath, 'phenotypesIDX.json');
fid = fopen(fname, 'Wb');
fprintf(fid, '%s', str);
fclose(fid);


%% write pheno names
fname = fullfile(outputPath, 'phenotypeNames.json');
fid = fopen(fname, 'Wb');
fprintf(fid,'%s','[');
for i=1:numel(DATA_IMG_DESCRIPTIONS)
    fprintf(fid,strcat('"',DATA_IMG_DESCRIPTIONS{i},'"'));
    if i~=numel(DATA_IMG_DESCRIPTIONS)
        fprintf(fid,',');
    end
end
fprintf(fid,'%s',']');
fclose(fid);


%% write gene expression profiles

for n = 1:(numel(gene_symbols)-1)    
    gene = gene_symbols{n};
    disp(gene);
    
    % load results from Null-spin
    spinfile = ['../output/null_spin/', gene, '_beta.txt'];
    if exist(spinfile, 'file')
        spindata = dlmread(spinfile);
    else
        spindata = zeros(size(DATA_IMG,2), 2);
    end
    
    % gene expression all regions
    datastruct = struct();    
    val = round(mDataGEall(:, n)', 4);
    val = fillmissing(val, 'constant', 0);
    datastruct.ahba = val;

    % gene expression cortical regions
    val = round(mDataGEctx(:, n)', 4);
    val = fillmissing(val, 'constant', 0);
    datastruct.ahba_ctx = val;
    
    % mean/std beta from Null-spin model
    val = round(spindata(:, 1)', 4);
    val = fillmissing(val, 'constant', 0);
    datastruct.spin_mbeta = val;
    
    val = round(spindata(:, 2)', 4);
    val = fillmissing(val, 'constant', 0);
    datastruct.spin_sbeta = val;

    str = jsonencode(datastruct);

    if ~exist(fullfile(outputPath, 'genes'), 'dir')
        mkdir(fullfile(outputPath, 'genes'));
    end

    fname = fullfile(outputPath, 'genes', [gene,'.json']);
    fid = fopen(fname,'Wb');
    fprintf(fid,'%s',str);
    fclose(fid);
end

disp('## Write all gene expression to: ');
disp([outputPath, 'genes_v3']);


%% write phenotype data
table1 = table(tbl.features(DATA_IMG_IDX), DATA_IMG_DESCRIPTIONS);
table1.Properties.VariableNames = {'Category', 'Item'};

table2 = array2table(DATA_IMG');
tmp = regionDescriptionCtx;
tmp = cellfun(@(X) strrep(X, '-', '_'), tmp, 'UniformOutput', false);
table2.Properties.VariableNames = tmp;

table3 = [table1, table2];
table3(contains(table3.Category, 'cross_disorder'), :) = [];
head(table3)
writetable(table3, fullfile(outputPath, 'GAMBA_pheno_all.csv'))


%% write region names
fname = fullfile(outputPath,'regionNames.json');
datastruct = struct();
datastruct.regionDescriptions_all = regionDescriptionAll;
datastruct.regionDescriptions_ctx = regionDescriptionCtx;
str = jsonencode(datastruct);
fid = fopen(fname, 'Wb');
fprintf(fid, '%s', str);
fclose(fid);


%% write gene_symbol json
fname = fullfile(outputPath, 'GENEsymbols.json');
fid = fopen(fname,'Wb');
fprintf(fid,'%s','[');
for i=1:numel(gene_symbols)
    fprintf(fid,strcat('"',gene_symbols{i},'"'));
    if i~=numel(gene_symbols)
        fprintf(fid,',');
    end
end
fprintf(fid,'%s',']');
fclose(fid);


%% write colors_blue.json
addpath('~/Documents/codes/cbrewer/'); % use cbrewer to make color maps
CT = cbrewer('seq','Blues',101)*255;% Reds, Oranges, Greens, Purples, Greys
CT = flipud(CT);
CT1 = [dec2hex(CT(:,1)),dec2hex(CT(:,2)),dec2hex(CT(:,3))];
fname = fullfile(outputPath, 'colors', 'colors_blue.json');
fid = fopen(fname,'Wb');
fprintf(fid,'%s','[');
for i=1:size(CT1,1)
    fprintf(fid,strcat('"#',CT1(i,:),'"'));
    if i~=size(CT1,1)
        fprintf(fid,',');
    end
end
fprintf(fid,'%s',']');
fclose(fid);

