clc, clear, close all

datapath = fullfile(fileparts(mfilename('fullpath')), 'output');

% load data
data = load(fullfile(datapath, 'LR_NULLSPIN_all_singlegene_gradients.mat'));
ngenes = size(data.pval, 1);
abonf = 0.05 ./ 7;


% fdr correction
for ii = 1:size(data.pval, 1)
    data.padj(ii, :) = mafdr(data.pval(ii, :), 'BHFDR', true);
end
% % % save(fullfile(datapath, 'LR_NULLSPIN_all_singlegene_gradients.mat'),...
% % %     '-struct', 'data');


% significant associations in LR
% -- uncorrected
numSig(1,1) = nnz(data.pval(:) < 0.05); 
ratSig(1,1) = numSig(1,1)./ numel(data.pval(:));
disp(['## LR, uncorrected: ', num2str(numSig(1,1)), ', ' , ...
    num2str(ratSig(1,1)*100), '%']);

% -- fdr
numSig(2,1) = nnz(data.padj(:) < 0.05);
ratSig(2,1) = numSig(2,1)./ numel(data.padj(:));
disp(['## LR, FDR corrected: ', num2str(numSig(2,1)), ', ' , ...
    num2str(ratSig(2,1)*100), '%']);

% -- bonf
numSig(3,1) = nnz(data.pval(:) < abonf);
ratSig(3,1) = numSig(3,1) ./ numel(data.pval(:));
disp(['## LR, Bonferroni corrected: ', num2str(numSig(3,1)), ', ' , ...
    num2str(ratSig(3,1)*100), '%']);


% significant associations in LR & null-spin
% -- uncorrected
numSig(1,2) = nnz((data.pval(:) < 0.05) & (data.pval_nullspin(:) < 0.05));
ratSig(1,2) = numSig(1,2)./ numSig(1,1); % divided by N from LR
disp(['## LR & null-spin, uncorrected: ', num2str(numSig(1,2)), ', ' , ...
    num2str(ratSig(1,2)*100), '%']);

% -- fdr
numSig(2,2) = nnz((data.padj(:) < 0.05) & (data.pval_nullspin(:) < 0.05));
ratSig(2,2) = numSig(2,2)./ numSig(2,1);
disp(['## LR & null-spin, FDR corrected: ', num2str(numSig(2,2)), ', ' , ...
    num2str(ratSig(2,2)*100), '%']);

% -- bonf
numSig(3,2) = nnz((data.pval(:) < abonf) & (data.pval_nullspin(:) < 0.05));
ratSig(3,2) = numSig(3,2)./ numSig(3,1);
disp(['## LR & null-spin , bonf corrected: ', num2str(numSig(3,2)), ', ' , ...
    num2str(ratSig(3,2)*100), '%']);


% significant associations in LR & null-spin & null-random
% -- uncorrected
numSig(1,3) = nnz((data.pval(:) < 0.05) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullrandom(:) < 0.05));
ratSig(1,3) = numSig(1,3)./ numSig(1,2); % divided by N from LR + null-spin
disp(['## LR & null-spin & null-random, uncorrected: ', num2str(numSig(1,3)), ...
    ', ' , num2str(ratSig(1,3)*100), '%']);

% -- fdr
numSig(2,3) = nnz((data.padj(:) < 0.05) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullrandom(:) < 0.05));
ratSig(2,3) = numSig(2,3)./ numSig(2,2);
disp(['## LR & null-spin & null-random, fdr corrected: ', num2str(numSig(2,3)), ...
    ', ' , num2str(ratSig(2,3)*100), '%']);

% -- bonf
numSig(3,3) = nnz((data.pval(:) < abonf) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullrandom(:) < 0.05));
ratSig(3,3) = numSig(3,3)./ numSig(3,2);
disp(['## LR & null-spin & null-random, bonf corrected: ', num2str(numSig(3,3)), ...
    ', ' , num2str(ratSig(3,3)*100), '%']);


% significant associations in LR & null-spin & null-brain
% -- uncorrected
numSig(1,4) = nnz((data.pval(:) < 0.05) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullbrain(:) < 0.05));
ratSig(1,4) = numSig(1,4)./ numSig(1,3); % divided by N from LR + null-random
disp(['## LR & null-spin & null-brain, uncorrected: ', num2str(numSig(1,4)), ...
    ', ' , num2str(ratSig(1,4)*100), '%']);

% -- fdr
numSig(2,4) = nnz((data.padj(:) < 0.05) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullbrain(:) < 0.05));
ratSig(2,4) = numSig(2,4)./ numSig(2,3);
disp(['## LR & null-spin & null-brain, fdr corrected: ', num2str(numSig(2,4)), ...
    ', ' , num2str(ratSig(2,4)*100), '%']);

% -- bonf
numSig(3,4) = nnz((data.pval(:) < abonf) & (data.pval_nullspin(:) < 0.05) ...
    & (data.pval_nullbrain(:) < 0.05));
ratSig(3,4) = numSig(3,4)./ numSig(3,3);
disp(['## LR & null-spin & null-brain, bonf corrected: ', num2str(numSig(3,4)), ...
    ', ' , num2str(ratSig(3,4)*100), '%']);


% make a table
tbl = array2table(numSig');
tbl.Properties.VariableNames = ...
    {'uncorrected', 'FDR-corrected', 'Bonferroni-corrected'};
tbl.Properties.RowNames = ...
    {'LR', 'LR & N-spin', 'LR & N-spin & N-rand', 'LR & N-spin & N-brain'};
disp(tbl);
% % % writetable(tbl, [datapath, '/tbl_single_gene_gradients.xlsx'], ...
% % %     'WriteRowNames', true);


% plot venn diagram
if ~exist('venn')
   warning(['Please install "venn" package to plot venn diagram. ',...
       'Darik (2021). venn (https://www.mathworks.com/matlabcentral/fileexchange/22282-venn), ',...
       'MATLAB Central File Exchange.']); 
% % %    addpath('~/Documents/codes/venn/');
else
    tmp1 = data.padj;
    tmp2 = data.pval_nullspin;
    tmp3 = data.pval_nullbrain;
    % overlap
    c1 = nnz(tmp1(:)<0.05);
    c2 = nnz(tmp2(:)<0.05);
    c3 = nnz(tmp3(:)<0.05);
    i12 = nnz((tmp1(:)<0.05) & (tmp2(:)<0.05));
    i13 = nnz((tmp1(:)<0.05) & (tmp3(:)<0.05));
    i23 = nnz((tmp2(:)<0.05) & (tmp3(:)<0.05));
    i123 = nnz((tmp1(:)<0.05) & (tmp2(:)<0.05) & (tmp3(:)<0.05));

    % color map
    cm = [0.8941, 0.1020, 0.1098; 0.2157, 0.4941, 0.7216; ...
        0.3020, 0.6863, 0.2902];

    figure(1); venn([c1,c2,c3], [i12,i13,i23,i123], ...
        'FaceColor', {cm(1,:), cm(2,:), cm(3,:)}, ...
        'EdgeAlpha', 0, 'FaceAlpha', 0.6, 'ErrMinMode', 'None');
    legend({'LR', 'null-spin', 'null-brain'}); legend boxoff
    axis equal
    axis off
    set(gcf,'units','centimeters','position',[10,10,12,15]);
    set(gca,'FontSize', 8, 'FontName', 'Arial');
% % % saveas(gcf, '../figures/venn_all_singlegene_gradients.svg');
end
