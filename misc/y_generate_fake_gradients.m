clc, clear, close all

fsPath = fullfile(fileparts(mfilename('fullpath')), '../data/MNI152_FS/');

[vertices, label, ctab] = read_annotation(...
    fullfile(fsPath, 'label', 'lh.lausanne120.annot'));

[coords, faces] = read_surf(fullfile(fsPath, 'surf', 'lh.pial'));

NN = size(coords, 1);
values = zeros(NN ,1);

% gradient x
[~, xsortidx] = sort(coords(:, 1));
values(xsortidx) = [1:NN] ./ NN;
for ii = 1:ctab.numEntries
    II = label == ctab.table(ii, 5);
    regionProperties(ii, 1) = nanmean(values(II));
end

% gradient y
[~, ysortidx] = sort(coords(:, 2));
values(ysortidx) = [1:NN] ./ NN;
for ii = 1:ctab.numEntries
    II = label == ctab.table(ii, 5);
    regionProperties(ii, 2) = nanmean(values(II));
end

% gradient z
[~, zsortidx] = sort(coords(:, 3));
values(zsortidx) = [1:NN] ./ NN;
for ii = 1:ctab.numEntries
    II = label == ctab.table(ii, 5);
    regionProperties(ii, 3) = nanmean(values(II));
end

regionDescriptions = strcat('ctx-lh-', ctab.struct_names);

II = contains(regionDescriptions, {'unknown', 'corpuscallosum'});
regionDescriptions(II) = '';
regionProperties(II, :) = [];

save('../data/gradients_dk114.mat', 'regionProperties', 'regionDescriptions');
