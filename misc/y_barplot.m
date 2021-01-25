function y_barplot(vals, sig_index, x_labels, y_labels, output_file, isylabel, W, H)
% This function produces plot that resembles GAMBA bar charts

if ~exist('W', 'Var')
    W = 6;
end

if ~exist('H', 'Var')
    H = 4;
end

vals = flipud(vals);

II_sig = zeros(numel(vals), 1);
II_sig(sig_index) = 1;
II_sig = flipud(II_sig);

figure
barh(vals, 'EdgeAlpha', 0, 'FaceColor', [209,226,243]./255); hold on;

tmpvals = vals;
tmpvals(~II_sig) = 0;
barh(tmpvals, 'EdgeAlpha', 0, 'FaceColor', [20,101,173]./255);

% set parameters
if ~strcmp(x_labels,'z score')
    xlim([-0.75,1]);
    xticks([-0.75:0.25:0.75]);
else
    xlim([-6,6]);
    xticks([-6:2:6]);
end

yticklabels(flipud(y_labels));
set(gca,'FontSize', 7, 'FontName', 'Arial');
xlabel(x_labels);
ax = gca;
ax.XGrid = 'on';
box off;

if isylabel==0
    ax.YAxis.Visible = 'off';
end

set(gcf,'units','centimeters','position',[10,10,W,H]);
print(output_file,'-dsvg')

end