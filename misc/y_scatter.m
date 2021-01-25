function y_scatter(x, y, x_label, y_label, dir_out, W, H)
% This function produces scatter plots

if ~exist('W', 'Var')
    W = 3.5;
end

if ~exist('H', 'Var')
    H = 2.8;
end

figure;
scatter(x, y, 12, [120,120,120]./255, 'filled'); 

lsline;
xlabel(x_label)
ylabel(y_label)
set(gca,'FontSize', 8, 'FontName', 'Calibri');
set(gcf,'units','centimeters','position',[10,10,W,H]);
print(dir_out,'-dsvg');

end