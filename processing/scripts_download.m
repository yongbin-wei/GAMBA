% Download necessary data for statistical analysis in GAMBA
clc, clear, close all

if ~exist('input', 'dir')
  mkdir input
end

cmd = 'wget "https://dl.dropbox.com/s/dzomugk2m5m5gc6/input.mat?dl=1" -O ./input/input.mat --no-check-certificate';
system(cmd)


