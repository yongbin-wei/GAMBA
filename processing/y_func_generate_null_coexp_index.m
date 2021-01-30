function y_func_generate_null_coexp_index(inputFile, outputDir, gs_size)

load(inputFile, 'mDataGEctx');

nPerm = 500;

targ_coexp = 0.05: 0.05: 0.3;
nCoexp = numel(targ_coexp);
gs_size = str2num(gs_size);

% initialize
disp(num2str(gs_size));
coexp = [];
gene_idx = [];

% for each coexp level
for jj = 1:nCoexp
    for kk = 1:nPerm %for each perm
        tmp_status = false;
        while tmp_status ~= true
            [rid, coexp(jj, kk), tmp_status] = y_rand_gs_coexp(...
                mDataGEctx, targ_coexp(jj), gs_size);
        end     
        gene_idx(jj, kk, :) = rid;
    end
end

save(fullfile(outputDir, ['idx_rand_',num2str(gs_size),'.mat']), ...
    'gene_idx', 'targ_coexp', 'nPerm', 'coexp');
