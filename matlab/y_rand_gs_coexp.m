function [gene_id, coexp, status] = y_rand_gs_coexp(G, T, NG, maxDiff)
% y_rand_gs_coexp() looks for random gene sets with similar coexpression
% level as the input genes of interest
% Input:
%   G -- gene expression matrix, regions by genes
%   T -- target mean coexpression level
%   NG -- number of genes within the gene sets
% Output:
%   gene_id -- index of random genes in the gene expression matrix G
%   coexp -- mean coexpression level among random genes
%   status -- indicator of success performance. false if the function fails
% by Yongbin Wei, 2020, VU University Amsterdam

if nargin == 3
    maxDiff = 0.025;
end

nGenes = size(G, 2);
nCount = 1;
status = true;

% initialize a set of random genes
gene_id = randperm(nGenes, NG); % random gene set S1
Grand = G(:, gene_id); % expression matrix of S1

% compute mean coexpression for each gene and for the whole gene set
[~, tmp_coexp_mean, tmp_coexp_mean_mean] = compute_mean_coexp(Grand); 
coexp = tmp_coexp_mean_mean;
delta_coexp = tmp_coexp_mean_mean - T; % difference from the targ

% compute coexpression between the gene set to the rest of genes
Gsub_G_coexp = compute_Gsub_G_coexp(Grand, G);
[~, Gsub_G_coexp_sorted_idx] = sort(Gsub_G_coexp, 'ascend'); % + -> 0
Gsub_G_coexp_sorted_idx(ismember(Gsub_G_coexp_sorted_idx, gene_id)) = [];

nRep = ceil(NG/100);

% do if delta_coexp > maxdiff
while abs(delta_coexp) > maxDiff
    nCount = nCount + 1;
    
    % in case we need to decrease coexpression level
    if delta_coexp > 0
        % find gene with max coexpression and drop it
        [~, idx_max] = maxk(tmp_coexp_mean, nRep);
        Grand(:, idx_max) = [];        
        gene_id(idx_max) = [];
        
        % add a new gene
        rid = Gsub_G_coexp_sorted_idx(1:nRep); % replace with the one differs the most
        Gsub_G_coexp_sorted_idx(1:nRep) = [];
        
        Grand = [Grand, G(:, rid)];
        gene_id = [gene_id, rid];
           
        % compute delta coexp after add this gene
        [~, tmp_coexp_mean, tmp_coexp_mean_mean] = compute_mean_coexp(Grand);
        
        delta_coexp = tmp_coexp_mean_mean - T;                
        coexp = tmp_coexp_mean_mean;
    % in case we need to incease coexpression level
    else
        % find gene with min coexpression and drop it
        [~, idx_min] = mink(tmp_coexp_mean, nRep);
        Grand(:, idx_min) = [];        
        gene_id(idx_min) = '';
        
        % add a new gene  
        rid = Gsub_G_coexp_sorted_idx((end-nRep+1):end); % replace with the one differs the most
        Gsub_G_coexp_sorted_idx((end-nRep+1):end) = [];
        Grand = [Grand, G(:, rid)];
        gene_id = [gene_id, rid];
           
        % compute delta coexp after add this gene
        [~, tmp_coexp_mean, tmp_coexp_mean_mean] = compute_mean_coexp(Grand);
        
        delta_coexp = tmp_coexp_mean_mean - T;            
        coexp = tmp_coexp_mean_mean;
    end
        
    if nCount == 500
%         disp('ERROR: fail to find random genes. Please use larger maxdiff and retry');
        status = false;
        break
    end    
end
% if status == true 
%     disp('## Finished without errors');
% end

end

% compute mean coexpression
function [coexp_mat, coexp_mean, coexp_mean_mean] = compute_mean_coexp(G_sub)
    coexp_mat = corr(G_sub, 'rows', 'pairwise');
    coexp_mean = nanmean(coexp_mat, 1);
    coexp_mean_mean = nanmean(coexp_mat(:));
end

% compute coexpression from geneset to all
function Gsub_G_coexp = compute_Gsub_G_coexp(Gsub, G)
    rtmp = corr(Gsub, G, 'rows', 'pairwise');    
    Gsub_G_coexp = nanmean(rtmp, 1);
end
