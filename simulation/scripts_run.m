%% all pheno
expressionData = '../data/GE_lausanne120_2mm_GAMBA_20200803.mat';
imgData = '../data/IMG_DATA_ALL_SINGLE_GENE_20200803.mat';
braingenesData = '../data/GTEx_brain_genes_0.05_updated.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_singlegene_all_pheno(expressionData, imgData, ...
    braingenesData, spinDir, outputDir)


%% gradients
expressionData = '../data/GE_lausanne120_2mm_GAMBA_20200803.mat';
gradientsData = '../data/gradients_dk114.mat';
braingenesData = '../data/GTEx_brain_genes_0.05_updated.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_singlegene_gradients(expressionData, gradientsData, ...
    braingenesData, spinDir, outputDir)


%% all pheno (GO)
expressionData = '../data/GE_lausanne120_2mm_GAMBA_20200803.mat';
imgData = '../data/IMG_DATA_ALL_SINGLE_GENE_20200803.mat';
goData = '../data/GOterms_BP.mat';
braingenesData = '../data/GTEx_brain_genes_0.05_updated.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_GO_all_pheno(expressionData, imgData, goData, ...
    braingenesData, spinDir, outputDir)
