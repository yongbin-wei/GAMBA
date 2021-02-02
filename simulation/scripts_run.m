%% all pheno
expressionData = '../data/gene_expression.mat';
imgData = '../data/img_data.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_singlegene_all_pheno(expressionData, imgData, spinDir, outputDir);


%% gradients
expressionData = '../data/gene_expression.mat';
gradientsData = '../data/gradients_dk114.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_singlegene_gradients(expressionData, gradientsData, spinDir, outputDir);


%% all pheno (GO)
expressionData = '../data/gene_expression.mat';
imgData = '../data/img_data.mat';
goData = '../data/GOterms_BP.mat';
spinDir = '../data/gene_expression_spinned';
outputDir = './output/';

scripts_all_GO_all_pheno(expressionData, imgData, goData, spinDir, outputDir);


