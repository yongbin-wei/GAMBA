% ======= Performing statistical analysis for GAMBA =======

% Step 1. Run null-random, null-brain, null-coexpression, null-spin models
scripts_statistics.m

% Step 2. Make json files for GAMBA
scripts_make_json.m

% NOTE: there are several output files in the 'output' folder. These files are used in '../examples'.

