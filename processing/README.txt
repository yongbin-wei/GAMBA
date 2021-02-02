% ======= Performing statistical analysis for GAMBA =======

% Step 1. Run null-random, null-brain, null-coexpression models
scripts_statistics.m

% Step 2. Run null-spin model
scripts_statistics_spin.m

% Step 3. Make json files for GAMBA
scripts_make_json.m

% NOTE: there are several output files in the 'output' folder. These files are used in '../examples'.

