% ======= Performing statistical analysis for GAMBA =======

% Step 1. Download necessary input data
scripts_download.m

% Step 2. Run null-random, null-brain, null-coexpression models
scripts_statistics.m

% Step 3. Run null-spin model
scripts_statistics_spin.m

% Step 4. Make json files for GAMBA
scripts_make_json.m



