% run_single.m
% Main entry script. Runs one simulation, displays results and plots.

close all;
clear all;

% --- Simulation Control ---
params = config_params_HF();

% --- Run Simulation ---
results = run_sim(params);

% --- Save Results (only for multi-trial runs) ---
if params.trials > 1
    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, build_result_filename(params));
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end

% --- Plots and detection table ---
plot_results(results);
