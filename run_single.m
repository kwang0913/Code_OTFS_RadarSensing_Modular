% run_single.m
% Main entry script. Runs one simulation, displays results and plots.

close all;
clear all;

% --- Simulation Control ---
params = config_params_HF();

% --- Run Simulation ---
results = run_sim_new(params);

% --- Save Results (only for multi-trial runs) ---
if params.trials > 1
    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, build_result_filename(params));
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end

% --- Plot MSE vs SNR (averaged across targets) ---
all_param_labels = {'Angle [deg]', 'Delay [taps]', 'Doppler'};
K = size(results.pct_err_coarse, 3);

if params.skip_aoa
    mse_dims = [2 3];  % delay, Doppler only
else
    mse_dims = [1 2 3];
end

% --- MSE (conditional on detection) ---
figure;
for si = 1:length(mse_dims)
    p = mse_dims(si);
    subplot(1, length(mse_dims), si);
    mse_coarse = squeeze(mean(mean(results.pct_err_coarse(:, :, :, p).^2, 1, 'omitnan'), 3, 'omitnan'));
    mse_ssr    = squeeze(mean(mean(results.pct_err_ssr(:, :, :, p).^2, 1, 'omitnan'), 3, 'omitnan'));
    semilogy(results.SNR_db_list, mse_coarse, '-o', 'LineWidth', 2); hold on;
    semilogy(results.SNR_db_list, mse_ssr, '-s', 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]');
    ylabel(['MSE | detected — ' all_param_labels{p}]);
    legend('Coarse', 'SSR');
    set(gca, 'FontSize', 11, 'FontWeight', 'bold');
end
sgtitle('Conditional MSE (detected targets only)');

%% --- Detection Probability ---
if params.skip_aoa
    % Nr=1: use delay error for detection (angle is unresolvable)
    delay_tol = params.ssr.range.r;
    coarse_err = results.pct_err_coarse(:, :, :, 2);  % delay error
    ssr_err    = results.pct_err_ssr(:, :, :, 2);
    coarse_detected = ~isnan(coarse_err) & (coarse_err < delay_tol);
    ssr_detected    = ~isnan(ssr_err)    & (ssr_err < delay_tol);
    pd_label = sprintf('Detection Probability (delay tol = %d taps)', delay_tol);
    pd_table_label = sprintf('--- Detection Probability (delay error < %d taps) ---', delay_tol);
else
    % MIMO: use angle error for detection
    ang_tol = params.ssr.range.a;
    coarse_err = results.pct_err_coarse(:, :, :, 1);  % angle error
    ssr_err    = results.pct_err_ssr(:, :, :, 1);
    coarse_detected = ~isnan(coarse_err) & (coarse_err < ang_tol);
    ssr_detected    = ~isnan(ssr_err)    & (ssr_err < ang_tol);
    pd_label = sprintf('Detection Probability (angle tol = %d°)', ang_tol);
    pd_table_label = sprintf('--- Detection Probability (angle error < %d°) ---', ang_tol);
end

% --- P_d vs SNR plot ---
pd_coarse = squeeze(mean(mean(coarse_detected, 1), 3));
pd_ssr    = squeeze(mean(mean(ssr_detected, 1), 3));

figure;
plot(results.SNR_db_list, pd_coarse, '-o', 'LineWidth', 2); hold on;
plot(results.SNR_db_list, pd_ssr, '-s', 'LineWidth', 2);
grid on; ylim([0 1.05]);
xlabel('SNR [dB]'); ylabel('P_d');
legend('Coarse', 'SSR');
title(pd_label);
set(gca, 'FontSize', 11, 'FontWeight', 'bold');

% --- Per-target table ---
fprintf('\n%s\n', pd_table_label);
fprintf('%-16s', 'SNR(dB)');
fprintf('%8d ', results.SNR_db_list);
fprintf('\n');
fprintf('  Coarse:\n');
for k = 1:K
    fprintf('    Target %-4d', k);
    fprintf('%8.2f ', mean(coarse_detected(:, :, k), 1));
    fprintf('\n');
end
fprintf('    %-10s', 'Overall');
fprintf('%8.2f ', mean(mean(coarse_detected, 3), 1));
fprintf('\n');
fprintf('  SSR:\n');
for k = 1:K
    fprintf('    Target %-4d', k);
    fprintf('%8.2f ', mean(ssr_detected(:, :, k), 1));
    fprintf('\n');
end
fprintf('    %-10s', 'Overall');
fprintf('%8.2f ', mean(mean(ssr_detected, 3), 1));
fprintf('\n');
