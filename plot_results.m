function plot_results(results)
%PLOT_RESULTS MSE and detection-probability plots/tables from run_sim results.
%   plot_results(results) accepts the struct returned by run_sim or loaded
%   from a saved .mat file. Failed trials (results.failed) are excluded
%   from P_d; results without a `failed` field are assumed to have none.

    params = results.params;
    K = size(results.pct_err_coarse, 3);
    snr = results.SNR_db_list;

    if isfield(results, 'failed')
        failed = results.failed;
    else
        failed = false(size(results.n_detected));
    end

    all_param_labels = {'Angle [deg]', 'Delay [taps]', 'Doppler'};
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
        semilogy(snr, mse_coarse, '-o', 'LineWidth', 2); hold on;
        semilogy(snr, mse_ssr, '-s', 'LineWidth', 2);
        grid on;
        xlabel('SNR [dB]');
        ylabel(['MSE | detected — ' all_param_labels{p}]);
        legend('Coarse', 'SSR');
        set(gca, 'FontSize', 11, 'FontWeight', 'bold');
    end
    sgtitle('Conditional MSE (detected targets only)');

    % --- Detection probability ---
    if params.skip_aoa
        % Nr=1: use delay error for detection (angle is unresolvable)
        dim = 2;  tol = params.ssr.range.r;
        pd_label = sprintf('Detection Probability (delay tol = %d taps)', tol);
        pd_table_label = sprintf('--- Detection Probability (delay error < %d taps) ---', tol);
    else
        % MIMO: use angle error for detection
        dim = 1;  tol = params.ssr.range.a;
        pd_label = sprintf('Detection Probability (angle tol = %d°)', tol);
        pd_table_label = sprintf('--- Detection Probability (angle error < %d°) ---', tol);
    end
    coarse_err = results.pct_err_coarse(:, :, :, dim);
    ssr_err    = results.pct_err_ssr(:, :, :, dim);
    coarse_detected = double(coarse_err < tol);  % NaN error -> not detected
    ssr_detected    = double(ssr_err < tol);
    coarse_detected(repmat(failed, 1, 1, K)) = NaN;  % failed trials are missing
    ssr_detected(repmat(failed, 1, 1, K))    = NaN;

    pd_coarse = mean(mean(coarse_detected, 3), 1, 'omitnan');
    pd_ssr    = mean(mean(ssr_detected, 3), 1, 'omitnan');

    figure;
    plot(snr, pd_coarse, '-o', 'LineWidth', 2); hold on;
    plot(snr, pd_ssr, '-s', 'LineWidth', 2);
    grid on; ylim([0 1.05]);
    xlabel('SNR [dB]'); ylabel('P_d');
    legend('Coarse', 'SSR');
    title(pd_label);
    set(gca, 'FontSize', 11, 'FontWeight', 'bold');

    % --- Per-target table ---
    fprintf('\n%s\n', pd_table_label);
    fprintf('%-16s', 'SNR(dB)');
    fprintf('%8d ', snr);
    fprintf('\n');
    print_pd_rows('Coarse', coarse_detected);
    print_pd_rows('SSR', ssr_detected);
end

function print_pd_rows(name, detected)
    fprintf('  %s:\n', name);
    for k = 1:size(detected, 3)
        fprintf('    Target %-4d', k);
        fprintf('%8.2f ', mean(detected(:, :, k), 1, 'omitnan'));
        fprintf('\n');
    end
    fprintf('    %-10s', 'Overall');
    fprintf('%8.2f ', mean(mean(detected, 3), 1, 'omitnan'));
    fprintf('\n');
end
