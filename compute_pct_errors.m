function [per_target_err, matched_ests] = compute_pct_errors(chanParams, ests)
%COMPUTE_PCT_ERRORS Per-target absolute errors (monostatic).
%   [per_target_err, matched_ests] = compute_pct_errors(chanParams, ests)
%
%   Computes per-target absolute errors between true and estimated
%   channel parameters. For monostatic DFRC, AoA = AoD, so 4 columns.
%
%   per_target_err: Kx4 matrix, columns = [angle, delay, doppler, gain_magnitude].
%                   Rows aligned to true targets (matched by nearest AoA).
%                   Unmatched targets have NaN.
%   matched_ests:   Kx4 matrix of estimated values in the same target order.
%                   Unmatched targets have NaN.
%
%   Greedy matching: each estimate is assigned to the closest unmatched
%   true target by AoA.

    n_est  = length(ests.AoAs);
    n_true = length(chanParams.AoAs);
    n_match = min(n_est, n_true);

    per_target_err = NaN(n_true, 4);
    matched_ests   = NaN(n_true, 4);

    if n_est == 0 || n_true == 0, return; end

    % Build distance matrix explicitly to avoid implicit expansion issues
    dist = zeros(n_true, n_est);
    for ti = 1:n_true
        for ei = 1:n_est
            dist(ti, ei) = abs(chanParams.AoAs(ti) - ests.AoAs(ei));
        end
    end
    [~, order] = sort(dist(:));
    used_true = false(1, n_true);
    used_est  = false(1, n_est);
    for ii = 1:length(order)
        [ti, ei] = ind2sub([n_true, n_est], order(ii));
        if used_true(ti) || used_est(ei), continue; end
        used_true(ti) = true;
        used_est(ei)  = true;
        per_target_err(ti, :) = [ ...
            abs(chanParams.AoAs(ti) - ests.AoAs(ei)), ...
            abs(chanParams.pathDelays(ti) - ests.pathDelays(ei)), ...
            abs(chanParams.pathDopplers(ti) - ests.pathDopplers(ei)), ...
            abs(abs(chanParams.pathGains(ti)) - abs(ests.pathGains(ei)))];
        matched_ests(ti, :) = [ ...
            ests.AoAs(ei), ests.pathDelays(ei), ...
            ests.pathDopplers(ei), abs(ests.pathGains(ei))];
        if sum(used_true) >= n_match, break; end
    end
end
