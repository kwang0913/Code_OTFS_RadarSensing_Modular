function [per_target_err, matched_ests] = compute_errors(chanParams, ests, params)
%COMPUTE_ERRORS Per-target absolute errors of estimated vs true parameters.
%   [per_target_err, matched_ests] = compute_errors(chanParams, ests, params)
%
%   per_target_err: Kx5 matrix, columns = [angle (AoA), delay, doppler, gain_magnitude, AoD].
%                   Rows aligned to the K true targets (matched jointly, see below).
%                   Unmatched (missed or mis-assigned) targets have NaN.
%                   The gain column is NaN unless both sides carry one gain per
%                   path (per-Tx gain matrices are not compared); the AoD column
%                   is NaN unless both sides carry AoDs (monostatic: AoD = AoA).
%   matched_ests:   Kx5 matrix of estimated values in the same target order
%                   (gain as magnitude). Unmatched targets have NaN.
%
%   Paths flagged is_los (the known direct path) are excluded on both sides.
%   Gated greedy matching: a pair is eligible only if |dAoA|, |dDelay|,
%   |dDoppler| are within params.ssr.range (a, r, v; angle ignored if
%   params.skip_aoa). Eligible pairs are assigned one-to-one in ascending
%   order of the sum of squared range-normalised errors.

    t = 1:numel(chanParams.AoAs);
    if isfield(chanParams, 'is_los'), t = find(~chanParams.is_los); end
    e = 1:numel(ests.AoAs);
    if isfield(ests, 'is_los'), e = find(~ests.is_los); end
    n_true = numel(t);
    n_est  = numel(e);
    n_match = min(n_est, n_true);
    use_aod  = isfield(chanParams, 'AoDs') && isfield(ests, 'AoDs');
    use_gain = isfield(chanParams, 'pathGains') && isfield(ests, 'pathGains') && ...
        isvector(chanParams.pathGains) && isvector(ests.pathGains);

    per_target_err = NaN(n_true, 5);
    matched_ests   = NaN(n_true, 5);

    if n_est == 0 || n_true == 0, return; end

    rg = params.ssr.range;
    col = @(x, i) reshape(x(i), [], 1);  row = @(x, i) reshape(x(i), 1, []);
    da = abs(col(chanParams.AoAs, t) - row(ests.AoAs, e)) / rg.a;
    dr = abs(col(chanParams.pathDelays, t) - row(ests.pathDelays, e)) / rg.r;
    dv = abs(col(chanParams.pathDopplers, t) - row(ests.pathDopplers, e)) / rg.v;
    if isfield(params, 'skip_aoa') && params.skip_aoa, da(:) = 0; end
    dist = da.^2 + dr.^2 + dv.^2;
    dist(da > 1 | dr > 1 | dv > 1) = Inf;  % outside gate
    [~, order] = sort(dist(:));
    used_true = false(1, n_true);
    used_est  = false(1, n_est);
    for ii = 1:length(order)
        [ti, ei] = ind2sub([n_true, n_est], order(ii));
        if isinf(dist(ti, ei)), break; end  % remaining pairs outside gate
        if used_true(ti) || used_est(ei), continue; end
        used_true(ti) = true;
        used_est(ei)  = true;
        T = t(ti);  E = e(ei);
        matched_ests(ti, 1:3) = [ests.AoAs(E), ests.pathDelays(E), ests.pathDopplers(E)];
        per_target_err(ti, 1:3) = abs(matched_ests(ti, 1:3) - ...
            [chanParams.AoAs(T), chanParams.pathDelays(T), chanParams.pathDopplers(T)]);
        if use_gain
            matched_ests(ti, 4) = abs(ests.pathGains(E));
            per_target_err(ti, 4) = abs(matched_ests(ti, 4) - abs(chanParams.pathGains(T)));
        end
        if use_aod
            matched_ests(ti, 5) = ests.AoDs(E);
            per_target_err(ti, 5) = abs(matched_ests(ti, 5) - chanParams.AoDs(T));
        end
        if sum(used_true) >= n_match, break; end
    end
end
