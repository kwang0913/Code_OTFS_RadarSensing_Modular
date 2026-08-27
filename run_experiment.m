function run_experiment(json_file)
%RUN_EXPERIMENT JSON-driven batch runner with sweep support.
%
%   run_experiment('experiments/my_config.json')
%
%   JSON format:
%     {
%       "overrides": { "pilot_power": 10, "trials": 50, "seed": 2 },
%       "sweep":     [ {"pilot_power": 1}, {"pilot_power": 5} ]
%     }
%
%   Results saved to results/<json_name>.mat (or <json_name>_run<N>.mat for sweeps)
%
%   Bash: matlab -batch "run_experiment('experiments/my_config.json')"

    cfg = jsondecode(fileread(json_file));

    if isfield(cfg, 'overrides')
        base_overrides = cfg.overrides;
    else
        base_overrides = struct();
    end

    % sweep: array of override structs; each element triggers a separate run
    if isfield(cfg, 'sweep')
        sweep_list = cfg.sweep;
        n_runs = numel(sweep_list);
    else
        sweep_list = {};
        n_runs = 1;
    end

    [~, json_name] = fileparts(json_file);
    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    fprintf('=== Experiment: %s (%d run(s)) ===\n', json_name, n_runs);

    for run_idx = 1:n_runs
        if ~isempty(sweep_list)
            if iscell(sweep_list)
                sweep_overrides = sweep_list{run_idx};
            else
                sweep_overrides = sweep_list(run_idx);
            end
            run_overrides = apply_overrides(base_overrides, sweep_overrides);
        else
            run_overrides = base_overrides;
        end

        params = config_params_HF(run_overrides);

        fprintf('\n--- Run %d/%d ---\n', run_idx, n_runs);
        override_fields = fieldnames(run_overrides);
        for f = 1:numel(override_fields)
            key = override_fields{f};
            val = run_overrides.(key);
            if isnumeric(val) && isscalar(val)
                fprintf('  %s = %g\n', key, val);
            elseif isnumeric(val)
                fprintf('  %s = [%s]\n', key, num2str(val));
            elseif ischar(val) || isstring(val)
                fprintf('  %s = %s\n', key, val);
            end
        end

        results = run_sim(params);

        if n_runs > 1
            fname = sprintf('%s_run%d.mat', json_name, run_idx);
        else
            fname = build_result_filename(params, json_name);
        end
        out_file = fullfile(out_dir, fname);
        save(out_file, 'results');
        fprintf('  Saved: %s\n', out_file);
    end

    fprintf('\n=== Done: %d run(s) completed ===\n', n_runs);
end
