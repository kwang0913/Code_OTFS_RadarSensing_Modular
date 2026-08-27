# OTFS Radar Sensing

Self-contained MATLAB implementation of monostatic OTFS radar/DFRC parameter estimation. The default workflow mirrors `OTFS_DFRC_Modula`: waveform generation, channel simulation, coarse angle-delay-Doppler estimation, successive-refinement estimation, and Monte Carlo reporting.

> **Maintainer space:** Replace this note with a one-paragraph project description before release.

## Contents

| Path | Purpose |
|---|---|
| `run_single.m` | Default simulation, reporting, and plots |
| `run_sim_new.m` | Current SNR-by-trial workflow used by `run_single.m` |
| `run_sim.m` | Original coarse-estimator workflow retained for comparison |
| `config_params_HF.m` | Default 24.25 GHz monostatic system configuration |
| `estimate_coarse_params_new.m` | Current coarse angle-delay-Doppler estimator |
| `estimate_ssr_params.m` | Successive-refinement estimator |
| `generate_TF.m` | OTFS waveform and private-pilot generation |
| `generate_channel_params.m` | Deterministic, random, or pool-based target generation |
| `targets/` | Packaged target pools |

All runtime source and default target data are included. Generated results, plots, local tooling state, archives, and MATLAB scratch files are excluded by `.gitignore`.

## Requirements

- MATLAB R2024b (verified)
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for `findpeaks`
- Parallel Computing Toolbox only when `trials > 1`


## Quick start

The package resolves its target pool relative to its own location, so the caller does not need to change the working directory.

```matlab
project_dir = '/path/to/OTFS_DFRC_Sensing';
addpath(project_dir);

params = config_params_HF(struct( ...
    'SNR_db_list', 20, ...
    'trials', 1, ...
    'show_figures', false, ...
    'verbose_output', true));

results = run_sim_new(params);
```

To run the default plotting script:

```matlab
run(fullfile(project_dir, 'run_single.m'));
```

## Configuration

`config_params_HF` supports three target modes, in priority order:

1. `targets`: explicit deterministic target structs.
2. `n_targets`: random selection from the packaged pool, optionally modified by angle, delay, or Doppler ranges.
3. `selected_targets`: fixed indices into the packaged pool.

A deterministic compact smoke configuration:

```matlab
target = struct('angle', 0, 'delay', 2, 'doppler', 1);
params = config_params_HF(struct( ...
    'M', 32, ...
    'N', 16, ...
    'Nt', 1, ...
    'Nr', 1, ...
    'n_private', 16, ...
    'targets', target, ...
    'n_targets', [], ...
    'SNR_db_list', 20, ...
    'trials', 1, ...
    'show_figures', false, ...
    'verbose_output', false));
results = run_sim_new(params);
```

`gain` may be a scalar applied to every target or one value per target. Leave it empty for random complex-Rayleigh gains.

## Outputs

`run_sim_new` returns:

- `pct_err_coarse`: coarse-estimator errors with dimensions `(trials, SNR, target, parameter)`.
- `pct_err_ssr`: refinement errors with the same dimensions.
- `n_detected`: target count per trial and SNR.
- `SNR_db_list`: simulated SNR values.
- `params`: resolved configuration.
- `chanParams`: final single-trial ground truth, when available.

The parameter dimension is `[angle, delay, Doppler]`.

## Release metadata

- **License:** BSD 3-Clause; see [`LICENSE`](LICENSE).
- **Citation:** [Add preferred paper citations and BibTeX.]
- **Authors/contact:** [Add public maintainer details.]
- **Repository:** <https://github.com/kwang0913/OTFS_DFRC_Sensing>
- **Third-party material:** Confirm redistribution rights before adding imported code or paper PDFs.
- **Validation:** MATLAB R2024b on Linux; compact single-target smoke workflow above.

## Citation

[Add the preferred citation text and BibTeX here.]

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
