# OTFS Radar Sensing

This is the code for [ISAC MIMO Systems With OTFS Waveforms and Virtual Arrays
](https://ieeexplore.ieee.org/abstract/document/11159304)

**Update (2026-09-14):** Added fractional delay support in channel simulation and SSR estimation, alongside fractional Doppler.

## Requirements

- MATLAB R2024b (verified) and above
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for `findpeaks`
- Parallel Computing Toolbox only when `trials > 1`


## Quick start

The package resolves its target pool relative to its own location, so the caller does not need to change the working directory.

```matlab
project_dir = '/path/to/Code_OTFS_RadarSensing_Modular';
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
2. `n_targets`: random selection from the packaged pool; delay/Doppler range overrides sample a `0.1`-tap grid.
3. `selected_targets`: fixed indices into the packaged pool.

SSR stops at `targets_num` by default; set `params.ssr.max_J` to a positive integer to override this limit (`[]` = default). Residual-based stopping may stop earlier.

Fractional delay: set `targets.delay` directly or use `params.epsilon` to override its fractional part (`[]` preserves input delays; scalar or per-target values in `[0,1)`). SSR delay/Doppler grid steps are `params.ssr.length.r/v` (both default to `0.1` tap).

## Outputs

`run_sim_new` returns:

- `pct_err_coarse`: coarse-estimator errors with dimensions `(trials, SNR, target, parameter)`.
- `pct_err_ssr`: refinement errors with the same dimensions.
- `n_detected`: target count per trial and SNR.
- `SNR_db_list`: simulated SNR values.
- `params`: resolved configuration.
- `chanParams`: final single-trial ground truth, when available.

The parameter dimension is `[angle, delay, Doppler]`.

## Citation

```bibtex
@article{wang2026isacmimo,
  author  = {Wang, Kailong and Petropulu, Athina},
  title   = {ISAC MIMO Systems With OTFS Waveforms and Virtual Arrays},
  journal = {IEEE Journal on Selected Areas in Communications},
  year    = {2026},
  volume  = {44},
  pages   = {229--244},
  doi     = {10.1109/JSAC.2025.3608761}
}
```

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
