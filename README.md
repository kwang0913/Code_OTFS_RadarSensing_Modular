# OTFS Radar Sensing

This is the code for [ISAC MIMO Systems With OTFS Waveforms and Virtual Arrays
](https://ieeexplore.ieee.org/abstract/document/11159304)

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
2. `n_targets`: random selection from the packaged pool, optionally modified by angle, delay, or Doppler ranges.
3. `selected_targets`: fixed indices into the packaged pool.

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
@article{11159304,
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
