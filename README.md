# OTFS Radar Sensing

This is the code for [ISAC MIMO Systems With OTFS Waveforms and Virtual Arrays
](https://ieeexplore.ieee.org/abstract/document/11159304)

**Update (2026-09-14):** Added fractional delay support in channel simulation and SSR estimation, alongside fractional Doppler.

**Update (2026-11-08):** Coarse xcorr is now FFT circular cross-correlation against the transmitted DD grid (pilots included).

## Requirements

- MATLAB R2024b (verified) and above
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for `findpeaks`
- Parallel Computing Toolbox only when `trials > 1` and `params.n_workers > 0`

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

results = run_sim(params);
```

To run the default plotting script:

```matlab
run(fullfile(project_dir, 'run_single.m'));
```

## Configuration

`config_params_HF` supports three target modes, in priority order:

1. `targets`: explicit deterministic target structs.
2. `n_targets`: random selection from the packaged pool; delay/Doppler range overrides sample the `delay_resolution`/`doppler_resolution` grid.
3. `selected_targets`: fixed indices into the packaged pool.

SSR stops at `targets_num` by default; set `params.ssr.max_J` to a positive integer to override this limit (`[]` = default). Residual-based stopping may stop earlier.

Fractional delay: set `targets.delay` directly or use `params.epsilon` to override its fractional part (`[]` preserves input delays; scalar or per-target values in `[0,1)`). `params.angle_resolution`, `delay_resolution`, `doppler_resolution` (default `0.1`) round the generated truth and are the default SSR stopping steps (`params.ssr.min_step.a/r/v` override). SSR stage 1 searches angle at `ssr.length.a` and delay/Doppler at the min steps; stage 2 refines each selected angle at the angle min step. Coarse xcorr stays on integer taps/degrees; the fused estimator keeps its 1/16-tap refinement.

`params.padLen` is the fixed guard length (samples, chosen by design; `generate_channel_params` errors if a delay exceeds it). Noise: `n0 = P_rx/(1+SNR)`; the Wiener pilot vector `y = Y*conj(X)/(|X|^2+n0)` follows `y = c*Phi*gains + noise` with `c` and the noise variance from `pilot_noise_model` (constant-modulus pilots), used by all gain fits.

Trial `i` is seeded with `rng(params.rng_seed_start + i - 1, 'threefry')` (same seeds at every SNR; serial and parallel runs are identical; `run_sim` restores the caller's RNG state). Trials run in parallel on `params.n_workers` workers when `trials > 1` (`0` = serial). BER (LMMSE, `detect_lmmse`) is computed in `run_sim` for `trials == 1`; for saved multi-trial results `run_ber_batch('results/<file>.mat')` regenerates each trial from its seed and equalizes with the stored estimates (no re-estimation), saving `results/ber_<file>.mat`.

Pilots (`config_pilots`): `pilot_blocks_def` gives structured delay/Doppler/random blocks per Tx, and every Tx without a block gets `n_pilotspertx` random pilots per trial (scalar or one value per Tx, `0` = none); an empty `pilot_blocks_def` makes all Tx random. A TF bin may carry only one Tx (`config_pilots` errors otherwise).

`coarse_method`: `'xcorr'` (default) or `'fused'` (`estimate_coarse_params_fused`, needs delay and Doppler blocks in `pilot_blocks_def`): xcorr pairs angle/delay/Doppler; CLEAN: candidates snapped on pilot-arm spectra, refined to 1/16 tap on a local grid by coherent LS over all pilots, kept if they reduce the residual.

## Outputs

`run_sim` returns:

- `err_coarse`: coarse-estimator errors with dimensions `(trials, SNR, target, parameter)`.
- `err_ssr`: refinement errors with the same dimensions.
- `rng_seeds`: `(trials, 1)` seed of each trial.
- `ests_coarse`, `ests_ssr`: `(trials, SNR)` cells of structs with the estimates (`AoAs`, `AoDs`, `pathDelays`, `pathDopplers`, `pathGains`, one entry per estimated path; SSR also `residual_nmse`); `[]` for failed trials or skipped SSR.
- `n_detected_coarse`, `n_detected_ssr`: estimated path count per trial and SNR.
- `n_matched_coarse`, `n_matched_ssr`: estimates matched to a true target (non-NaN error rows).
- `BER`: `(trials, SNR)`, only for `trials == 1` (SSR estimates, coarse if SSR is skipped), else `NaN`.
- `failed`: logical `(trials, SNR)` flag for trials that errored in the parallel path.
- `SNR_db_list`: simulated SNR values.
- `params`: resolved configuration.
- `chanParams`: `(trials, SNR)` cells of structs with the true channel; `[]` for failed trials.

The parameter dimension is `[angle, delay, Doppler]`. Replot saved results with `plot_results(results)`.

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
