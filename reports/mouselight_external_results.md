# MouseLight external test: frozen soma-only geometry on new animals

## Summary

* **What was tested:** frozen 3D axon-geometry methods, predicting complete axons from the soma position only,
  on 97 hypothalamic neurons from 22 MouseLight animals excluded from this pipeline's recorded training and
  development inputs.
* **Result:** neither fixed learned candidate (E4 or E6) beat every fixed control on both mean distance and mean
  F1. **There is no learned win, and nothing is promoted.**
* **Best method:** the simple local control E1 (the animal-balanced medoid of the 16 nearest training somas) had
  both the lowest distance and the highest F1.
* **Verification:** an independent replay verified the source and forecast pins, all 97 lossless trace
  conversions, admission and equal-animal score aggregation. Differences in recomputed means were roundoff only.

## Results

All 97 neurons were admitted and none were excluded. Each neuron's score is averaged within its animal, then
the 22 animals are weighted equally. Lower distance and higher F1 are better.

| method | kind | mean distance (mm) | F1 within 200 um |
|---|---|---:|---:|
| nearest copy | control | 1.051023 | 0.157367 |
| global medoid copy | control | 0.832989 | 0.168743 |
| **E1 local uniform medoid** | control | **0.784120** | **0.208989** |
| E2 local kernel medoid | control | 0.798112 | 0.202019 |
| E4 local kernel regional | learned candidate | 0.815964 | 0.199773 |
| E6 local consensus | learned candidate | 0.871021 | 0.205837 |

**How each candidate did:**

* **E4:** better than both plain copies on distance and F1, but worse than E1 and E2 on both. Gate: not met.
* **E6:** better than nearest copy on both, worse than E1 on both. Against E2 it is a trade-off: higher F1 but
  worse distance. Gate: not met.

The predeclared gate required strictly lower distance AND strictly higher F1 than every control. The
comparison is descriptive, with no significance test.

## Order of events (each step before the next)

| step | time (UTC) or reference |
|---|---|
| 1. Source code and protocol frozen | commit `169b06b7c2902cc912f163f5b0f55111399f9943` |
| 2. Predictions made (582 forecasts = 6 methods x 97 neurons) | 2026-11-12 20:00:03.81 |
| 3. Forecast hashes locked publicly | commit `88493d55ef690d0c624fcf6c795a5326cc99705a`, 20:01:31 |
| 4. Trace download started (official route; 178,607,296 bytes of raw responses) | 20:02:18.79 |
| 5. Single scoring run | 20:03:52.63 to 20:04:02.15 |

* No MouseLight evaluation-query axon existed locally before the lock. Training axons and the separate
  non-test calibration neuron were already available.
* Every trace request was bound to the publicly locked forecasts.
* The score ran once.

## What this does and does not show

* **Inputs:**
  * The queries used only the measured soma and its hemisphere.
  * Missing molecular labels stayed empty.
  * No query axon, target profile or image feature chose a donor.
* **Not a gene-model test:** MouseLight has no measured gene profiles. This tests a 3D geometry module only and
  does not validate any gene-based model or gene-to-circuit inference.
* **No circuit claim:** F1 is the harmonic mean of predicted-cable precision and measured-cable recall within
  200 um. It is not the share of circuits reconstructed. Nothing here shows complete or accurate circuit
  reconstruction, new connections or synapses.
* **Shared atlas:** MouseLight and the training data share the same atlas frame. One separate cortical neuron
  from a non-test animal checked the format conversion and the physical basis. That check did not establish
  how accurate each laboratory's registration is.
* **Exposure:** the 97 neurons and 22 animals are now exposed (`docs/mouselight_external_exposure.json`). They
  can be used for development, but no longer as an untouched validation set.

## Related raw-image work

A separate pilot (`docs/raw_brain_image_pilot_protocol.json`, `analysis/raw_brain_image_features.py`, commit
`5937437c1e0a904dded7c0de2d268ae189c0e36c`) read real MouseLight fluorescence imagery.

**What it did:**

* read a whole-brain root-level TIFF and a high-resolution PVH-region texture with a strict decoder;
* measured intensity, segmentation uncertainty, scale and orientation in physical units.

**Limits:**

* At high resolution the dominant direction followed the axial (z) axis, which is consistent with optical blur
  rather than fibre direction.
* The pilot supplies measurements only, not traced cells, genes or validated connections.
* It was not an input to this test.

## Files

| file | SHA-256 |
|---|---|
| `reports/mouselight_external_metrics.json` (exact run metrics) | `6b760fe7de833adb14005ae056250ab60e41104972336e280a5e56622c7e5482` |
| `reports/mouselight_external_replay.json` (independent replay, PASS) | `2dff8ed6c81869419138a6b90ce0070fd6ccf1e4a93ed28a504e929892b55a94` |
| `reports/mouselight_external_results.json` | `c2cff732426a0525ccbb01a74c2a894aa950f7203db6e6fac79a600be6f77a1d` |
| `docs/mouselight_external_exposure.json` | `bf09ae93203e98687aa52effd7462e8df6da864290d0dfa708fce3bff3d5e943` |

Pinned run records: freeze `b3d40c26…`, prediction manifest `ecad76d9…`, public prediction commitment
`72956c2f…`.

The traces stay local, with the licence and provenance of each neuron. The figshare terms of each neuron apply.

## Next step

Use the strongest local-medoid control (E1) as the bar for further geometry development. Keep every failure on
record. Do not tune on this cohort or retest it as if it were fresh.
