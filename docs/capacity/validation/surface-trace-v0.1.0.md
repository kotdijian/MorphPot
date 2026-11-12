# SELF TEST — PotterySurfaceTraceAnalyzer v0.1.0

This file records tests that were actually executed against the distributed code.

## Scope

The program has two separate layers:

- `pottery_multiscale_features.py` — common physical-mm multiscale geometry feature layer
- `pottery_surface_trace.py` — CLI and preliminary morphometric state classification

The classifier does **not** assign archaeological identities such as sherd join,
hake-me, crack, repair, or use-wear. It classifies scale/sign/persistence of geometry.

## 1. Syntax / CLI

Actually executed:

```text
python3 -m py_compile pottery_multiscale_features.py pottery_surface_trace.py
PASS

python3 pottery_surface_trace.py --version
pottery_surface_trace.py 0.1.0
PASS

python3 pottery_surface_trace.py --diagnose-env
PASS

python3 pottery_surface_trace.py --help
PASS
```

Test container:

```text
Python  3.13.5
NumPy   2.3.5
SciPy   1.17.0
Trimesh 4.11.1
```

## 2. Real sample — 0015Jinmen_small.ply (700,000 faces)

Actually executed:

```bash
python3 pottery_surface_trace.py 0015Jinmen_small.ply --unit m
```

Observed source geometry:

```text
vertices             : 350,070
faces                : 700,000
boundary edges       : 0
non-manifold edges   : 1
geometry components  : 40
active faces         : 699,614
```

Surface classes:

```text
OUTER      : 235,278
INNER      : 404,730
TRANSITION : 59,606
IGNORED    : 386
```

Physical scale calibration:

```text
median same-surface face step : 0.565160 mm
0.5 mm -> 1 iterations -> 0.565160 mm
1.0 mm -> 3 iterations -> 0.978887 mm
2.0 mm -> 13 iterations -> 2.037715 mm
4.0 mm -> 50 iterations -> 3.996288 mm
```

No requested scale was clipped by the 160-iteration safety cap.

Actual timing recorded by the program:

```text
multiscale feature extraction : 2.063 s
total internal time           : 4.805 s
```

External `/usr/bin/time` wall time in the same run was 5.92 s.

The 0.5 mm request is below one median face step on this 700k-face mesh
(scale/step = 0.885).
It is therefore exported as a diagnostic fine-scale response, not as evidence
that 0.5 mm physical widths can be measured accurately.

PASS.

## 3. Synthetic geometry groove shell

The existing geometry-only synthetic hollow shell from the FragmentBoundaryExtractor
self-test was run through the new common feature layer. It contains no UV evidence.

Observed scale calibration:

```text
median face step : 0.371632 mm
0.5 mm realized  : 0.525566 mm
1.0 mm realized  : 0.983245 mm
2.0 mm realized  : 2.001298 mm
4.0 mm realized  : 4.002595 mm
```

Thus the requested physical scales remain approximately stable despite a different
mesh face spacing from the real sample.

For a conservative geometry mask around the known synthetic groove, measured from
the generated feature PLY:

```text
groove-mask faces                  : 376
dominant response median           : 0.907834
background median                  : 0.003267
background 99 percentile           : 0.179914
persistent trace score median      : 0.818794
background persistent median       : 0.000539
background persistent 99 percentile: 0.072306
```

All 376 groove-mask faces were classified as either:

```text
step_edge_like                    : 188
multiscale_persistent_anomaly     : 188
```

This validates that a known geometric groove produces a concentrated multiscale
response without UV information.

PASS.

## 4. Output read-back / structure

The real-sample output was generated with:

```text
features/outer_multiscale_features_face_centers.ply : 235,278 points
features/inner_multiscale_features_face_centers.ply : 404,730 points
classification/outer_trace_classes_colored.ply      : 235,278 faces
classification/inner_trace_classes_colored.ply      : 404,730 faces
```

The scalar PLY header was inspected and contains, among others:

```text
trace_class
dominant_scale_mm
dominant_response
persistence_fraction
fine_trace_score
medium_trace_score
broad_trace_score
persistent_trace_score
step_edge_score
normal_dev_0p5mm / 1p0mm / 2p0mm / 4p0mm
residual_signed_0p5mm / 1p0mm / 2p0mm / 4p0mm
residual_abs_0p5mm / 1p0mm / 2p0mm / 4p0mm
response_0p5mm / 1p0mm / 2p0mm / 4p0mm
```

PASS.

## 5. Interpretation boundary

The current rule-based states are morphometric descriptions only:

- fine / medium / broad ridge-like
- fine / medium / broad valley-like
- step-edge-like
- multiscale persistent anomaly
- smooth/low response

They are deliberately not archaeological categories. Parallel comparison with
FragmentBoundaryExtractor should determine which combinations of these descriptors
correspond to known sherd joins, hake-me traces, cracks, or other surface phenomena.

## SHA-256

```text
pottery_multiscale_features.py  2f3333aa484c670b07ebba4ba79c75980384caa17e514eec3ca4a8b8296dde1e
pottery_surface_trace.py        3948dabf5a14a3db33eeb1537aef1653b511a63cef499b629ef3a90545a136c3
```
