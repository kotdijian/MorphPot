# MorphPot: Pottery Capacity Estimation

MorphPot is the maintained home of the former PotteryVolumeCalculator tools. Run commands from the MorphPot checkout. The standalone OBJ/PLY command remains available; for normalized PLY, use `python MorphPot.py volume MODEL.ply --unit auto` to read the coordinate unit from metadata without reapplying the pose.

A set of experimental Python tools for estimating pottery vessel capacity from either 3D meshes or 2D archaeological drawings.

The repository currently contains three complementary workflows:

1. **PotteryVolumeCalculator (voxel-fluid)** — estimates retained liquid directly from voxelized 3D free-space connectivity.
2. **PotteryRadialSections** — derives radial section profiles from a 3D mesh and estimates capacity using section-based numerical models.
3. **PotteryDrawingCapacity** — digitizes left/right inner profiles from a 2D drawing and evaluates them with the same `single` integration core used by PotteryRadialSections.

The voxel-fluid workflow runs independently of the section and drawing tools. See the [MorphPot README](../../README.md) for the scripts and documentation included in this checkout.

## Contents

- [PotteryVolumeCalculator — voxel-fluid method](#potteryvolumecalculator-v140)
- [Japanese README](README.md)
- [Repository guide and extension tools](../../README.md)
- [AOM posture comparison](#compare-postures-exported-by-aom)
- [Validation report v1.3.1 (PDF)](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf)
- [Radial section example output](https://github.com/kotdijian/PotteryVolumeCalculator/tree/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/examples/0015Jinmen_small_RadialSections_30deg)

---


## PotteryVolumeCalculator v1.4.0

Experimental Python CLI for estimating the **maximum retained liquid capacity** of a pottery 3D mesh using a voxel-based fluid-space method.

The target quantity is the volume retained immediately before liquid first reaches the exterior, assuming the vessel is upright with **+Z upward** and gravity in the -Z direction.

### Overview

PotteryVolumeCalculator operates directly on the 3D free space defined by a mesh. It does not require an explicit inner-surface extraction or an artificial horizontal cap at the rim.

Core workflow:

```text
PLY / OBJ
  -> geometry-preserving preprocessing and topology QA
  -> exact-coordinate vertex weld
  -> conversion to internal millimetres
  -> whole-mesh surface voxelization
  -> cavity-seed validation
  -> height-limited 3-D flood fill
  -> binary search for safe / first-spill levels
  -> fluid voxel count
  -> retained volume
```

The voxel-fluid implementation is independent of the section- and drawing-based capacity tools.


### Compare postures exported by AOM

In [ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker#pottery-pose), export the same input mesh separately for each posture you want to compare (Rim, Base, or Slice). Pass two or more export folders containing `item_manifest.json`, `transform.json`, and the normalized `*_rev.ply`.

```bash
python3 MorphPot.py volume --compare-poses "AOM/output/cup" "AOM/output/cup/revisions/BASE_REV" "AOM/output/cup/revisions/SLICE_REV" --pitch 1.0 --pose-baseline rim --output-dir "cup_PoseComparison"
```

Use the same original mesh and any identical reconstruction for every posture. The comparison reads the unit from AOM and verifies the two exported files against their manifest digests. It uses the normalized PLY directly. The existing `--validate` workflow checks voxel-pitch convergence separately.

An accepted Slice is included in the primary comparison. Add `--include-reference-slice` to calculate a separately labeled reference Slice. Held or unassessed Slice exports are skipped. The default baseline is the first successful primary mode in Rim, Base, accepted Slice order. Use a fresh output directory. The JSON and CSV summaries record the baseline, angles, spill levels, capacity differences, and any skips or failures. Each posture has its own output subfolder.


### Capacity definition

The reported capacity is:

> The maximum liquid volume that can be retained by the existing geometry immediately before the connected fluid region first reaches the exterior.

For a complete vessel this approximates the vessel's maximum retained capacity. For a partially preserved vessel with a continuous bottom and walls, the result is the capacity retained by the **surviving geometry**, not a reconstruction of the original complete vessel.

### Intended input

- PLY or OBJ triangular mesh
- known coordinate unit: `mm`, `cm`, or `m`
- vessel oriented with +Z upward
- continuous geometry sufficient to define a fluid-retaining cavity

Large holes, missing walls, or openings are treated as potential leakage paths. The program does not automatically reconstruct missing vessel geometry.

### Installation

Python 3 is required. Install either the core dependencies:

```bash
python3 -m pip install -r requirements.txt
```

or the full dependency set including the optional PyMeshLab cross-check:

```bash
python3 -m pip install -r requirements-qa.txt
```

Primary dependencies are NumPy, SciPy, Trimesh, and Pillow. PyMeshLab is optional and is not required for the capacity calculation.

### Usage

#### Single-pitch calculation

```bash
python3 vessel_voxel_volume.py model.ply --unit m --pitch 1.0
```

`--pitch` is always specified in millimetres, regardless of the input coordinate unit.

#### Multi-pitch validation

```bash
python3 vessel_voxel_volume.py model.ply --unit m --validate
```

Validation mode runs three independent full pipelines:

```text
2.0 mm
1.0 mm
0.5 mm
```

Use `--output-dir` for reproducible runs that should not share files with earlier executions.

```bash
python3 vessel_voxel_volume.py model.ply --unit m --pitch 1.0 --output-dir run_01
```

Environment diagnostics:

```bash
python3 vessel_voxel_volume.py --diagnose-env
```

### Method

#### 1. Geometry-preserving preprocessing

The processing path removes unreferenced vertices and merges vertices only when their XYZ coordinates are **exactly identical**. No tolerance-based welding, automatic hole filling, non-manifold repair, or vertex displacement is applied.

#### 2. Surface voxelization

The complete mesh surface is voxelized at pitch `p`. Surface voxels are treated as barriers; candidate liquid occupies connected free voxels.

Because the barrier has finite voxel thickness, the result is resolution-dependent. Multi-pitch validation is therefore part of the intended research workflow.

#### 3. Cavity-seed validation

A candidate seed must represent the vessel cavity rather than a free-space region between wall surfaces. The current validation logic requires a candidate to be enclosed in a local horizontal section, connected to the exterior when the full height is allowed, and not already leaking when movement is restricted to the seed elevation.

#### 4. First-spill search

For each candidate water level, movement is allowed only through free voxels at or below that level. A 6-connected flood fill starts from the validated cavity seed.

Connectivity is monotonic as the permitted maximum Z increases, so a binary search identifies:

- the highest non-leaking level (`safe level`), and
- the first leaking level (`spill level`).

#### 5. Volume

At the safe level, if `N` connected fluid voxels are retained at voxel pitch `p` mm:

```text
V = N * p^3  [mm^3]
```

Results are also reported in mL and L.

### QA and audit outputs

The workflow records, among other diagnostics:

- vertex/face counts before and after exact-coordinate welding
- topology summaries and non-manifold/boundary diagnostics
- voxel grid dimensions and surface-voxel counts
- cavity-seed candidates and selected seed
- safe and spill levels
- fluid-surface and spill-region QC point clouds
- per-run JSON results
- multi-pitch CSV/JSON summaries

QC point clouds can be overlaid on the source mesh in software such as CloudCompare.

### Validation summary

Validation dataset:

**Human-face ink-inscribed pottery, Nara period**\
Goten / Ninomiya Site, Iwata City, Shizuoka, Japan\
Collection: Iwata City Archaeological Center\
Legashizu3D model LS0015\
Source: https://lega-shizu.com/legashizu3d/archives/data/117

Input used for validation: `0015Jinmen_small.ply`, oriented with +Z upward.

| voxel pitch | safe level | first-spill upper level | fluid voxels | capacity | total time |
|---:|---:|---:|---:|---:|---:|
| 2.0 mm | 240.0 mm | 242.0 mm | 725,340 | **5.802720 L** | 14.84 s |
| 1.0 mm | 241.0 mm | 242.0 mm | 5,917,481 | **5.917481 L** | 33.25 s |
| 0.5 mm | 240.5 mm | 241.0 mm | 47,583,957 | **5.947995 L** | 132.08 s |

Observed refinement:

```text
5.802720 -> 5.917481 -> 5.947995 L
```

The 1.0 mm to 0.5 mm difference is **30.514 mL**, or approximately **0.513% of the 0.5 mm result**. The spill upper level varied by 1 mm across the three resolutions.

A three-point empirical convergence diagnostic gave an estimated order of approximately 1.91 and an extrapolated value of approximately 5.959 L. This extrapolation is a diagnostic model value, **not a measured capacity**. The finest directly computed value, 5.947995 L, is the primary detailed result for this dataset.

### Texture-UV seam note

The validation audit showed that Trimesh's normal texture handling can re-index vertices along face-corner UV discontinuities. In this dataset that produced apparent pre-weld boundary edges that could be mistaken for archaeological sherd joins.

For v1.3.1, `archaeological/raw_*` derivatives and raw-boundary/spill proximity diagnostics are therefore disabled. This change does **not** alter the voxel-fluid capacity algorithm. Post-weld topology outputs remain QA products and must not be interpreted as automatic sherd identification.

### Output structure

Typical validation output:

```text
<stem>_PotteryVolume_v1/
├── archaeological/
│   ├── after_exact_weld_boundary_stats.json
│   ├── after_exact_weld_components.csv
│   ├── after_exact_weld_components.json
│   └── after_exact_weld_components_colored.ply
├── processed/
├── qa/
├── pitch_2mm/
├── pitch_1mm/
├── pitch_0p5mm/
├── validation_summary.csv
└── validation_summary.json
```

### Limitations

- Capacity depends on voxel pitch; convergence should be checked on representative data.
- Input orientation directly controls the gravity/spill direction.
- Large geometric gaps may cause early leakage or failure to identify a stable cavity.
- Very fine pitch values can require substantial memory because the voxel grid is three-dimensional.
- The method estimates capacity of the represented geometry; it does not infer missing original vessel geometry.
- Archaeological sherd-boundary detection is outside the current capacity algorithm.

### Relationship to the section-based tools

```text
Voxel-fluid family
3-D mesh
  -> surface voxels / free-space topology
  -> first-spill connectivity
  -> retained voxel volume

Section/profile family
3-D mesh or 2-D drawing
  -> sectional profiles / areas
  -> numerical integration
```

Agreement between the two families is treated as validation evidence, not as a calibration objective.
