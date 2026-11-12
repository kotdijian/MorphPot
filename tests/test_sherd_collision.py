import json
from pathlib import Path
from dataclasses import replace
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import trimesh

import pottery_guidefit_reconstruction as guide
from morphpot.sherd_collision import Placement, SherdMesh, load_sherds


def box(center, size=(10, 10, 10)):
    return trimesh.creation.box(extents=size).apply_translation(center)


def scene(*meshes):
    return Placement({i: SherdMesh(mesh, Path(f"{i}.ply")) for i, mesh in enumerate(meshes, 1)})


@pytest.mark.parametrize(("center", "size", "status"), [
    ((12, 0, 0), (10, 10, 10), "clear"),
    ((10, 0, 0), (10, 10, 10), "clear"),
    ((9.98, 0, 0), (10, 10, 10), "clear"),
    ((9, 0, 0), (10, 10, 10), "penetration"),
    ((0, 0, 0), (2, 2, 2), "penetration"),
    ((0, 0, 0), (10, 10, 10), "penetration"),
])
def test_solid_contact_overlap_and_containment(center, size, status):
    assert scene(box((0, 0, 0)), box(center, size)).audit()["status"] == status


def test_open_triangle_crossing_is_not_coplanar_contact():
    def triangle(points):
        return trimesh.Trimesh(vertices=points, faces=[[0, 1, 2]], process=False)
    a = triangle([[-3, -2, 0], [4, -2, 0], [0, 4, 0]])
    b = triangle([[0, -3, -3], [0, 3, -3], [0, 0, 4]])
    assert scene(a, b).pair(1, 2) == "surface_crossing"
    assert scene(a, a.copy()).audit()["status"] == "clear_surface_only"


def test_touching_files_keep_identity_and_cli_reports_initial_overlap(tmp_path):
    a, b, joined = [tmp_path / f"{name}.ply" for name in ("a", "b", "joined")]
    ma = box((0, 0, 0))
    mb = box((9, 0, 0))
    ma.export(a)
    mb.export(b)
    trimesh.util.concatenate([ma, mb]).export(joined)
    assert len(load_sherds([joined], 1)) == 2
    target = tmp_path / "report"
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(root / "pottery_reconstruction_3d.py"), str(a), str(b),
         "--collision-only", "--input-unit", "mm", "--output-dir", str(target)],
        capture_output=True, text=True,
    )
    assert result.returncode == 2, result.stderr
    report = json.loads((target / "sherd_collision_report.json").read_text())
    assert report["conflicting_ids"] == [1, 2]
    assert [f["source"] for f in report["fragments"]] == [str(a), str(b)]
    assert np.allclose(trimesh.load_mesh(b).bounds, mb.bounds)


def test_guidefit_rejects_penetrating_best_fit_and_uses_safe_shorter_move(monkeypatch, tmp_path):
    # A sherd spanning X=28..32 would prefer X=20..24. A fixed neighbour
    # spanning X=19..23 blocks that fit but not every improving move.
    moving = box((30, 0, 0), (4, 4, 4))
    blocker = box((21, 0, 0), (4, 4, 4))
    assembly = scene(moving, blocker)
    axis = SimpleNamespace(axis_point_mm=np.zeros(3), axis_direction=np.array([0., 0., 1.]))
    profile = SimpleNamespace(z_mm=np.array([-10., 10.]),
                              outer_model_mm=np.full(2, 24.), inner_model_mm=np.full(2, 20.))
    vessel = SimpleNamespace(axis=axis, profile=profile)
    z = np.linspace(-2, 2, 40)
    outer = np.column_stack((np.full(40, 32.), np.zeros(40), z))
    inner = np.column_stack((np.full(40, 28.), np.zeros(40), z))
    empty = np.empty((0, 3))
    clusters = [SimpleNamespace(fragment_id=i) for i in (1, 2)]
    result = SimpleNamespace(axis=axis, input_unit="mm", sherd_meshes=assembly.sherds,
                             clusters=clusters, metadata={"sherd_collision": assembly.audit()})
    monkeypatch.setattr(guide, "extract_fragment_outer_inner_samples",
                        lambda result, fid, **kw: (
                            guide.SurfaceSamples(outer, inner, empty, empty, 40, 40) if fid == 1
                            else guide.SurfaceSamples(empty, empty, empty, empty, 0, 0)))
    unconstrained = guide.fit_fragment_to_guide_vessel(result, vessel, 1, tilt_limit_deg=0,
                                                      radial_search_limit_mm=12)
    assert unconstrained.radial_shift_mm == pytest.approx(-8, abs=1e-3)
    assert assembly.conflicts(1, unconstrained.matrix_mm) == [2]

    fits = guide.fit_all_fragments(result, vessel, tilt_limit_deg=0, radial_search_limit_mm=12)
    fit = fits[1]
    assert fit.status == "ok"
    assert fits[2].status == "insufficient_dual_surface"
    assert fit.collision_blocked_by == (2,)
    assert -8 < fit.radial_shift_mm < 0
    assert fit.combined_after_mm < fit.combined_before_mm
    assembly.poses[1] = fit.matrix_mm
    assert assembly.audit()["status"] == "clear"

    result.export_unit = "mm"
    result.mesh_vertices_mm = vessel.mesh_vertices_mm = moving.vertices
    result.mesh_faces = vessel.mesh_faces = moving.faces
    result.display_points_mm = np.vstack([moving.vertices, blocker.vertices])
    result.display_labels = np.repeat([1, 2], [len(moving.vertices), len(blocker.vertices)])
    profile.inner_observed_mm = profile.inner_model_mm
    profile.outer_observed_mm = profile.outer_model_mm
    profile.thickness_model_mm = profile.outer_model_mm - profile.inner_model_mm
    vessel.current_residual_mm = np.zeros(len(result.display_labels))
    vessel.current_nearest_surface = np.zeros(len(result.display_labels), dtype=np.uint8)
    vessel.source_ring_count, vessel.outer_ring_count = 2, 1
    vessel.inner_ring_count, vessel.unknown_ring_count = 1, 0
    reg = guide.base.GuideRegistration("ok", 1, np.eye(3), np.zeros(3), np.eye(4),
                                      0, 0, 0, 0, 0, 0, 2, 2, 20, -10, 10, -10, 10, 0, 0)
    session = guide.GuideFitSession(result, None, reg, vessel, fits)
    target = tmp_path / "export"
    guide.export_guidefit_session(session, target)
    meshes = sorted((target / "fitted_sherd_meshes").glob("*.ply"))
    assert len(meshes) == 2
    assert Placement(load_sherds(meshes, 1)).audit()["status"] == "clear"
    assert json.loads((target / "sherd_collision_report.json").read_text())["fitted"]["status"] == "clear"
    unsafe = replace(session, sherd_fits={1: unconstrained, 2: fits[2]})
    with pytest.raises(ValueError, match="書き出しません"):
        guide.export_guidefit_session(unsafe, tmp_path / "unsafe")
    assert not (tmp_path / "unsafe").exists()

    touching = scene(moving, box((22, 0, 0), (12, 4, 4)))
    result.sherd_meshes = touching.sherds
    result.metadata["sherd_collision"] = touching.audit()
    blocked = guide.fit_all_fragments(result, vessel, tilt_limit_deg=0, radial_search_limit_mm=12)[1]
    assert blocked.status == "blocked_collision"
    assert np.array_equal(blocked.matrix_mm, np.eye(4))
    result.metadata["sherd_collision"] = scene(moving, moving.copy()).audit()
    with pytest.raises(ValueError, match="初期配置"):
        guide.fit_all_fragments(result, vessel)
