"""Full-resolution, rigid sherd placement checks in millimetres."""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from itertools import combinations
from pathlib import Path

import numpy as np

DEFAULT_TOLERANCE_MM = 0.05


def move(points, matrix):
    return points @ matrix[:3, :3].T + matrix[:3, 3]


@dataclass
class SherdMesh:
    mesh: object
    source: Path

    @cached_property
    def solid(self):
        return bool(self.mesh.is_volume)

    @cached_property
    def face_bounds(self):
        tris = self.mesh.triangles
        return np.column_stack((tris.min(axis=1), tris.max(axis=1)))


def load_sherds(paths, scale):
    """Separate files define sherd IDs. A single mesh can use connected parts."""
    import trimesh

    paths = [Path(p) for p in paths]
    sherds = {}
    for path in paths:
        mesh = trimesh.load_scene(path, process=False).to_geometry()
        if not isinstance(mesh, trimesh.Trimesh) or not len(mesh.faces):
            raise ValueError(f"Input must contain triangle faces: {path}")
        if not np.isfinite(mesh.vertices).all():
            raise ValueError(f"Input has non-finite vertex coordinates: {path}")
        mesh.apply_scale(scale)
        mesh.update_faces(mesh.nondegenerate_faces())
        mesh.remove_unreferenced_vertices()
        if not len(mesh.faces):
            raise ValueError(f"Input contains no nondegenerate triangle faces: {path}")
        parts = mesh.split(only_watertight=False, repair=False) if len(paths) == 1 else [mesh]
        for part in parts:
            # Do not weld across distinct sherds. Coincident UV/normal seam
            # vertices in one identified sherd must not hide a closed volume.
            part.merge_vertices(merge_tex=True, merge_norm=True, digits_vertex=9)
            if len(part.faces):
                sherds[len(sherds) + 1] = SherdMesh(part, path)
    if not sherds:
        raise ValueError("No usable triangle mesh in the sherd input.")
    return sherds


def _intersection_intervals(triangles, distances, direction):
    # Intersect the edges with the other plane, and project the resulting
    # segment onto the shared line. Vertices on the plane are included.
    following = np.roll(distances, -1, axis=1)
    edge = np.roll(triangles, -1, axis=1) - triangles
    denom = distances - following
    valid = ((distances <= 0) & (following >= 0) |
             (distances >= 0) & (following <= 0)) & (denom != 0)
    fraction = np.divide(distances, denom, out=np.zeros_like(distances), where=valid)
    projected = np.einsum("nkj,nj->nk", triangles + fraction[..., None] * edge, direction)
    return (np.where(valid, projected, np.inf).min(axis=1),
            np.where(valid, projected, -np.inf).max(axis=1))


def _proper_crossings(a, b, tolerance):
    an = np.cross(a[:, 1] - a[:, 0], a[:, 2] - a[:, 0])
    bn = np.cross(b[:, 1] - b[:, 0], b[:, 2] - b[:, 0])
    an /= np.linalg.norm(an, axis=1)[:, None]
    bn /= np.linalg.norm(bn, axis=1)[:, None]
    da = np.einsum("nkj,nj->nk", a - b[:, :1], bn)
    db = np.einsum("nkj,nj->nk", b - a[:, :1], an)
    keep = ((da.min(axis=1) < -tolerance) & (da.max(axis=1) > tolerance) &
            (db.min(axis=1) < -tolerance) & (db.max(axis=1) > tolerance))
    if not keep.any():
        return False
    direction = np.cross(an[keep], bn[keep])
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    alo, ahi = _intersection_intervals(a[keep], da[keep], direction)
    blo, bhi = _intersection_intervals(b[keep], db[keep], direction)
    return bool(np.any(np.minimum(ahi, bhi) - np.maximum(alo, blo) > tolerance))


def _inside(probe, target, transform, tolerance):
    import trimesh

    # Vertices catch containment and many aligned penetrations. Inward face
    # probes distinguish coincident solids from opposing, touching faces.
    bounds = target.mesh.bounds
    candidates = [probe.mesh.vertices]
    if probe.solid:
        candidates.append(probe.mesh.triangles_center -
                          probe.mesh.face_normals * (2 * tolerance))
    for source in candidates:
        points = move(source, transform)
        points = points[np.all((points > bounds[0] + tolerance) &
                               (points < bounds[1] - tolerance), axis=1)]
        for start in range(0, len(points), 256):
            distance = trimesh.proximity.signed_distance(target.mesh, points[start:start + 256])
            if np.any(distance > tolerance):
                return True
    return False


def penetration(a, b, b_to_a, tolerance=DEFAULT_TOLERANCE_MM):
    """Return evidence type or None. Tangency and coplanar contact are allowed."""
    bounds = a.mesh.bounds
    bbox = b.mesh.bounds
    centre = move(bbox.mean(axis=0), b_to_a)
    radius = np.abs(b_to_a[:3, :3]) @ (np.ptp(bbox, axis=0) / 2)
    overlap = np.minimum(bounds[1], centre + radius) - np.maximum(bounds[0], centre - radius)
    if np.any(overlap < -tolerance) or (np.any(overlap <= tolerance) and a.solid and b.solid):
        return None
    points = move(b.mesh.vertices, b_to_a)
    indices = np.arange(len(b.mesh.faces))
    if np.allclose(b_to_a[:3, :3], np.eye(3), atol=1e-12, rtol=0):
        lo, hi = np.split(b.face_bounds, 2, axis=1)
        lo, hi = lo + b_to_a[:3, 3], hi + b_to_a[:3, 3]
        indices = indices[np.all((hi >= bounds[0]) & (lo <= bounds[1]), axis=1)]
    tree = None
    for start in range(0, len(indices), 256):
        tris = points[b.mesh.faces[indices[start:start + 256]]]
        lows, highs = tris.min(axis=1), tris.max(axis=1)
        tri_ids, other_ids = [], []
        for idx in np.flatnonzero(np.all((highs >= bounds[0]) & (lows <= bounds[1]), axis=1)):
            if tree is None:
                tree = a.mesh.triangles_tree
            matches = list(tree.intersection((*lows[idx], *highs[idx])))
            other_ids.extend(matches)
            tri_ids.extend([idx] * len(matches))
        for j in range(0, len(tri_ids), 2048):
            if _proper_crossings(a.mesh.triangles[other_ids[j:j + 2048]],
                                 tris[tri_ids[j:j + 2048]], tolerance):
                return "surface_crossing"
    if a.solid and _inside(b, a, b_to_a, tolerance):
        return "solid_interior"
    if b.solid and _inside(a, b, np.linalg.inv(b_to_a), tolerance):
        return "solid_interior"
    return None


class Placement:
    def __init__(self, sherds, tolerance=DEFAULT_TOLERANCE_MM):
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("Collision tolerance must be finite and > 0 mm.")
        self.sherds = sherds
        self.tolerance = tolerance
        self.poses = {fid: np.eye(4) for fid in sherds}

    def pair(self, a, b, pose_a=None, pose_b=None):
        relative = np.linalg.solve(self.poses[a] if pose_a is None else pose_a,
                                   self.poses[b] if pose_b is None else pose_b)
        return penetration(self.sherds[a], self.sherds[b], relative, self.tolerance)

    def conflicts(self, fid, candidate):
        return [other for other in self.sherds if other != fid and (
            self.pair(other, fid, pose_b=candidate) if other < fid
            else self.pair(fid, other, pose_a=candidate))]

    def audit(self):
        pairs = []
        for a, b in combinations(sorted(self.sherds), 2):
            evidence = self.pair(a, b)
            pairs.append({"fragment_a": a, "fragment_b": b,
                          "status": "penetration" if evidence else "clear",
                          "evidence": evidence,
                          "surface_only": not (self.sherds[a].solid and self.sherds[b].solid)})
        conflicts = sorted({p[key] for p in pairs if p["status"] == "penetration"
                            for key in ("fragment_a", "fragment_b")})
        status = ("penetration" if conflicts else "not_applicable" if not pairs else
                  "clear_surface_only" if any(p["surface_only"] for p in pairs) else "clear")
        return {"status": status, "tolerance_mm": self.tolerance, "conflicting_ids": conflicts,
                "fragments": [{"fragment_id": fid, "source": str(s.source),
                               "triangle_count": len(s.mesh.faces), "closed_solid": s.solid}
                              for fid, s in self.sherds.items()],
                "pairs": pairs}
