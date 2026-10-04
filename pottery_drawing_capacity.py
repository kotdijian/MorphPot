#!/usr/bin/env python3
"""
PotteryDrawingCapacity
Version 0.6.0

Interactive GUI for calculating pottery capacity from a 2-D longitudinal-section drawing.

Assumptions for v0.6.0:
- The drawing is already upright: image horizontal = vessel horizontal, image vertical = vessel vertical.
- Both left and right inner-wall profiles are drawn.
- The user calibrates scale by clicking two points of known real-world distance.
- The user specifies the vertical rotation axis by one click (X position only).
- The user digitizes the left and right inner profiles manually with ordered clicks.
- Capacity uses the same two-sided 'single longitudinal section' definition as the 3-D tool:

      A(z) = pi/2 * (r_left(z)^2 + r_right(z)^2)

  which is equivalent to rotating the area-preserving equivalent radius

      r_eq(z) = sqrt((r_left(z)^2 + r_right(z)^2)/2).

Supported interactive source formats:
- Raster images: PNG / JPEG / TIFF / BMP / WEBP
- SVG (including Adobe Illustrator exported SVG): rendered with PyMuPDF for digitizing
- PDF: selected page rendered with PyMuPDF for digitizing

SVG/PDF vector semantics are intentionally not trusted for scale; GUI calibration is always used.
All calculation inputs are exported to CSV/JSON for audit and reproducibility.
Numerical capacity integration is delegated to pottery_volume_core.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pottery_volume_core as volume_core
from PIL import Image, ImageDraw, ImageFont, ImageTk

__version__ = "0.6.0"

RASTER_EXTENSIONS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
VECTOR_DOCUMENT_EXTENSIONS = {".svg", ".pdf"}
SUPPORTED_EXTENSIONS = RASTER_EXTENSIONS | VECTOR_DOCUMENT_EXTENSIONS


@dataclass
class Calibration:
    p1_px: tuple[float, float]
    p2_px: tuple[float, float]
    known_length_mm: float
    pixel_length: float
    mm_per_pixel: float


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_source_image(path: Path, page_number: int = 1, render_dpi: float = 180.0) -> tuple[Image.Image, dict]:
    """Load a source into a pixel image used by the digitizing GUI."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported drawing format: {ext}. Supported: {sorted(SUPPORTED_EXTENSIONS)}")

    if ext in RASTER_EXTENSIONS:
        img = Image.open(path).convert("RGB")
        return img, {
            "source_kind": "raster_image",
            "source_extension": ext,
            "page_number": None,
            "render_dpi": None,
            "vector_semantics_used": False,
        }

    try:
        import pymupdf
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "SVG/PDF input requires PyMuPDF. Install with: python3 -m pip install pymupdf"
        ) from exc

    doc = pymupdf.open(str(path))
    if doc.page_count < 1:
        raise ValueError("SVG/PDF contains no renderable page.")
    page_number = int(page_number)
    if not (1 <= page_number <= doc.page_count):
        raise ValueError(f"page_number must be 1..{doc.page_count}")
    page = doc[page_number - 1]
    zoom = float(render_dpi) / 72.0
    pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    kind = "svg_render" if ext == ".svg" else "pdf_page_render"
    return img, {
        "source_kind": kind,
        "source_extension": ext,
        "page_number": page_number,
        "document_page_count": int(doc.page_count),
        "render_dpi": float(render_dpi),
        "vector_semantics_used": False,
        "note": "The rendered pixels are digitized. Physical scale is determined only by GUI calibration.",
    }


def _dedupe_profile_z(z_mm: np.ndarray, r_mm: np.ndarray, tolerance_mm: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(z_mm)
    z = np.asarray(z_mm, dtype=float)[order]
    r = np.asarray(r_mm, dtype=float)[order]
    out_z: list[float] = []
    out_r: list[float] = []
    i = 0
    while i < len(z):
        j = i + 1
        while j < len(z) and abs(float(z[j] - z[i])) <= tolerance_mm:
            j += 1
        out_z.append(float(np.mean(z[i:j])))
        out_r.append(float(np.median(r[i:j])))
        i = j
    return np.asarray(out_z, dtype=float), np.asarray(out_r, dtype=float)


def profile_from_pixels(
    points_px: Iterable[tuple[float, float]],
    side: str,
    axis_x_px: float,
    mm_per_pixel: float,
    y_bottom_px: float,
) -> dict:
    pts = np.asarray(list(points_px), dtype=float)
    if len(pts) < 2:
        raise ValueError(f"{side} profile requires at least 2 points.")
    x = pts[:, 0]
    y = pts[:, 1]
    z_mm = (float(y_bottom_px) - y) * float(mm_per_pixel)
    if side == "left":
        r_mm = (float(axis_x_px) - x) * float(mm_per_pixel)
    elif side == "right":
        r_mm = (x - float(axis_x_px)) * float(mm_per_pixel)
    else:
        raise ValueError("side must be left or right")
    if np.any(r_mm < -1e-6):
        raise ValueError(f"Some {side} profile points cross the specified rotation axis.")
    r_mm = np.maximum(r_mm, 0.0)
    z_u, r_u = _dedupe_profile_z(z_mm, r_mm)
    if len(z_u) < 2 or np.ptp(z_u) <= 0:
        raise ValueError(f"{side} profile does not span a measurable vertical range.")
    return {
        "side": side,
        "points_px": pts,
        "z_mm": z_u,
        "radius_mm": r_u,
    }


def calculate_two_sided_single(
    left: dict,
    right: dict,
    z_step_mm: float = 0.5,
) -> dict:
    """Delegate Drawing-Single numerical integration to PotteryVolumeCore."""
    calc = volume_core.calculate_single_two_sided(
        left["z_mm"], left["radius_mm"],
        right["z_mm"], right["radius_mm"],
        z_step_mm=z_step_mm,
        method="single",
    )
    # The GUI/export layer keeps the historical flat dictionary interface,
    # while all numerical integration is performed in pottery_volume_core.py.
    return {**calc["summary"], **calc["profile"]}

def write_profile_csv(path: Path, left_points: list[tuple[float, float]], right_points: list[tuple[float, float]],
                      axis_x_px: float, mm_per_pixel: float, y_bottom_px: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["side", "point_order", "x_source_px", "y_source_px", "signed_x_from_axis_mm", "radius_mm", "z_mm"])
        for side, pts in (("left", left_points), ("right", right_points)):
            for i, (x, y) in enumerate(pts, start=1):
                signed = (x - axis_x_px) * mm_per_pixel
                radius = -signed if side == "left" else signed
                z = (y_bottom_px - y) * mm_per_pixel
                w.writerow([side, i, x, y, signed, radius, z])


def write_resampled_csv(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["z_mm", "left_radius_mm", "right_radius_mm", "equivalent_radius_mm", "area_mm2"])
        for row in zip(result["z_mm"], result["left_radius_mm"], result["right_radius_mm"],
                       result["equivalent_radius_mm"], result["area_mm2"]):
            w.writerow([float(v) for v in row])


def jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


def save_qc_overlay(source_image: Image.Image, path: Path, calibration: Calibration, axis_x_px: float,
                    left_points, right_points, y_bottom_px: float, spill_z_mm: float, mm_per_pixel: float) -> None:
    im = source_image.copy().convert("RGB")
    draw = ImageDraw.Draw(im)
    width = max(2, int(round(max(im.size) / 900)))
    r = max(3, 2 * width)

    # calibration
    draw.line([calibration.p1_px, calibration.p2_px], fill=(180, 40, 40), width=width)
    for x, y in (calibration.p1_px, calibration.p2_px):
        draw.ellipse((x-r, y-r, x+r, y+r), outline=(180, 40, 40), width=width)
    # vertical axis
    draw.line([(axis_x_px, 0), (axis_x_px, im.height)], fill=(30, 100, 210), width=width)
    # profiles
    if len(left_points) >= 2:
        draw.line(left_points, fill=(20, 150, 80), width=max(width, 2))
    if len(right_points) >= 2:
        draw.line(right_points, fill=(190, 120, 20), width=max(width, 2))
    for pts, col in ((left_points, (20,150,80)), (right_points, (190,120,20))):
        for x, y in pts:
            draw.ellipse((x-r, y-r, x+r, y+r), fill=col)
    # spill line
    spill_y = y_bottom_px - spill_z_mm / mm_per_pixel
    draw.line([(0, spill_y), (im.width, spill_y)], fill=(120, 40, 160), width=width)
    try:
        font = ImageFont.load_default()
        draw.text((10, 10), f"scale: {mm_per_pixel:.6g} mm/px", fill=(0,0,0), font=font)
        draw.text((10, 28), f"spill: {spill_z_mm:.3f} mm", fill=(0,0,0), font=font)
    except Exception:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    im.save(path)


def save_profile_plot(path: Path, result: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    z = np.asarray(result["z_mm"], float)
    l = np.asarray(result["left_radius_mm"], float)
    r = np.asarray(result["right_radius_mm"], float)
    fig, ax = plt.subplots(figsize=(7.0, 8.0))
    ax.plot(-l, z, label="left inner")
    ax.plot(r, z, label="right inner")
    ax.axvline(0.0, linestyle="--", linewidth=1.0, label="rotation axis")
    ax.axhline(float(result["spill_z_mm"]), linestyle=":", linewidth=1.0, label="spill level")
    ax.set_xlabel("Radius from axis [mm]")
    ax.set_ylabel("Z [mm]")
    ax.set_title("Digitized longitudinal inner profiles")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def export_drawing_result(output_dir: Path, source_path: Path, source_image: Image.Image, source_meta: dict,
                          calibration: Calibration, axis_x_px: float, left_points, right_points,
                          z_step_mm: float) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_image.save(output_dir / "source_render.png")

    all_y = [p[1] for p in left_points] + [p[1] for p in right_points]
    if not all_y:
        raise ValueError("No profile points.")
    y_bottom_px = float(max(all_y))
    left = profile_from_pixels(left_points, "left", axis_x_px, calibration.mm_per_pixel, y_bottom_px)
    right = profile_from_pixels(right_points, "right", axis_x_px, calibration.mm_per_pixel, y_bottom_px)
    result = calculate_two_sided_single(left, right, z_step_mm=z_step_mm)

    write_profile_csv(output_dir / "drawing_profile_raw.csv", list(left_points), list(right_points),
                      axis_x_px, calibration.mm_per_pixel, y_bottom_px)
    write_resampled_csv(output_dir / "drawing_profile_resampled.csv", result)

    summary_fields = [
        "method", "status", "core_version", "bottom_z_mm", "spill_z_mm", "common_height_mm", "z_step_mm",
        "volume_mm3", "volume_ml", "volume_l", "left_only_volume_l", "right_only_volume_l",
        "left_right_difference_l", "left_right_difference_percent_of_mean", "integration",
        "area_definition", "equivalent_radius_definition",
    ]
    with (output_dir / "drawing_volume_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerow({k: result.get(k) for k in summary_fields})

    metadata = {
        "program": "PotteryDrawingCapacity",
        "version": __version__,
        "volume_core_version": volume_core.__version__,
        "calculation_engine": "pottery_volume_core.calculate_single_two_sided",
        "source_file": str(source_path),
        "source_sha256": sha256_file(source_path),
        "source": source_meta,
        "rendered_image_size_px": [int(source_image.width), int(source_image.height)],
        "assumptions": {
            "drawing_horizontal_vertical_fixed": True,
            "both_left_and_right_inner_profiles_required": True,
            "rotation_axis_vertical": True,
            "axis_x_manual_click": True,
            "scale_gui_calibration": True,
        },
        "calibration": asdict(calibration),
        "axis_x_px": float(axis_x_px),
        "profile_z_zero_definition": "lowest digitized profile point in rendered image",
        "spill_definition": "lower of the two digitized inner-profile top elevations (common upper Z limit)",
        "calculation": {k: jsonable(v) for k, v in result.items() if k not in {"z_mm","left_radius_mm","right_radius_mm","equivalent_radius_mm","area_mm2"}},
    }
    (output_dir / "drawing_metadata.json").write_text(json.dumps(jsonable(metadata), ensure_ascii=False, indent=2), encoding="utf-8")
    (output_dir / "drawing_volume_summary.json").write_text(json.dumps(jsonable(metadata["calculation"]), ensure_ascii=False, indent=2), encoding="utf-8")

    save_qc_overlay(source_image, output_dir / "drawing_qc_overlay.png", calibration, axis_x_px,
                    list(left_points), list(right_points), y_bottom_px, float(result["spill_z_mm"]), calibration.mm_per_pixel)
    save_profile_plot(output_dir / "drawing_profile_plot.png", result)
    return {"result": result, "metadata": metadata, "output_dir": output_dir}


class DrawingCapacityGUI:
    def __init__(self, root, source_path: Path | None, page_number: int, render_dpi: float,
                 z_step_mm: float, output_dir: Path | None):
        import tkinter as tk
        from tkinter import filedialog, messagebox, simpledialog, ttk

        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.simpledialog = simpledialog
        self.ttk = ttk
        self.root = root
        self.root.title(f"PotteryDrawingCapacity v{__version__}")
        self.root.geometry("1240x860")

        self.source_path = Path(source_path) if source_path else None
        self.page_number = int(page_number)
        self.render_dpi = float(render_dpi)
        self.z_step_mm = float(z_step_mm)
        self.output_dir_arg = Path(output_dir) if output_dir else None

        self.source_image: Image.Image | None = None
        self.source_meta: dict = {}
        self.tk_image = None
        self.zoom = 1.0
        self.mode: str | None = None
        self.scale_points: list[tuple[float,float]] = []
        self.calibration: Calibration | None = None
        self.axis_x_px: float | None = None
        self.left_points: list[tuple[float,float]] = []
        self.right_points: list[tuple[float,float]] = []

        self.status_var = tk.StringVar(value="図面を開いてください。")
        self.scale_var = tk.StringVar(value="scale: 未設定")
        self.axis_var = tk.StringVar(value="axis: 未設定")

        top = ttk.Frame(root)
        top.pack(side="top", fill="x")
        for text, cmd in [
            ("図面を開く", self.open_source_dialog),
            ("1 縮尺 (2点)", self.start_scale),
            ("2 回転軸X", self.start_axis),
            ("3 左内面", lambda: self.start_profile("left")),
            ("4 右内面", lambda: self.start_profile("right")),
            ("1点戻す", self.undo),
            ("計算・保存", self.calculate_and_save),
        ]:
            ttk.Button(top, text=text, command=cmd).pack(side="left", padx=3, pady=4)
        ttk.Button(top, text="−", width=3, command=lambda: self.change_zoom(0.8)).pack(side="left", padx=2)
        ttk.Button(top, text="＋", width=3, command=lambda: self.change_zoom(1.25)).pack(side="left", padx=2)
        ttk.Label(top, textvariable=self.scale_var).pack(side="left", padx=8)
        ttk.Label(top, textvariable=self.axis_var).pack(side="left", padx=8)

        body = ttk.Frame(root)
        body.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(body, background="white")
        xbar = ttk.Scrollbar(body, orient="horizontal", command=self.canvas.xview)
        ybar = ttk.Scrollbar(body, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=xbar.set, yscrollcommand=ybar.set)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        ybar.grid(row=0, column=1, sticky="ns")
        xbar.grid(row=1, column=0, sticky="ew")
        body.rowconfigure(0, weight=1); body.columnconfigure(0, weight=1)

        ttk.Label(root, textvariable=self.status_var, anchor="w").pack(side="bottom", fill="x", padx=6, pady=3)
        self.canvas.bind("<Button-1>", self.on_click)
        self.canvas.bind("<Button-3>", lambda e: self.finish_current())
        self.root.bind("<Return>", lambda e: self.finish_current())
        self.root.bind("<Control-z>", lambda e: self.undo())
        self.canvas.bind("<MouseWheel>", self.on_mousewheel)
        self.canvas.bind("<Button-4>", lambda e: self.change_zoom(1.12))
        self.canvas.bind("<Button-5>", lambda e: self.change_zoom(0.89))

        if self.source_path:
            self.load_source(self.source_path)

    def open_source_dialog(self):
        filename = self.filedialog.askopenfilename(
            title="土器縦断面図を選択",
            filetypes=[
                ("対応図面", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.webp *.svg *.pdf"),
                ("画像", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.webp"),
                ("SVG", "*.svg"), ("PDF", "*.pdf"), ("すべて", "*.*")
            ])
        if filename:
            self.load_source(Path(filename))

    def load_source(self, path: Path):
        page = self.page_number
        if path.suffix.lower() == ".pdf":
            try:
                import pymupdf
                d = pymupdf.open(str(path))
                if d.page_count > 1:
                    val = self.simpledialog.askinteger("PDFページ", f"ページ番号 (1-{d.page_count})", initialvalue=min(page,d.page_count), minvalue=1, maxvalue=d.page_count)
                    if val is None:
                        return
                    page = int(val)
            except Exception:
                pass
        try:
            image, meta = load_source_image(path, page_number=page, render_dpi=self.render_dpi)
        except Exception as exc:
            self.messagebox.showerror("読み込みエラー", str(exc)); return
        self.source_path = path
        self.page_number = page
        self.source_image = image
        self.source_meta = meta
        self.scale_points=[]; self.calibration=None; self.axis_x_px=None; self.left_points=[]; self.right_points=[]
        self.scale_var.set("scale: 未設定"); self.axis_var.set("axis: 未設定")
        # Fit initial view to a typical 1100 x 700 viewport, never enlarge above 1x initially.
        self.zoom = min(1.0, 1080/max(1,image.width), 690/max(1,image.height))
        self.zoom = max(self.zoom, 0.05)
        self.redraw()
        self.status_var.set("1) 縮尺を指定: 『1 縮尺 (2点)』を押し、既知距離の両端をクリック。")

    def display_xy_to_source(self, event) -> tuple[float,float]:
        x = self.canvas.canvasx(event.x) / self.zoom
        y = self.canvas.canvasy(event.y) / self.zoom
        return float(x), float(y)

    def source_to_display(self, p):
        return float(p[0]*self.zoom), float(p[1]*self.zoom)

    def on_click(self, event):
        if self.source_image is None:
            return
        x,y = self.display_xy_to_source(event)
        if not (0 <= x < self.source_image.width and 0 <= y < self.source_image.height):
            return
        if self.mode == "scale":
            self.scale_points.append((x,y))
            if len(self.scale_points) == 2:
                p1,p2=self.scale_points
                pix=math.hypot(p2[0]-p1[0],p2[1]-p1[1])
                known=self.simpledialog.askfloat("実寸", "選択した2点間の実寸 [mm]", minvalue=1e-9)
                if known is None or pix <= 0:
                    self.scale_points=[]
                else:
                    self.calibration=Calibration(p1,p2,float(known),float(pix),float(known/pix))
                    self.scale_var.set(f"scale: {self.calibration.mm_per_pixel:.6g} mm/px")
                    self.status_var.set("2) 回転軸を指定: 『2 回転軸X』を押して縦の中心軸上をクリック。")
                self.mode=None
        elif self.mode == "axis":
            self.axis_x_px=x
            self.axis_var.set(f"axis X: {x:.2f} px")
            self.mode=None
            self.status_var.set("3) 左内面を底部から口縁までクリック。Enter/右クリックで終了。次に右内面。")
        elif self.mode == "left":
            self.left_points.append((x,y))
        elif self.mode == "right":
            self.right_points.append((x,y))
        self.redraw()

    def start_scale(self):
        if self.source_image is None: return
        self.mode="scale"; self.scale_points=[]
        self.status_var.set("縮尺: 既知距離の1点目・2点目を順にクリックしてください。")
        self.redraw()

    def start_axis(self):
        if self.source_image is None: return
        self.mode="axis"
        self.status_var.set("回転軸: 縦中心軸上の任意の1点をクリックしてください（X座標のみ使用）。")

    def start_profile(self, side):
        if self.source_image is None: return
        if self.axis_x_px is None:
            self.messagebox.showwarning("回転軸未設定", "先に回転軸Xを指定してください。")
            return
        self.mode=side
        self.status_var.set(f"{ '左' if side=='left' else '右' }内面: 底部から口縁へ順にクリック。Enter/右クリックで終了。Ctrl+Zで1点戻す。")

    def finish_current(self):
        if self.mode in {"left","right"}:
            side=self.mode
            n=len(self.left_points if side=="left" else self.right_points)
            self.status_var.set(f"{ '左' if side=='left' else '右' }内面を{n}点で確定。")
            self.mode=None

    def undo(self):
        target=None
        if self.mode=="left" or (self.mode is None and self.left_points and not self.right_points): target=self.left_points
        elif self.mode=="right" or (self.mode is None and self.right_points): target=self.right_points
        elif self.mode=="scale" and self.scale_points: target=self.scale_points
        if target:
            target.pop(); self.redraw()

    def change_zoom(self, factor):
        if self.source_image is None: return
        self.zoom=max(0.03,min(8.0,self.zoom*float(factor))); self.redraw()

    def on_mousewheel(self, event):
        self.change_zoom(1.12 if event.delta>0 else 0.89)

    def redraw(self):
        self.canvas.delete("all")
        if self.source_image is None:
            return
        w=max(1,int(round(self.source_image.width*self.zoom))); h=max(1,int(round(self.source_image.height*self.zoom)))
        display=self.source_image.resize((w,h),Image.Resampling.LANCZOS)
        self.tk_image=ImageTk.PhotoImage(display)
        self.canvas.create_image(0,0,image=self.tk_image,anchor="nw")
        self.canvas.configure(scrollregion=(0,0,w,h))
        # scale line
        if len(self.scale_points)>=1:
            p=[self.source_to_display(x) for x in self.scale_points]
            for x,y in p: self.canvas.create_oval(x-4,y-4,x+4,y+4,outline="red",width=2)
            if len(p)==2: self.canvas.create_line(*p[0],*p[1],fill="red",width=2)
        # axis
        if self.axis_x_px is not None:
            x=self.axis_x_px*self.zoom; self.canvas.create_line(x,0,x,h,fill="blue",width=2,dash=(6,4))
        # profiles
        for pts,col in ((self.left_points,"green"),(self.right_points,"orange")):
            dp=[self.source_to_display(p) for p in pts]
            if len(dp)>=2:
                flat=[v for p in dp for v in p]; self.canvas.create_line(*flat,fill=col,width=2)
            for x,y in dp: self.canvas.create_oval(x-3,y-3,x+3,y+3,fill=col,outline=col)

    def calculate_and_save(self):
        if self.source_path is None or self.source_image is None:
            self.messagebox.showwarning("未設定", "図面を開いてください。"); return
        if self.calibration is None:
            self.messagebox.showwarning("未設定", "縮尺を指定してください。"); return
        if self.axis_x_px is None:
            self.messagebox.showwarning("未設定", "回転軸Xを指定してください。"); return
        if len(self.left_points)<2 or len(self.right_points)<2:
            self.messagebox.showwarning("未設定", "左右の内面プロファイルを各2点以上指定してください。"); return
        out=self.output_dir_arg or (self.source_path.parent / f"{self.source_path.stem}_DrawingCapacity")
        try:
            exported=export_drawing_result(out,self.source_path,self.source_image,self.source_meta,self.calibration,
                                           self.axis_x_px,self.left_points,self.right_points,self.z_step_mm)
        except Exception as exc:
            self.messagebox.showerror("計算エラー",str(exc)); return
        r=exported["result"]
        self.status_var.set(f"容量 {r['volume_l']:.6f} L - 保存: {out}")
        self.messagebox.showinfo(
            "計算完了",
            f"容量: {r['volume_l']:.6f} L\n"
            f"左のみ: {r['left_only_volume_l']:.6f} L\n"
            f"右のみ: {r['right_only_volume_l']:.6f} L\n"
            f"spill Z: {r['spill_z_mm']:.3f} mm\n\n保存先:\n{out}"
        )


def self_test() -> None:
    """Headless numerical self-test using a 50 mm radius, 100 mm high cylinder."""
    mm_per_px = 0.5
    axis_x = 200.0
    y_bottom = 300.0
    # z=0 at y=300, z=100 at y=100; radius=50 mm = 100 px
    left_pts=[(100.0,300.0),(100.0,250.0),(100.0,200.0),(100.0,150.0),(100.0,100.0)]
    right_pts=[(300.0,300.0),(300.0,250.0),(300.0,200.0),(300.0,150.0),(300.0,100.0)]
    left=profile_from_pixels(left_pts,"left",axis_x,mm_per_px,y_bottom)
    right=profile_from_pixels(right_pts,"right",axis_x,mm_per_px,y_bottom)
    result=calculate_two_sided_single(left,right,z_step_mm=0.5)
    expected=math.pi*50.0**2*100.0/1_000_000.0
    err=abs(result["volume_l"]-expected)
    print(f"expected : {expected:.9f} L")
    print(f"computed : {result['volume_l']:.9f} L")
    print(f"error    : {err:.12g} L")
    if err > 1e-9:
        raise SystemExit("SELF TEST FAILED")
    print("SELF TEST PASSED")


def main():
    parser=argparse.ArgumentParser(description="GUI capacity calculation from a 2-D pottery longitudinal-section drawing.")
    parser.add_argument("input",nargs="?",type=Path,help="PNG/JPEG/TIFF/BMP/WEBP/SVG/PDF drawing. Omit to choose in GUI.")
    parser.add_argument("--page",type=int,default=1,help="PDF page number, 1-based (default: 1; GUI asks when multipage PDF is opened)")
    parser.add_argument("--render-dpi",type=float,default=180.0,help="SVG/PDF rendering DPI for GUI digitizing (default: 180)")
    parser.add_argument("--z-step-mm",type=float,default=0.5,help="Z resampling step for capacity integration [mm] (default: 0.5)")
    parser.add_argument("--output-dir",type=Path,help="Output directory (default: <source>_DrawingCapacity)")
    parser.add_argument("--self-test",action="store_true",help="Run a headless numerical self-test and exit")
    parser.add_argument("--version",action="version",version=f"%(prog)s {__version__}")
    args=parser.parse_args()
    if args.self_test:
        self_test(); return
    if args.z_step_mm <= 0 or args.render_dpi <= 0:
        parser.error("--z-step-mm and --render-dpi must be positive")
    try:
        import tkinter as tk
    except ModuleNotFoundError as exc:
        parser.exit(2,"ERROR: tkinter is required for GUI mode. Use a standard Python.org macOS/Windows Python build.\n")
    root=tk.Tk()
    DrawingCapacityGUI(root,args.input,args.page,args.render_dpi,args.z_step_mm,args.output_dir)
    root.mainloop()


if __name__ == "__main__":
    main()
