#!/usr/bin/env python3
"""
Septin ring size quantification pipeline (Python port of Macro_SeptinRing.ijm).

Measures upright septin rings (GFP) in CZI Z-stacks:
  - XY diameter (PRIMARY): object width = equivalent diameter of the
    segmented object on the MIP, sqrt(4*area/pi). The deconvolved MIPs
    show FILLED septin structures (bars/rods at the bud neck), not the
    hollow rings the macro's two-peak logic assumes, so the object width
    is the robust size measure.
  - XY circumference: pi * diameter
  - Legacy two-peak diameter (macro logic) kept as xy_diameter_two_peak_um
  - Z profile at the lobe-1 pixel: shown in the QC gallery only. The
    Z-EXTENT measurement was DROPPED (user decision): the single-pixel
    profile cuts through other structures/PSF tails, producing multi-peak
    profiles whose 50%-span massively overestimates the true ring height.

Preprocessing follows the nuclear/histone volume pipeline (the reference
pipeline for this project):
  background_subtract (median) -> theoretical PSF -> Richardson-Lucy
  deconvolution (40 iterations) -> Gaussian smoothing -> MIP.

Detection runs on the 2D MIP (like the macro) with the SAME adaptive
threshold mechanism as the volume pipeline: median + k * sigma_robust
(one-sided MAD below the median; see adaptive_seed_threshold).

NOTE: this script is STANDALONE -- no imports from other pipeline files.
All shared steps (CZI reading, background subtraction, PSF generation,
deconvolution, adaptive thresholding, Excel output) are implemented
directly below.

Deviation from the ImageJ macro (documented, deliberate):
  The macro measures the Z profile on the RAW 16-bit stack. Here the Z
  profile is taken from the DECONVOLVED (unsmoothed) stack, because
  deconvolution is part of the reference preprocessing. The profile is
  shown in the QC gallery only; no Z-extent is derived from it.

Outputs per run:
  all_rings.xlsx                 one row per detected ring + QC flags
  ring_diameter_per_strain.xlsx  one column per strain, ring diameters (um)
                                 in measurement order below (via --strain-map)
  qc_gallery_<image>.png/.svg    one panel per ring: XY crop with lobe line
                                 + Z profile (visual QC only)
  overview_mip_<image>.png/.svg  MIP with numbered ring ROIs
  crops_<image>/ring_*.png       one standalone panel PNG per ring (input
                                 for the interactive septin_ring_review.py)
  all_qc_metrics.xlsx            per-image QC metrics
  failed_images.tsv              batch failure log

Example:
  python septin_ring_pipeline.py --input ./czi_folder/ --output ./results \
      --channel-gfp 2 --strain-map strains.csv

  python septin_ring_pipeline.py --input W21417_1.czi --output ./out \
      --channel-gfp 2 --detection-threshold auto --detection-k 12
"""

import argparse
import os
import sys
import csv
import re
import time
import traceback
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Rectangle

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

from scipy.ndimage import gaussian_filter
from skimage.measure import regionprops, label
from skimage.restoration import richardson_lucy
from scipy import ndimage as ndi

rcParams['font.family'] = ['DejaVu Sans']
rcParams['svg.fonttype'] = 'none'

GFP_EMISSION_NM = 509


# ============================================================
# CZI single-channel reader
# ============================================================

def read_czi_single(filepath, channel_index):
    """Read one channel of a .czi file + metadata.

    channel_index: 0-based index into the squeezed channel axis.
    Returns (stack, (dz, dy, dx), na, immersion_ri, immersion).
    """
    import czifile
    with czifile.CziFile(filepath) as czi:
        data = czi.asarray().squeeze()
        metadata = czi.metadata()

    # Extract voxel sizes from metadata
    scaling_match = re.search(r'<Distance Id="X">\s*<Value>([\d.eE+-]+)</Value>', metadata)
    dx = float(scaling_match.group(1)) * 1e6 if scaling_match else 0.1024  # m -> um
    scaling_match = re.search(r'<Distance Id="Y">\s*<Value>([\d.eE+-]+)</Value>', metadata)
    dy = float(scaling_match.group(1)) * 1e6 if scaling_match else 0.1024
    scaling_match = re.search(r'<Distance Id="Z">\s*<Value>([\d.eE+-]+)</Value>', metadata)
    dz = float(scaling_match.group(1)) * 1e6 if scaling_match else 0.2

    # Extract NA and immersion
    na_match = re.search(r'<NumericalAperture>([\d.]+)</NumericalAperture>', metadata)
    # Find the 63x objective specifically
    obj_match = re.search(r'<Objective Name="([^"]*63[^"]*)"[^>]*>.*?<NumericalAperture>([\d.]+)</NumericalAperture>.*?<Immersions>([^<]+)</Immersions>', metadata, re.DOTALL)
    if obj_match:
        na = float(obj_match.group(2))
        immersion = obj_match.group(3)
    elif na_match:
        na = float(na_match.group(1))
        immersion = 'Oil'
    else:
        na = 1.4
        immersion = 'Oil'

    immersion_ri = {'Oil': 1.518, 'Water': 1.33, 'Air': 1.0, 'Glycerol': 1.47, 'Silicone': 1.40}.get(immersion, 1.518)

    channel = data[channel_index].astype(np.float32)
    return channel, (dz, dy, dx), na, immersion_ri, immersion


# ============================================================
# Background subtraction
# ============================================================

def background_subtract(stack):
    """Subtract stack-wide median background."""
    bg = np.median(stack)
    return np.clip(stack - bg, 0, None), bg


# ============================================================
# Theoretical PSF generation
# ============================================================

def theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength=607, psf_shape=(15, 31, 31),
                    modality='widefield'):
    """
    Generate a theoretical Gaussian PSF.
    emission_wavelength in nm (default 607 for mRFP1.2).

    modality='widefield': FWHM_lat = 0.51*lam/NA, FWHM_ax = 2.0*lam*n/NA^2
    modality='confocal' (spinning disk): FWHM_lat = 0.4*lam/NA,
        FWHM_ax = 1.4*lam*n/NA^2  (approx. confocal spot improvement)
    """
    lam = emission_wavelength / 1000.0  # um
    if modality == 'confocal':
        lateral_fwhm = 0.4 * lam / na
        axial_fwhm = 1.4 * lam * immersion_ri / (na ** 2)
    else:
        lateral_fwhm = 0.51 * lam / na
        axial_fwhm = 2.0 * lam * immersion_ri / (na ** 2)

    sigma_z = axial_fwhm / 2.355 / voxel_size[0]
    sigma_y = lateral_fwhm / 2.355 / voxel_size[1]
    sigma_x = lateral_fwhm / 2.355 / voxel_size[2]

    nz, ny, nx = psf_shape
    zz, yy, xx = np.indices((nz, ny, nx))
    cz, cy, cx = (nz - 1) / 2.0, (ny - 1) / 2.0, (nx - 1) / 2.0
    psf = np.exp(-(
        (zz - cz) ** 2 / (2 * sigma_z ** 2) +
        (yy - cy) ** 2 / (2 * sigma_y ** 2) +
        (xx - cx) ** 2 / (2 * sigma_x ** 2)
    )).astype(np.float32)
    psf /= psf.sum()
    return psf


# ============================================================
# Richardson-Lucy deconvolution
# ============================================================

def deconvolve(stack, psf, iterations=40):
    """Richardson-Lucy deconvolution via skimage."""
    stack_norm = stack / stack.max()
    deconv = richardson_lucy(stack_norm, psf, num_iter=iterations)
    return (deconv * stack.max()).astype(np.float32)


# ============================================================
# Adaptive seed threshold
# ============================================================

def adaptive_seed_threshold(smoothed, k=5.0, n_samples=5_000_000, rng_seed=0):
    """
    Noise-adaptive seed threshold: median + k * sigma_robust.

    sigma_robust is the ONE-SIDED MAD of the smoothed stack below its median
    (sigma = 1.4826 * median(median - x) over x < median). Bright objects sit
    above the median and therefore cannot inflate the noise estimate. Computed
    on a reproducible random subsample so it stays fast and memory-safe on
    large stacks. Because sigma is measured on the same (smoothed, possibly
    deconvolved) representation the threshold is applied to, k is
    scale-independent across microscopes, exposure times, and deconvolution.
    """
    flat = smoothed.ravel()
    if flat.size > n_samples:
        rng = np.random.default_rng(rng_seed)
        flat = flat[rng.integers(0, flat.size, size=n_samples)]
    med = float(np.median(flat))
    low = flat[flat < med]
    if low.size < 100:
        # Degenerate lower tail (e.g. heavy zero-clipping): fall back to the
        # upper 84.13th percentile as a 1-sigma proxy for the noise hump.
        return med + k * (float(np.percentile(flat, 84.134)) - med)
    sigma = 1.4826 * float(np.median(med - low))
    return med + k * sigma


# ============================================================
# Figure helpers
# ============================================================

def safe_tight_layout(fig):
    """Wrap tight_layout so a FreeType error doesn't crash everything."""
    try:
        fig.tight_layout()
    except Exception as e:
        print(f"    [warn] tight_layout failed ({type(e).__name__}), skipping.")


def safe_savefig(fig, path, format=None, dpi=150, bbox_inches='tight'):
    """
    Save a figure as SVG and/or PNG.

    SVG: saved via matplotlib (vector text, no FreeType rasterization).
    PNG: converted from SVG via cairosvg, bypassing matplotlib's Agg
    rasterizer entirely. This avoids the FT_Render_Glyph raster overflow
    that occurs on systems with FreeType < 2.13.2.

    Requires: pip install cairosvg
    """
    fmt = format or os.path.splitext(path)[1].lstrip('.').lower()

    if fmt == 'svg':
        try:
            fig.savefig(path, format='svg', bbox_inches=bbox_inches)
        except Exception:
            fig.savefig(path, format='svg', bbox_inches=None)
        return

    if fmt == 'png':
        # Save as SVG first (always works - vector, no rasterization)
        svg_temp = os.path.splitext(path)[0] + '_tmp.svg'
        try:
            fig.savefig(svg_temp, format='svg', bbox_inches=bbox_inches)
        except Exception:
            fig.savefig(svg_temp, format='svg', bbox_inches=None)

        # Convert SVG -> PNG via cairosvg (no matplotlib rasterizer)
        try:
            import cairosvg
            cairosvg.svg2png(url=svg_temp, write_to=path, scale=dpi/72.0)
        except ImportError:
            print(f"    [ERROR] cairosvg not installed. Install with: pip install cairosvg")
            print(f"            SVG saved instead: {os.path.basename(svg_temp)}")
            final_svg = os.path.splitext(path)[0] + '.svg'
            if not os.path.exists(final_svg):
                os.rename(svg_temp, final_svg)
            return
        except Exception as e:
            print(f"    [ERROR] PNG conversion failed: {e}")
            return
        finally:
            if os.path.exists(svg_temp):
                os.remove(svg_temp)
        return


# ============================================================
# Excel output helper
# ============================================================

# Column-specific number formats (Excel format strings).
# Any key not listed here falls back to FLOAT_FORMAT_DEFAULT for floats.
XLSX_NUMBER_FORMATS = {
    # Integer columns
    'pair_id':                 '0',
    'nuclear_object_id':       '0',
    'histone_object_id':       '0',
    'n_voxels_nuclear':        '0',
    'n_voxels_histone':        '0',
    'n_matched_pairs':         '0',
    # Volume / geometry (4 decimals)
    'nuclear_volume_um3':      '0.0000',
    'histone_volume_um3':      '0.0000',
    'ratio_percent':           '0.0000',
    'nuclear_mesh_volume_um3': '0.0000',
    'histone_mesh_volume_um3': '0.0000',
    'nuclear_surface_area_um2':'0.0000',
    'histone_surface_area_um2':'0.0000',
    'nuclear_sphericity':      '0.0000',
    'histone_sphericity':      '0.0000',
    'overlap_fraction':        '0.0000',
    'nuclear_z_span_um':       '0.0000',
    'histone_z_span_um':       '0.0000',
    'nuclear_xy_equiv_diameter_um': '0.0000',
    'histone_xy_equiv_diameter_um': '0.0000',
    'nuclear_z_xy_ratio':      '0.0000',
    'histone_z_xy_ratio':      '0.0000',
    'axial_correction_factor': '0.0000',
    # Intensities (2 decimals)
    'nuclear_mean_intensity':  '0.00',
    'histone_mean_intensity':  '0.00',
    'nuclear_max_intensity':   '0.00',
    'histone_max_intensity':   '0.00',
    'nuclear_median_intensity':'0.00',
    'histone_median_intensity':'0.00',
    'nuclear_total_intensity': '0.00',
    'histone_total_intensity': '0.00',
    # Summary statistics (4 decimals)
    'ratio_median_percent':    '0.0000',
    'ratio_mean_percent':      '0.0000',
    'ratio_std_percent':       '0.0000',
    'nuclear_volume_median_um3':'0.0000',
    'histone_volume_median_um3':'0.0000',
    'intensity_vol_pearson_r_nuclear': '0.0000',
    'intensity_vol_pearson_r_histone': '0.0000',
}
FLOAT_FORMAT_DEFAULT = '0.0000'

HEADER_FONT = Font(bold=True)
HEADER_FILL = PatternFill(start_color='D9D9D9', end_color='D9D9D9', fill_type='solid')
HEADER_ALIGN = Alignment(horizontal='center', vertical='center')


def save_xlsx(records, path, sheet_name):
    """Save a list of dict records as a formatted .xlsx file.

    - Bold header row with light-grey background, centered
    - Frozen top row
    - Column widths auto-fit to content (capped at 50 chars)
    - Numeric columns get context-appropriate number formats

    Empty records: writes an empty file with no header, so downstream tools
    still see a file exists (matches previous CSV behavior of writing nothing).
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet_name

    if not records:
        try:
            wb.save(path)
        except PermissionError:
            raise PermissionError(
                f"Cannot write '{path}'. File is likely open in Excel. "
                f"Please close it and re-run the pipeline."
            )
        return

    keys = list(records[0].keys())

    # Header row
    for col_idx, key in enumerate(keys, start=1):
        cell = ws.cell(row=1, column=col_idx, value=key)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = HEADER_ALIGN

    # Data rows
    for row_idx, record in enumerate(records, start=2):
        for col_idx, key in enumerate(keys, start=1):
            value = record.get(key)

            # Booleans -> "True"/"False" as text (openpyxl would otherwise
            # write TRUE/FALSE which reads odd in a scientific table)
            if isinstance(value, (bool, np.bool_)):
                cell = ws.cell(row=row_idx, column=col_idx, value=str(bool(value)))
            elif isinstance(value, (np.integer,)):
                cell = ws.cell(row=row_idx, column=col_idx, value=int(value))
            elif isinstance(value, (np.floating,)):
                cell = ws.cell(row=row_idx, column=col_idx, value=float(value))
            else:
                cell = ws.cell(row=row_idx, column=col_idx, value=value)

            # Apply number format
            if isinstance(value, (bool, np.bool_)):
                pass  # string, no format
            elif isinstance(value, (int, float, np.integer, np.floating)):
                fmt = XLSX_NUMBER_FORMATS.get(key)
                if fmt is None:
                    fmt = '0' if isinstance(value, (int, np.integer)) else FLOAT_FORMAT_DEFAULT
                cell.number_format = fmt

    # Freeze header row
    ws.freeze_panes = 'A2'

    # Auto-size columns based on content length
    for col_idx, key in enumerate(keys, start=1):
        max_len = len(str(key))
        for record in records:
            value = record.get(key)
            if value is None:
                continue
            # Approximate displayed length using the formatted number
            if isinstance(value, (bool, np.bool_)):
                s = str(bool(value))
            elif isinstance(value, float) or isinstance(value, np.floating):
                fmt = XLSX_NUMBER_FORMATS.get(key, FLOAT_FORMAT_DEFAULT)
                decimals = len(fmt.split('.')[-1]) if '.' in fmt else 0
                s = f"{float(value):.{decimals}f}"
            else:
                s = str(value)
            if len(s) > max_len:
                max_len = len(s)
        # Cap at 50 chars, minimum 8 for readability, +2 padding
        width = min(max(max_len + 2, 8), 50)
        ws.column_dimensions[get_column_letter(col_idx)].width = width

    try:
        wb.save(path)
    except PermissionError:
        raise PermissionError(
            f"Cannot write '{path}'. File is likely open in Excel. "
            f"Please close it and re-run the pipeline."
        )


# ============================================================
# Detection on the MIP (macro-equivalent)
# ============================================================

def detect_rings_mip(mip, voxel_yx, detection_threshold='auto', detection_k=5.0,
                     min_ring_area_um2=None, max_ring_area_um2=None):
    """Detect ring candidates on the 2D MIP.

    Macro equivalent: threshold -> convert to mask -> analyze particles
    (size filter). Deliberate deviation: NO watershed split (the ImageJ
    watershed splits hollow rings into their two lobes; connected
    components keep each ring whole, and merged neighbors are caught by
    the 'suspicious_aspect' QC flag in the gallery review).

    detection_threshold: float (absolute intensity on the MIP) or 'auto'
        (median + k * sigma_robust via adaptive_seed_threshold, the same
        mechanism as the volume pipeline).
    min/max_ring_area_um2: size filter in physical units. Defaults match
        the macro's 2-5000 px at the W21616 pixel size (0.10238 um/px):
        2 px = 0.021 um2, 5000 px = 52.4 um2.

    Returns (labels, n_objects, effective_threshold, qc dict).
    """
    dy, dx = voxel_yx
    px_area = dy * dx
    if min_ring_area_um2 is None:
        # Calibrated on W21417/W21418 septin-GFP: real septin structures have
        # median area ~200 px; the macro's 2 px floor admits only noise
        # speckles. 20 px at 0.10238 um/px = 0.21 um2 removes the speckles
        # while keeping all real structures.
        min_ring_area_um2 = 20 * 0.10238 ** 2     # 20 px at 0.10238 um/px
    if max_ring_area_um2 is None:
        max_ring_area_um2 = 5000 * 0.10238 ** 2   # macro: 5000 px
    min_px = max(1, int(round(min_ring_area_um2 / px_area)))
    max_px = max(min_px, int(round(max_ring_area_um2 / px_area)))

    if isinstance(detection_threshold, str):
        if detection_threshold.lower() != 'auto':
            raise ValueError(
                f"detection_threshold must be a number or 'auto', got "
                f"'{detection_threshold}'")
        thr = adaptive_seed_threshold(mip, k=detection_k)
    else:
        thr = float(detection_threshold)

    mask = mip > thr
    # NOTE on watershed: the ImageJ macro runs Watershed on the binary mask,
    # which splits a hollow ring into its two lobes whenever the lobe
    # separation exceeds the lobe width (each lobe is an ultimate eroded
    # point). The macro's two-peak diameter measurement then still works
    # because its 41x41 px crop is centered on the ROI and spans both halves
    # -- but every ring appears TWICE in the review. Here we keep connected
    # components as ring candidates (no watershed split), so each ring is
    # detected and measured once. Touching rings that share one component
    # are flagged 'suspicious_aspect' for the manual gallery review instead
    # of being force-split.
    labels = label(mask) if mask.any() else np.zeros_like(mask, dtype=np.int32)

    # Size filter (macro: Analyze Particles size=min-max)
    out = np.zeros_like(labels)
    n = 0
    for prop in regionprops(labels):
        if min_px <= prop.area <= max_px:
            n += 1
            out[labels == prop.label] = n

    qc = {
        '_detection_threshold_effective': thr,
        'detection_threshold_mode': ('auto_k%.1f' % detection_k)
        if isinstance(detection_threshold, str) else 'fixed',
        'min_ring_area_px': min_px,
        'max_ring_area_px': max_px,
        'n_candidates_raw': int(labels.max()),
        'n_candidates_size_filtered': n,
    }
    return out, n, thr, qc


# ============================================================
# Per-ring measurements (macro-equivalent)
# ============================================================

def measure_ring_xy(mip, labels, ring_id, voxel_yx, crop_half_width=20,
                    min_peak_sep_px=3, detection_thr=0.0):
    """XY diameter via the macro's two-brightest-lobes logic.

    Crop (2*crop_half_width+1)^2 around the ROI centroid on the smoothed
    MIP; find the global maximum (lobe 1) and the second maximum at least
    min_peak_sep_px away (lobe 2). Diameter = distance * pixel size.

    Returns dict with positions, diameter, circumference, validity flag.
    """
    dy, dx = voxel_yx
    props = [p for p in regionprops(labels) if p.label == ring_id]
    if not props:
        raise ValueError(f"ring_id {ring_id} not in labels")
    prop = props[0]
    cy, cx = prop.centroid  # row, col (float)
    cy_i, cx_i = int(round(cy)), int(round(cx))

    ny, nx = mip.shape
    y1 = max(0, cy_i - crop_half_width)
    y2 = min(ny - 1, cy_i + crop_half_width)
    x1 = max(0, cx_i - crop_half_width)
    x2 = min(nx - 1, cx_i + crop_half_width)
    crop = mip[y1:y2 + 1, x1:x2 + 1]

    # Pass 1: global maximum (macro scans the 8-bit crop; here the
    # smoothed MIP crop -- same argmax).
    flat_idx = int(np.argmax(crop))
    lobe1_y, lobe1_x = np.unravel_index(flat_idx, crop.shape)
    max1_val = float(crop[lobe1_y, lobe1_x])

    # Pass 2: second maximum at least min_peak_sep_px from lobe 1.
    yy, xx = np.indices(crop.shape)
    dist = np.sqrt((yy - lobe1_y) ** 2 + (xx - lobe1_x) ** 2)
    eligible = dist >= min_peak_sep_px
    if eligible.any():
        masked = np.where(eligible, crop, -np.inf)
        flat2 = int(np.argmax(masked))
        lobe2_y, lobe2_x = np.unravel_index(flat2, crop.shape)
        max2_val = float(crop[lobe2_y, lobe2_x])
    else:
        lobe2_y, lobe2_x, max2_val = lobe1_y, lobe1_x, -np.inf

    # Macro validity: both peaks above the detection threshold.
    xy_valid = (max1_val > detection_thr) and (max2_val > detection_thr)
    if xy_valid:
        dist_px = float(np.sqrt((lobe1_x - lobe2_x) ** 2 +
                                (lobe1_y - lobe2_y) ** 2))
        xy_diam_um = dist_px * dy
        xy_circ_um = np.pi * xy_diam_um
    else:
        dist_px = np.nan
        xy_diam_um = np.nan
        xy_circ_um = np.nan

    # Bright-pixel count in the crop (macro quality indicator).
    n_bright = int((crop > detection_thr).sum())

    # ---- Object-width diameter (user-chosen primary measure) ----
    # The deconvolved MIPs show FILLED septin structures (bars/short rods at
    # the bud neck), not the hollow rings the macro's two-peak logic assumes.
    # The primary diameter is therefore the segmented object's width:
    #   - equivalent diameter (sqrt(4*area/pi)) for compact objects
    #   - major/minor axis length of the fitted ellipse for elongated ones
    # All in physical units (um).
    area_um2 = prop.area * dy * dx
    equiv_diam_um = float(np.sqrt(4.0 * area_um2 / np.pi))
    try:
        # skimage >= 0.26 names (old ones deprecated, removed in 2.0)
        major_axis_um = float(prop.axis_major_length * dy)
        minor_axis_um = float(prop.axis_minor_length * dy)
    except AttributeError:
        major_axis_um = float(prop.major_axis_length * dy)
        minor_axis_um = float(prop.minor_axis_length * dy)

    return {
        'crop_origin_yx': (y1, x1),
        'crop': crop,
        'lobe1_crop_yx': (int(lobe1_y), int(lobe1_x)),
        'lobe2_crop_yx': (int(lobe2_y), int(lobe2_x)),
        'lobe1_global_yx': (y1 + int(lobe1_y), x1 + int(lobe1_x)),
        'lobe2_global_yx': (y1 + int(lobe2_y), x1 + int(lobe2_x)),
        'lobe1_peak': max1_val,
        'lobe2_peak': max2_val,
        'xy_valid': xy_valid,
        'xy_diameter_um': xy_diam_um,          # macro two-peak (legacy)
        'xy_circumference_um': xy_circ_um,     # macro two-peak (legacy)
        'lobe_distance_px': dist_px,
        'object_equiv_diameter_um': equiv_diam_um,   # PRIMARY (user choice)
        'object_major_axis_um': major_axis_um,
        'object_minor_axis_um': minor_axis_um,
        'object_area_um2': area_um2,
        'n_bright_pixels': n_bright,
        'centroid_y_um': cy * dy,
        'centroid_x_um': cx * dx,
        'bbox': prop.bbox,  # (min_row, min_col, max_row, max_col)
        'area_px': prop.area,
    }


def measure_ring_z(deconv_stack, lobe1_global_yx, dz):
    """Z profile at the lobe-1 pixel -- kept for the QC GALLERY only.

    The Z-EXTENT (ring height) measurement was DROPPED (user decision):
    the single-pixel Z profile cuts through other structures and axial PSF
    tails, producing multi-peak profiles whose first-to-last-slice span
    massively overestimates the true ring height (~40% of rings affected,
    outliers to 8+ um for ~0.5 um objects). The profile is still shown in
    the gallery as a visual QC aid, but no height is derived from it.
    """
    gy, gx = lobe1_global_yx
    z_profile = deconv_stack[:, gy, gx].astype(np.float64)
    peak = float(z_profile.max())
    return {
        'z_profile': z_profile,
        'z_peak': peak,
    }


# ============================================================
# QC flags
# ============================================================

def qc_flags_for_ring(xy, z, mip_shape, n_z_slices, aspect_ratio_max=3.0,
                      min_bright_px=3):
    """Per-ring QC flags (review aid for the gallery click-through)."""
    flags = []
    if not xy['xy_valid']:
        flags.append('xy_not_measurable')
    if xy['n_bright_pixels'] < min_bright_px:
        flags.append('few_bright_pixels')
    min_r, min_c, max_r, max_c = xy['bbox']
    ny, nx = mip_shape
    if min_r == 0 or min_c == 0 or max_r == ny or max_c == nx:
        flags.append('touches_border')
    h = max_r - min_r
    w = max_c - min_c
    if min(h, w) > 0 and max(h, w) / min(h, w) > aspect_ratio_max:
        flags.append('suspicious_aspect')
    return ';'.join(flags) if flags else 'none'


# ============================================================
# QC gallery figure
# ============================================================

def generate_qc_gallery(image_id, mip, ring_records, output_dir,
                        max_panels_per_fig=24):
    """One panel per ring: XY crop with lobe line + Z profile.

    Layout mirrors the macro's review window: left = annotated XY crop
    (cyan lobe line, yellow peak crosses), right = Z profile (green),
    threshold (red), extent boundaries (blue).
    """
    if not ring_records:
        return
    n = len(ring_records)
    n_figs = int(np.ceil(n / max_panels_per_fig))
    for fi in range(n_figs):
        chunk = ring_records[fi * max_panels_per_fig:(fi + 1) * max_panels_per_fig]
        nrows = len(chunk)
        fig, axes = plt.subplots(nrows, 2, figsize=(8, 2.2 * nrows),
                                 squeeze=False)
        for ri, rec in enumerate(chunk):
            ax_xy, ax_z = axes[ri, 0], axes[ri, 1]
            crop = rec['xy']['crop']
            # contrast stretch like the macro's Enhance Contrast
            lo, hi = np.percentile(crop, [0.35, 99.65])
            ax_xy.imshow(crop, cmap='gray', vmin=lo, vmax=hi,
                         interpolation='nearest')
            if rec['xy']['xy_valid']:
                (y1, x1), (y2, x2) = (rec['xy']['lobe1_crop_yx'],
                                      rec['xy']['lobe2_crop_yx'])
                ax_xy.plot([x1, x2], [y1, y2], color='cyan', lw=1.2)
                for (py, px) in [(y1, x1), (y2, x2)]:
                    ax_xy.plot(px, py, marker='+', color='yellow',
                               markersize=8, markeredgewidth=1.5)
            diam = rec['xy']['xy_diameter_um']
            title = f"Ring {rec['ring_id']}"
            if np.isfinite(diam):
                title += f"  d={diam:.3f} um"
            else:
                title += "  d=n/a"
            ax_xy.set_title(title, fontsize=8)
            ax_xy.set_xticks([]); ax_xy.set_yticks([])

            zp = rec['z']['z_profile']
            zs = np.arange(1, len(zp) + 1)
            ax_z.plot(zs, zp, color='green', lw=1.2)
            ztitle = 'Z profile (QC only)'
            if rec['qc_flags'] != 'none':
                ztitle += f"  [{rec['qc_flags']}]"
            ax_z.set_title(ztitle, fontsize=8)
            ax_z.set_xlabel('Z slice', fontsize=7)
            ax_z.set_ylabel('Intensity', fontsize=7)
            ax_z.tick_params(labelsize=7)
        suffix = f'_p{fi + 1}' if n_figs > 1 else ''
        base = os.path.join(output_dir, f'qc_gallery_{image_id}{suffix}')
        safe_tight_layout(fig)
        safe_savefig(fig, base + '.svg')
        safe_savefig(fig, base + '.png')
        plt.close(fig)


def save_ring_crops(image_id, ring_records, output_dir):
    """Save one standalone panel PNG per ring for the interactive review
    tool (septin_ring_review.py).

    Each panel is identical to the QC-gallery panel: XY crop with lobe
    line + peak crosses (left) and the Z profile (right, visual QC only).
    Files go to crops_<image_id>/ring_<id>.png.
    """
    if not ring_records:
        return
    crops_dir = os.path.join(output_dir, f'crops_{image_id}')
    os.makedirs(crops_dir, exist_ok=True)
    for rec in ring_records:
        fig, (ax_xy, ax_z) = plt.subplots(1, 2, figsize=(8, 2.2))
        crop = rec['xy']['crop']
        lo, hi = np.percentile(crop, [0.35, 99.65])
        ax_xy.imshow(crop, cmap='gray', vmin=lo, vmax=hi,
                     interpolation='nearest')
        if rec['xy']['xy_valid']:
            (y1, x1), (y2, x2) = (rec['xy']['lobe1_crop_yx'],
                                  rec['xy']['lobe2_crop_yx'])
            ax_xy.plot([x1, x2], [y1, y2], color='cyan', lw=1.2)
            for (py, px) in [(y1, x1), (y2, x2)]:
                ax_xy.plot(px, py, marker='+', color='yellow',
                           markersize=8, markeredgewidth=1.5)
        diam = rec['xy']['xy_diameter_um']
        title = f"Ring {rec['ring_id']}"
        if np.isfinite(diam):
            title += f"  d={diam:.3f} um"
        else:
            title += "  d=n/a"
        ax_xy.set_title(title, fontsize=8)
        ax_xy.set_xticks([]); ax_xy.set_yticks([])

        zp = rec['z']['z_profile']
        zs = np.arange(1, len(zp) + 1)
        ax_z.plot(zs, zp, color='green', lw=1.2)
        ztitle = 'Z profile (QC only)'
        if rec['qc_flags'] != 'none':
            ztitle += f"  [{rec['qc_flags']}]"
        ax_z.set_title(ztitle, fontsize=8)
        ax_z.set_xlabel('Z slice', fontsize=7)
        ax_z.set_ylabel('Intensity', fontsize=7)
        ax_z.tick_params(labelsize=7)

        safe_tight_layout(fig)
        safe_savefig(fig, os.path.join(
            crops_dir, f"ring_{rec['ring_id']:03d}.png"))
        plt.close(fig)


def generate_overview_mip(image_id, mip, labels, ring_records, output_dir):
    """MIP with numbered ring ROIs (orientation map for the gallery)."""
    fig, ax = plt.subplots(figsize=(10, 10))
    lo, hi = np.percentile(mip, [0.35, 99.65])
    ax.imshow(mip, cmap='gray', vmin=lo, vmax=hi, interpolation='nearest')
    for rec in ring_records:
        min_r, min_c, max_r, max_c = rec['xy']['bbox']
        ax.add_patch(Rectangle((min_c, min_r), max_c - min_c, max_r - min_r,
                               fill=False, edgecolor='cyan', lw=1.0))
        # keep the label inside the canvas (clip-on) and clear of the title
        label_r = max(min_r - 2, 8)
        ax.text(min_c, label_r, str(rec['ring_id']), color='yellow',
                fontsize=8, clip_on=True)
    ax.set_title(f'{image_id}: {len(ring_records)} rings')
    ax.set_xticks([]); ax.set_yticks([])
    base = os.path.join(output_dir, f'overview_mip_{image_id}')
    safe_tight_layout(fig)
    safe_savefig(fig, base + '.svg')
    safe_savefig(fig, base + '.png')
    plt.close(fig)


# ============================================================
# Per-image processing
# ============================================================

def process_image(filepath, output_dir, channel_index,
                  sample_ri=1.35, emission_wavelength=GFP_EMISSION_NM,
                  deconv_iterations=40, strain=None, modality=None,
                  detection_threshold='auto', detection_k=5.0,
                  min_ring_area_um2=None, max_ring_area_um2=None,
                  crop_half_width=20, min_peak_sep_px=3,
                  dz_um=None, max_z_slices=None,
                  no_deconvolution=False):
    """Full per-image flow: read -> preprocess -> MIP -> detect -> measure."""
    t0 = time.time()
    os.makedirs(output_dir, exist_ok=True)
    image_id = os.path.splitext(os.path.basename(filepath))[0]
    print(f"\n{'=' * 60}\nProcessing: {image_id}\n{'=' * 60}")

    stack, voxel, na, immersion_ri, immersion = read_czi_single(
        filepath, channel_index)
    if max_z_slices is not None and stack.shape[0] > max_z_slices:
        print(f"  [trim] Z-stack {stack.shape[0]} -> {max_z_slices} slices")
        stack = stack[:max_z_slices]
    dz, dy, dx = voxel
    if dz_um is not None:
        print(f"  [override] dz = {dz_um} um (metadata PhysicalSizeZ ignored)")
        dz = dz_um
        voxel = (dz, dy, dx)
    print(f"  Channel shape: {stack.shape}, voxel: ({dz}, {dy}, {dx}), "
          f"NA={na}, immersion={immersion}")

    # --- preprocessing (reference volume-pipeline steps) ---
    stack_bg, bg = background_subtract(stack)
    print(f"  Background (median): {bg:.1f}")

    modality_eff = modality or 'widefield'  # CZI = Zeiss widefield
    if no_deconvolution:
        deconv = stack_bg
        print("  Deconvolution: SKIPPED (--no-deconvolution)")
    else:
        psf = theoretical_psf(voxel, na, immersion_ri,
                              emission_wavelength=emission_wavelength,
                              modality=modality_eff)
        print(f"  PSF: {psf.shape} ({modality_eff}, {emission_wavelength} nm)")
        td = time.time()
        deconv = deconvolve(stack_bg, psf, iterations=deconv_iterations)
        print(f"  Deconvolution: {time.time() - td:.1f}s "
              f"({deconv_iterations} iterations)")

    # Gaussian smoothing (same physical sigma as the volume pipeline)
    smooth_sigma_um = (0.2, 0.20476190476190476, 0.20476190476190476)
    smooth_sigma = tuple(s / v for s, v in zip(smooth_sigma_um, voxel))
    smoothed = gaussian_filter(deconv.astype(np.float32), sigma=smooth_sigma)

    # --- MIP + detection ---
    mip = smoothed.max(axis=0)
    labels, n_rings, thr, det_qc = detect_rings_mip(
        mip, (dy, dx), detection_threshold=detection_threshold,
        detection_k=detection_k, min_ring_area_um2=min_ring_area_um2,
        max_ring_area_um2=max_ring_area_um2)
    if isinstance(detection_threshold, str):
        print(f"  Adaptive detection ('auto', k={detection_k}): "
              f"threshold={thr:.2f}")
    else:
        print(f"  Fixed detection threshold: {thr:.2f}")
    print(f"  Rings detected: {n_rings} "
          f"(raw {det_qc['n_candidates_raw']}, size-filtered)")

    # --- per-ring measurement ---
    records = []
    gallery_records = []
    for rid in range(1, n_rings + 1):
        xy = measure_ring_xy(mip, labels, rid, (dy, dx),
                             crop_half_width=crop_half_width,
                             min_peak_sep_px=min_peak_sep_px,
                             detection_thr=thr)
        z = measure_ring_z(deconv, xy['lobe1_global_yx'], dz)
        flags = qc_flags_for_ring(xy, z, mip.shape, deconv.shape[0])
        rec = {
            'strain': strain if strain else 'unassigned',
            'image_id': image_id,
            'ring_id': rid,
            'x_center_um': xy['centroid_x_um'],
            'y_center_um': xy['centroid_y_um'],
            'xy_diameter_um': xy['object_equiv_diameter_um'],  # PRIMARY: object width
            'xy_circumference_um': (np.pi * xy['object_equiv_diameter_um']
                                    if np.isfinite(xy['object_equiv_diameter_um'])
                                    else np.nan),
            'xy_diameter_two_peak_um': xy['xy_diameter_um'],   # legacy macro value
            'object_major_axis_um': xy['object_major_axis_um'],
            'object_minor_axis_um': xy['object_minor_axis_um'],
            'object_area_um2': xy['object_area_um2'],
            'n_bright_pixels': xy['n_bright_pixels'],
            'lobe1_peak_intensity': xy['lobe1_peak'],
            'lobe2_peak_intensity': xy['lobe2_peak'],
            'z_peak_intensity': z['z_peak'],
            'detection_threshold_counts': thr,
            'area_px': xy['area_px'],
            'qc_flags': flags,
        }
        records.append(rec)
        gallery_records.append({'ring_id': rid, 'xy': xy, 'z': z,
                                'qc_flags': flags})

    # --- figures ---
    if n_rings:
        generate_overview_mip(image_id, mip, labels, gallery_records,
                              output_dir)
        generate_qc_gallery(image_id, mip, gallery_records, output_dir)
        save_ring_crops(image_id, gallery_records, output_dir)

    qc = {
        'image_id': image_id,
        'strain': strain if strain else 'unassigned',
        'n_rings': n_rings,
        'background': bg,
        'modality': modality_eff,
        'deconvolved': not no_deconvolution,
        **det_qc,
    }
    print(f"  Done in {time.time() - t0:.1f}s — {n_rings} rings")
    return records, qc


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Septin ring size quantification (Python port of '
                    'Macro_SeptinRing.ijm) with volume-pipeline '
                    'preprocessing (deconvolution).')
    parser.add_argument('--input', required=True,
                        help='Input .czi file or folder of .czi files')
    parser.add_argument('--output', required=True,
                        help='Output directory for results')
    parser.add_argument('--channel-gfp', type=int, default=0,
                        help='0-based index of the GFP channel in the .czi '
                             '(default: 0 = first channel). For the old '
                             'dual-channel CZIs GFP is typically index 1 '
                             '(second channel).')
    parser.add_argument('--sample-ri', type=float, default=1.35,
                        help='Sample refractive index (default: 1.35 yeast)')
    parser.add_argument('--deconv-iterations', type=int, default=40,
                        help='RL deconvolution iterations (default: 40)')
    parser.add_argument('--no-deconvolution', action='store_true',
                        help='Skip deconvolution (faster, less accurate)')
    parser.add_argument('--modality', choices=['widefield', 'confocal'],
                        default=None,
                        help='PSF model for deconvolution (default: '
                             'widefield for .czi).')
    parser.add_argument('--emission-nm', type=float, default=GFP_EMISSION_NM,
                        help='Emission wavelength of the septin marker '
                             '(default: 509 nm GFP).')
    parser.add_argument('--detection-threshold', default='auto',
                        metavar='FLOAT|auto',
                        help="Ring detection threshold on the smoothed MIP: "
                             "'auto' (default; median + k*sigma_robust, the "
                             "same adaptive mechanism as the volume "
                             "pipeline) or an absolute intensity value.")
    parser.add_argument('--detection-k', type=float, default=12.0,
                        help="Factor k for --detection-threshold auto "
                             "(default: 12.0, calibrated on W21417/W21418 "
                             "septin-GFP data: isolates the real septin "
                             "structures, k=5 is dominated by noise speckles). "
                             "Ignored otherwise.")
    parser.add_argument('--min-ring-area-um2', type=float, default=None,
                        help='Minimum ring area on the MIP in um^2 '
                             '(default: 20 px at 0.10238 um/px = 0.21 um^2, '
                             'calibrated on W21417/W21418 to remove noise '
                             'speckles; the macro default 2 px admits only '
                             'noise on this data).')
    parser.add_argument('--max-ring-area-um2', type=float, default=None,
                        help='Maximum ring area on the MIP in um^2 '
                             '(default: macro 5000 px at 0.10238 um/px = '
                             '52.4 um^2).')
    parser.add_argument('--crop-half-width-px', type=int, default=20,
                        help='Half-width of the XY crop around each ring '
                             '(default: 20 px, macro setting).')
    parser.add_argument('--min-peak-sep-px', type=int, default=3,
                        help='Minimum separation of the two lobe peaks '
                             '(default: 3 px, macro setting).')
    parser.add_argument('--dz-um', type=float, default=None,
                        help='Override the z voxel size in um when the '
                             'metadata value is wrong.')
    parser.add_argument('--max-z-slices', type=int, default=None,
                        help='Trim the Z-stack to the first N slices '
                             '(e.g. 53 to match the 3n/4n acquisition).')
    parser.add_argument('--strain-map',
                        help='CSV mapping filenames to strain labels '
                             '(columns: filename, strain).')
    args = parser.parse_args()

    # --detection-threshold accepts 'auto' or a number
    if isinstance(args.detection_threshold, str):
        if args.detection_threshold.lower() == 'auto':
            args.detection_threshold = 'auto'
        else:
            try:
                args.detection_threshold = float(args.detection_threshold)
            except ValueError:
                raise SystemExit(
                    "Error: --detection-threshold must be a number or "
                    f"'auto', got '{args.detection_threshold}'.")

    os.makedirs(args.output, exist_ok=True)

    # Build file list
    if os.path.isdir(args.input):
        files = sorted([os.path.join(args.input, f)
                        for f in os.listdir(args.input)
                        if f.lower().endswith('.czi')])
    else:
        files = [args.input]
    if not files:
        print(f"No .czi files found in {args.input}")
        sys.exit(1)
    print(f"Found {len(files)} .czi file(s) to process.")

    # Strain map
    strain_map = {}
    if args.strain_map:
        with open(args.strain_map) as f:
            for row in csv.DictReader(f):
                strain_map[row['filename'].strip()] = row['strain'].strip()
        print(f"Loaded strain map: {len(strain_map)} entries from "
              f"{args.strain_map}")

    def lookup_strain(*keys):
        for k in keys:
            if k is None:
                continue
            k = k.strip()
            if k in strain_map:
                return strain_map[k]
        return None

    failed_path = f'{args.output}/failed_images.tsv'
    with open(failed_path, 'w') as f:
        f.write('filename\terror_type\terror_message\ttimestamp\n')

    all_records = []
    all_qc = []
    n_success = n_failed = 0
    for filepath in files:
        fname = os.path.basename(filepath)
        strain = lookup_strain(fname, os.path.splitext(fname)[0])
        if strain is None:
            print(f"  [WARN] No strain assignment for '{fname}' -> "
                  f"'unassigned'.")
        try:
            records, qc = process_image(
                filepath, args.output,
                channel_index=args.channel_gfp,
                sample_ri=args.sample_ri,
                emission_wavelength=args.emission_nm,
                deconv_iterations=args.deconv_iterations,
                strain=strain, modality=args.modality,
                detection_threshold=args.detection_threshold,
                detection_k=args.detection_k,
                min_ring_area_um2=args.min_ring_area_um2,
                max_ring_area_um2=args.max_ring_area_um2,
                crop_half_width=args.crop_half_width_px,
                min_peak_sep_px=args.min_peak_sep_px,
                dz_um=args.dz_um,
                max_z_slices=args.max_z_slices,
                no_deconvolution=args.no_deconvolution,
            )
            all_records.extend(records)
            all_qc.append(qc)
            n_success += 1
        except KeyboardInterrupt:
            print(f"\n  [ABORT] KeyboardInterrupt during {fname}. "
                  f"Stopping batch.")
            raise
        except Exception as e:
            n_failed += 1
            tb = traceback.format_exc()
            err_type = type(e).__name__
            err_msg = str(e).replace('\n', ' ').replace('\t', ' ')
            ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(f"\n  [FAILED] {fname}: {err_type}: {err_msg}")
            print(tb)
            with open(failed_path, 'a') as f:
                f.write(f'{fname}\t{err_type}\t{err_msg}\t{ts}\n')
            continue

    # Combined results
    if all_records:
        combined = f'{args.output}/all_rings.xlsx'
        save_xlsx(all_records, combined, sheet_name='all_rings')
        print(f"\nCombined results: {combined} "
              f"({len(all_records)} rings total)")

        # Ring diameters per strain: one column per strain, diameters (um)
        # in measurement order below. Shorter columns stay empty at the
        # bottom. Only finite (measurable) diameters are included.
        strains = sorted({r['strain'] for r in all_records})
        diam_by_strain = {
            s: [r['xy_diameter_um'] for r in all_records
                if r['strain'] == s and np.isfinite(r['xy_diameter_um'])]
            for s in strains
        }
        n_max = max((len(v) for v in diam_by_strain.values()), default=0)
        diam_rows = [{s: (diam_by_strain[s][i] if i < len(diam_by_strain[s])
                          else None)
                      for s in strains}
                     for i in range(n_max)]
        diam_path = f'{args.output}/ring_diameter_per_strain.xlsx'
        save_xlsx(diam_rows, diam_path, sheet_name='ring_diameter_per_strain')
        print(f"Ring diameters per strain: {diam_path}")

    if all_qc:
        save_xlsx(all_qc, f'{args.output}/all_qc_metrics.xlsx',
                  sheet_name='qc')

    print(f"\nPipeline complete. {n_success} of {len(files)} images "
          f"successful, {n_failed} failed, {len(all_records)} rings total.")


if __name__ == '__main__':
    main()
