"""
Yeast Nuclear/Histone Volume Ratio Quantification Pipeline
==========================================================
Processes widefield fluorescence z-stacks (.czi) to measure yeast nuclear
and histone volumes separately, match corresponding objects, and compute
the histone/nuclear volume ratio per cell.

Pipeline:
  1. Read .czi and extract both channels + metadata
  2. Background subtraction (median) per channel
  3. Theoretical PSF generation per channel (different emission wavelengths)
  4. Richardson-Lucy deconvolution per channel (skimage, 40 iterations)
  5. Axial RI correction (oil -> aqueous sample)
  6. 3D Gaussian smoothing per channel
  7. Seed detection + size/z-span filtering per channel
  8. Watershed separation of touching nuclei per channel
  9. Per-object adaptive Otsu thresholding per channel
  10. Object matching between channels (histone centroid containment in nucleus)
  11. Volume + shape measurement per matched pair
  12. QC: intensity-volume correlation, threshold sensitivity, ratio distribution

Usage:
  python nuclear_histone_ratio_pipeline_1n_2n.py --input <folder> --output <results> [options]

  # Process a single file
  python nuclear_histone_ratio_pipeline_1n_2n.py --input image.czi --output ./results --channel-x 1 --channel-y 2

  # Process all .czi files in a folder
  python nuclear_histone_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x 1 --channel-y 2

  # With strain mapping CSV (columns: filename, strain)
  python nuclear_histone_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x 1 --channel-y 2 --strain-map strains.csv

Channel selection (1-based, as in ZEN):
  --channel-x = nuclear marker channel (1 = mRFP1.2, 2 = EGFP)
  --channel-y = histone marker channel (1 = mRFP1.2, 2 = EGFP)
  The emission wavelength follows automatically from the channel choice
  (607 nm for mRFP1.2, 509 nm for EGFP) -- no extra flags needed.

  Example: nuclear marker = mRFP1.2 (channel 1), histone = EGFP (channel 2):
    python nuclear_histone_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x 1 --channel-y 2

  Example with swapped fluorophores (nuclear = EGFP, histone = mRFP1.2):
    python nuclear_histone_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x 2 --channel-y 1

Requirements:
  pip install czifile scikit-image scipy numpy matplotlib openpyxl cairosvg
"""

import argparse
import os
import sys
import csv
import time
import traceback
from datetime import datetime
import numpy as np
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr, spearmanr
from skimage.restoration import richardson_lucy
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.measure import marching_cubes, mesh_surface_area

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = ['DejaVu Sans']
rcParams['svg.fonttype'] = 'none'


# ============================================================
# Safe figure saving (avoids the FreeType raster-overflow crash)
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
# Step 1: Read image and metadata
# ============================================================

def read_czi_dual(filepath, channel_x_index, channel_y_index):
    """Read a .czi file and extract both channels + metadata."""
    import czifile
    with czifile.CziFile(filepath) as czi:
        data = czi.asarray().squeeze()
        metadata = czi.metadata()

    # Extract voxel sizes from metadata
    import re
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

    channel_x = data[channel_x_index].astype(np.float32)
    channel_y = data[channel_y_index].astype(np.float32)
    return channel_x, channel_y, (dz, dy, dx), na, immersion_ri, immersion


# ============================================================
# Step 2: Background subtraction
# ============================================================

def background_subtract(stack):
    """Subtract stack-wide median background."""
    bg = np.median(stack)
    return np.clip(stack - bg, 0, None), bg


# ============================================================
# Step 3: Theoretical PSF generation
# ============================================================

def theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength=607, psf_shape=(15, 31, 31)):
    """
    Generate a widefield Gaussian PSF.
    emission_wavelength in nm (default 607 for mRFP1.2).
    """
    lam = emission_wavelength / 1000.0  # um
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
# Step 4: Richardson-Lucy deconvolution
# ============================================================

def deconvolve(stack, psf, iterations=40):
    """Richardson-Lucy deconvolution via skimage."""
    stack_norm = stack / stack.max()
    deconv = richardson_lucy(stack_norm, psf, num_iter=iterations)
    return (deconv * stack.max()).astype(np.float32)


# ============================================================
# Steps 5-9: Segmentation
# ============================================================

def segment_nuclei(
    stack,
    voxel_size=(0.2, 0.10238, 0.10238),
    smooth_sigma=(1.0, 2.0, 2.0),
    seed_threshold=8.0,
    min_voxels=500,
    max_z_span_slices=25,
    watershed_min_distance=10,
    otsu_dilation_iters=3,
):
    """
    Segment yeast nuclei: smoothing -> seeds -> watershed -> per-object Otsu.
    Returns labeled array, smoothed stack, n_objects, qc dict.
    """
    nz = stack.shape[0]
    smoothed = gaussian_filter(stack.astype(np.float32), sigma=smooth_sigma)

    # Seed detection
    binary_seed = smoothed > seed_threshold
    labeled_seed, n_seed = ndi.label(binary_seed)
    sizes_seed = ndi.sum(binary_seed, labeled_seed, range(1, n_seed + 1)).astype(int)

    # Size + z-span filter
    keep_ids = []
    for obj_id in np.where(sizes_seed > min_voxels)[0] + 1:
        mask = labeled_seed == obj_id
        z_idx = np.where(mask.any(axis=(1, 2)))[0]
        z_span = int(z_idx.max() - z_idx.min() + 1)
        if z_span <= max_z_span_slices:
            keep_ids.append(obj_id)

    binary_filtered = np.isin(labeled_seed, keep_ids)

    # Watershed
    distance = ndi.distance_transform_edt(binary_filtered)
    coords = peak_local_max(distance, min_distance=watershed_min_distance, labels=binary_filtered)
    mask_peaks = np.zeros(binary_filtered.shape, dtype=bool)
    mask_peaks[tuple(coords.T)] = True
    markers, _ = ndi.label(mask_peaks)
    labeled_ws = watershed(-distance, markers, mask=binary_filtered)

    # Per-object adaptive Otsu
    n_ws = labeled_ws.max()
    sizes_ws = ndi.sum(binary_filtered, labeled_ws, range(1, n_ws + 1)).astype(int)

    final_labels = np.zeros_like(labeled_ws, dtype=np.int32)
    qc = {}
    new_id = 0

    for obj_id in range(1, n_ws + 1):
        if sizes_ws[obj_id - 1] < 300:
            continue
        mask_ws = labeled_ws == obj_id
        region = ndi.binary_dilation(mask_ws, iterations=otsu_dilation_iters)
        signal_vals = smoothed[region]
        signal_vals = signal_vals[signal_vals > 0]

        if len(signal_vals) < 50 or signal_vals.max() < 3:
            continue

        try:
            t_otsu = threshold_otsu(signal_vals)
        except Exception:
            t_otsu = signal_vals.mean()

        final_mask = (smoothed > t_otsu) & region
        final_size = final_mask.sum()
        if final_size < 200:
            continue

        z_idx_final = np.where(final_mask.any(axis=(1, 2)))[0]
        touches_boundary = (z_idx_final.min() == 0 or z_idx_final.max() == nz - 1) if len(z_idx_final) > 0 else False

        new_id += 1
        final_labels[final_mask] = new_id
        qc[new_id] = {
            'otsu_thr': float(t_otsu),
            'final_voxels': int(final_size),
            'touches_boundary': touches_boundary,
        }

    # ---- Remove nuclei touching any image border (XY + Z) ----
    # Kernels that touch z=0/z=nz-1, y=0/y=ny-1, or x=0/x=nx-1 are geometrically
    # incomplete (cut off by the field of view). Drop them BEFORE measurement so
    # they don't appear in all_nuclei.xlsx or fig4_annotated_* overlays.
    nz_dim, ny_dim, nx_dim = final_labels.shape
    n_before = new_id
    kept_ids = []
    for lid in range(1, new_id + 1):
        mask = final_labels == lid
        if not mask.any():
            continue
        on_border = (
            mask[0, :, :].any() or mask[nz_dim - 1, :, :].any() or   # Z
            mask[:, 0, :].any() or mask[:, ny_dim - 1, :].any() or   # Y
            mask[:, :, 0].any() or mask[:, :, nx_dim - 1].any()      # X
        )
        if not on_border:
            kept_ids.append(lid)

    # Re-label kept nuclei to 1..N without gaps and keep qc dict in sync
    relabeled = np.zeros_like(final_labels)
    new_qc = {}
    for new_lid, old_lid in enumerate(kept_ids, start=1):
        relabeled[final_labels == old_lid] = new_lid
        if old_lid in qc:
            new_qc[new_lid] = qc[old_lid]
    final_labels = relabeled
    qc = new_qc
    qc['_n_removed_border'] = n_before - len(kept_ids)  # sentinel for logging
    new_id = len(kept_ids)

    return final_labels, smoothed, new_id, qc


# ============================================================
# Step 10: Object matching between channels
# ============================================================

def match_objects(labels_nuc, labels_his):
    """
    Match histone objects to nuclear objects using containment criterion.

    A histone object is matched to a nuclear object if the histone's centroid
    lies within the nuclear mask. This is biologically motivated: histones
    should be localized inside the nucleus.

    If multiple histone centroids fall within the same nucleus, the one with
    the largest overlap volume is chosen. If a histone centroid falls outside
    all nuclei, it is not matched.

    Returns list of (nuc_id, his_id, overlap_fraction) tuples.
    overlap_fraction = fraction of histone voxels that overlap with the matched
    nuclear object (can be >1.0 if histone extends beyond nucleus).
    """
    pairs = []
    n_nuc = int(labels_nuc.max())
    n_his = int(labels_his.max())

    if n_nuc == 0 or n_his == 0:
        return pairs

    # Pre-compute nuclear masks and centroids
    nuc_masks = {}
    nuc_centroids = {}
    for nuc_id in range(1, n_nuc + 1):
        mask = labels_nuc == nuc_id
        if mask.sum() == 0:
            continue
        nuc_masks[nuc_id] = mask
        nuc_centroids[nuc_id] = ndi.center_of_mass(mask)

    # For each histone object, check if its centroid is inside any nucleus
    for his_id in range(1, n_his + 1):
        mask_his = labels_his == his_id
        n_his_voxels = mask_his.sum()
        if n_his_voxels == 0:
            continue

        com_his = ndi.center_of_mass(mask_his)
        com_his_int = tuple(int(round(c)) for c in com_his)

        # Find which nucleus contains the histone centroid
        containing_nuc = None
        for nuc_id, mask_nuc in nuc_masks.items():
            if mask_nuc[com_his_int]:
                containing_nuc = nuc_id
                break

        if containing_nuc is None:
            continue

        # Calculate overlap fraction (histone voxels inside nucleus / total histone voxels)
        mask_nuc = nuc_masks[containing_nuc]
        overlap_voxels = (mask_his & mask_nuc).sum()
        overlap_frac = overlap_voxels / n_his_voxels

        pairs.append((containing_nuc, his_id, float(overlap_frac)))

    return pairs


# ============================================================
# Step 11: Measurement
# ============================================================

def measure_matched_pairs(pairs, labels_nuc, labels_his, smoothed_nuc, smoothed_his,
                          voxel_size, sample_ri=1.35, immersion_ri=1.518,
                          qc_nuc=None, qc_his=None):
    """Measure volume, shape, and intensity for each matched pair with RI correction."""
    axial_correction = sample_ri / immersion_ri
    dz_c = voxel_size[0] * axial_correction
    dy, dx = voxel_size[1], voxel_size[2]
    voxel_vol = dz_c * dy * dx

    results = []

    for pair_idx, (nuc_id, his_id, overlap_frac) in enumerate(pairs, start=1):
        mask_nuc = labels_nuc == nuc_id
        mask_his = labels_his == his_id

        # Nuclear measurements
        n_voxels_nuc = int(mask_nuc.sum())
        vol_voxel_nuc = n_voxels_nuc * voxel_vol

        spacing = (dz_c, dy, dx)
        try:
            verts, faces, _, _ = marching_cubes(mask_nuc.astype(np.float32), level=0.5, spacing=spacing)
            v_mesh_nuc = abs(sum(np.dot(verts[f[0]], np.cross(verts[f[1]], verts[f[2]])) / 6.0 for f in faces))
            surf_area_nuc = mesh_surface_area(verts, faces)
            sphericity_nuc = (np.pi ** (1/3) * (6 * v_mesh_nuc) ** (2/3)) / surf_area_nuc if surf_area_nuc > 0 and v_mesh_nuc > 0 else 0
        except Exception:
            v_mesh_nuc, surf_area_nuc, sphericity_nuc = 0, 0, 0

        z_idx_nuc = np.where(mask_nuc.any(axis=(1, 2)))[0]
        z_span_um_nuc = (z_idx_nuc.max() - z_idx_nuc.min() + 1) * dz_c if len(z_idx_nuc) > 0 else 0
        mip_obj_nuc = mask_nuc.any(axis=0)
        xy_diam_nuc = 2 * np.sqrt(mip_obj_nuc.sum() / np.pi) * dx if mip_obj_nuc.sum() > 0 else 0
        z_xy_ratio_nuc = z_span_um_nuc / xy_diam_nuc if xy_diam_nuc > 0 else 0

        intensities_nuc = smoothed_nuc[mask_nuc]
        total_int_nuc = float(intensities_nuc.sum())

        # Histone measurements
        n_voxels_his = int(mask_his.sum())
        vol_voxel_his = n_voxels_his * voxel_vol

        try:
            verts, faces, _, _ = marching_cubes(mask_his.astype(np.float32), level=0.5, spacing=spacing)
            v_mesh_his = abs(sum(np.dot(verts[f[0]], np.cross(verts[f[1]], verts[f[2]])) / 6.0 for f in faces))
            surf_area_his = mesh_surface_area(verts, faces)
            sphericity_his = (np.pi ** (1/3) * (6 * v_mesh_his) ** (2/3)) / surf_area_his if surf_area_his > 0 and v_mesh_his > 0 else 0
        except Exception:
            v_mesh_his, surf_area_his, sphericity_his = 0, 0, 0

        z_idx_his = np.where(mask_his.any(axis=(1, 2)))[0]
        z_span_um_his = (z_idx_his.max() - z_idx_his.min() + 1) * dz_c if len(z_idx_his) > 0 else 0
        mip_obj_his = mask_his.any(axis=0)
        xy_diam_his = 2 * np.sqrt(mip_obj_his.sum() / np.pi) * dx if mip_obj_his.sum() > 0 else 0
        z_xy_ratio_his = z_span_um_his / xy_diam_his if xy_diam_his > 0 else 0

        intensities_his = smoothed_his[mask_his]
        total_int_his = float(intensities_his.sum())

        # Ratio
        ratio_percent = (vol_voxel_his / vol_voxel_nuc * 100) if vol_voxel_nuc > 0 else 0

        # QC flags
        flags = []
        if sphericity_nuc < 0.70: flags.append('low_sphericity_nuclear')
        if sphericity_his < 0.70: flags.append('low_sphericity_histone')
        if z_xy_ratio_nuc > 2.5: flags.append('high_z_xy_ratio_nuclear')
        if z_xy_ratio_his > 2.5: flags.append('high_z_xy_ratio_histone')
        if z_idx_nuc.min() == 0 or z_idx_nuc.max() == labels_nuc.shape[0] - 1: flags.append('touches_z_boundary_nuclear')
        if z_idx_his.min() == 0 or z_idx_his.max() == labels_his.shape[0] - 1: flags.append('touches_z_boundary_histone')
        if vol_voxel_nuc < 1.0: flags.append('very_small_nuclear')
        if vol_voxel_nuc > 20.0: flags.append('very_large_nuclear')
        if vol_voxel_his < 0.5: flags.append('very_small_histone')
        if vol_voxel_his > 15.0: flags.append('very_large_histone')
        if overlap_frac < 0.5: flags.append('low_containment')
        if ratio_percent > 100: flags.append('ratio_over_100')
        if ratio_percent < 10: flags.append('ratio_under_10')

        seg_qc_nuc = qc_nuc.get(nuc_id, {}) if qc_nuc else {}
        seg_qc_his = qc_his.get(his_id, {}) if qc_his else {}
        com_nuc = ndi.center_of_mass(mask_nuc)
        com_his = ndi.center_of_mass(mask_his)

        results.append({
            'pair_id': pair_idx,
            'nuclear_object_id': nuc_id,
            'histone_object_id': his_id,
            'nuclear_volume_um3': round(vol_voxel_nuc, 4),
            'histone_volume_um3': round(vol_voxel_his, 4),
            'ratio_percent': round(ratio_percent, 4),
            'nuclear_mesh_volume_um3': round(v_mesh_nuc, 4),
            'histone_mesh_volume_um3': round(v_mesh_his, 4),
            'nuclear_surface_area_um2': round(surf_area_nuc, 4),
            'histone_surface_area_um2': round(surf_area_his, 4),
            'nuclear_sphericity': round(sphericity_nuc, 4),
            'histone_sphericity': round(sphericity_his, 4),
            'overlap_fraction': round(overlap_frac, 4),
            'n_voxels_nuclear': n_voxels_nuc,
            'n_voxels_histone': n_voxels_his,
            'nuclear_z_span_um': round(z_span_um_nuc, 4),
            'histone_z_span_um': round(z_span_um_his, 4),
            'nuclear_xy_equiv_diameter_um': round(xy_diam_nuc, 4),
            'histone_xy_equiv_diameter_um': round(xy_diam_his, 4),
            'nuclear_z_xy_ratio': round(z_xy_ratio_nuc, 4),
            'histone_z_xy_ratio': round(z_xy_ratio_his, 4),
            'nuclear_mean_intensity': round(float(intensities_nuc.mean()), 2),
            'histone_mean_intensity': round(float(intensities_his.mean()), 2),
            'nuclear_max_intensity': round(float(intensities_nuc.max()), 2),
            'histone_max_intensity': round(float(intensities_his.max()), 2),
            'nuclear_median_intensity': round(float(np.median(intensities_nuc)), 2),
            'histone_median_intensity': round(float(np.median(intensities_his)), 2),
            'nuclear_total_intensity': round(total_int_nuc, 2),
            'histone_total_intensity': round(total_int_his, 2),
            'nuclear_centroid_z_um': round(com_nuc[0] * dz_c, 4),
            'nuclear_centroid_y_um': round(com_nuc[1] * dy, 4),
            'nuclear_centroid_x_um': round(com_nuc[2] * dx, 4),
            'histone_centroid_z_um': round(com_his[0] * dz_c, 4),
            'histone_centroid_y_um': round(com_his[1] * dy, 4),
            'histone_centroid_x_um': round(com_his[2] * dx, 4),
            'nuclear_otsu_threshold': seg_qc_nuc.get('otsu_thr'),
            'histone_otsu_threshold': seg_qc_his.get('otsu_thr'),
            'qc_flags': ';'.join(flags) if flags else 'none',
            'axial_correction_factor': round(axial_correction, 4),
        })
    return results


# ============================================================
# Step 12: QC figures
# ============================================================

def _scatter_regplot(ax, x, y, scatter_color, line_color,
                     xlabel, ylabel, title_fmt, p_fmt='.3f', n_nuc=None):
    """
    Draw scatter + linear regression line + Pearson r/p in title.

    Guarded against low-n cases:
      - n < 2: no scatter/regression; placeholder text; returns (nan, nan)
      - x or y is constant: pearsonr is undefined -> returns (nan, nan),
        scatter is drawn but no regression line

    `title_fmt` is a format string accepting `r` and `p` placeholders using
    Python .format(), so callers control the exact title layout (matches the
    original per-figure titles bit-identically when n >= 2 with variation).

    Returns (r, p) tuple.
    """
    if n_nuc is None:
        n_nuc = len(x)

    # Guard 1: need >=2 points for correlation/fit
    if n_nuc < 2:
        ax.text(0.5, 0.5, f'Insufficient nuclei (n={n_nuc})\nfor correlation analysis',
                ha='center', va='center', transform=ax.transAxes, fontsize=11)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        return float('nan'), float('nan')

    ax.scatter(x, y, c=scatter_color, edgecolors='black', linewidth=0.5, s=60, alpha=0.8)

    # Guard 2: pearsonr and polyfit both require variation in x (and y)
    x_var = float(np.ptp(x))
    y_var = float(np.ptp(y))
    if x_var == 0 or y_var == 0:
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        return float('nan'), float('nan')

    r, pval = pearsonr(x, y)
    z = np.polyfit(x, y, 1); p = np.poly1d(z)
    xl = np.linspace(x.min(), x.max(), 100)
    ax.plot(xl, p(xl), line_color, linewidth=2)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(title_fmt.format(r=r, p=pval))
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    return r, pval


def generate_qc_figures(results, labels_nuc, labels_his, smoothed_nuc, smoothed_his,
                        image_id, output_dir, voxel_size, sample_ri, immersion_ri):
    """
    Generate QC figures for an image.

    Robust to very low nucleus counts (n=0, 1, or many): correlations,
    polyfits, and histogram bin counts are all guarded. When n<2 the
    correlation-related plots show a placeholder and pearson r/p is NaN.
    """
    dz_c = voxel_size[0] * (sample_ri / immersion_ri)
    voxel_vol = dz_c * voxel_size[1] * voxel_size[2]

    vols_nuc = np.array([r['nuclear_volume_um3'] for r in results])
    vols_his = np.array([r['histone_volume_um3'] for r in results])
    ratios = np.array([r['ratio_percent'] for r in results])
    mean_ints_nuc = np.array([r['nuclear_mean_intensity'] for r in results])
    mean_ints_his = np.array([r['histone_mean_intensity'] for r in results])
    total_ints_nuc = np.array([r['nuclear_total_intensity'] for r in results])
    total_ints_his = np.array([r['histone_total_intensity'] for r in results])
    overlaps = np.array([r['overlap_fraction'] for r in results])

    n_pairs = len(results)

    # Fig 1: Intensity vs volume for both channels
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    r_p_nuc, p_p_nuc = _scatter_regplot(
        axes[0, 0], mean_ints_nuc, vols_nuc, '#0279EE', '#FF9400',
        'Mean nuclear intensity (a.u.)', 'Nuclear volume (µm³)',
        'Nuclear: Mean intensity vs volume\n(artifact test: r={r:.3f}, p={p:.3f})',
        n_nuc=n_pairs,
    )
    r_t_nuc, p_t_nuc = _scatter_regplot(
        axes[0, 1], total_ints_nuc, vols_nuc, '#75A025', '#FF9400',
        'Total nuclear intensity (a.u.)', 'Nuclear volume (µm³)',
        'Nuclear: Total intensity vs volume\n(biological signal: r={r:.3f}, p={p:.1e})',
        n_nuc=n_pairs,
    )
    r_p_his, p_p_his = _scatter_regplot(
        axes[1, 0], mean_ints_his, vols_his, '#FD9BED', '#FF9400',
        'Mean histone intensity (a.u.)', 'Histone volume (µm³)',
        'Histone: Mean intensity vs volume\n(artifact test: r={r:.3f}, p={p:.3f})',
        n_nuc=n_pairs,
    )
    r_t_his, p_t_his = _scatter_regplot(
        axes[1, 1], total_ints_his, vols_his, '#E9ED4C', '#FF9400',
        'Total histone intensity (a.u.)', 'Histone volume (µm³)',
        'Histone: Total intensity vs volume\n(biological signal: r={r:.3f}, p={p:.1e})',
        n_nuc=n_pairs,
    )
    safe_tight_layout(fig)
    safe_savefig(fig, f'{output_dir}/fig1_intensity_vs_volume_{image_id}.svg', format='svg')
    safe_savefig(fig, f'{output_dir}/fig1_intensity_vs_volume_{image_id}.png', format='png', dpi=150)
    plt.close()

    # Fig 2: Threshold sensitivity + distribution for both channels
    thresholds = [6, 8, 10, 12, 15, 18, 20]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Nuclear channel
    sens_med_nuc, sens_n_nuc = [], []
    for thr in thresholds:
        binary = smoothed_nuc > thr
        labeled_t, n_t = ndi.label(binary)
        sizes_t = ndi.sum(binary, labeled_t, range(1, n_t + 1)).astype(int)
        vols_t = sizes_t[sizes_t > 500] * voxel_vol
        sens_med_nuc.append(np.median(vols_t) if len(vols_t) > 0 else 0)
        sens_n_nuc.append(len(vols_t))

    ax = axes[0, 0]
    ax.plot(thresholds, sens_med_nuc, 'o-', c='#0279EE', linewidth=2, markersize=8)
    ax2 = ax.twinx()
    ax2.plot(thresholds, sens_n_nuc, 's--', c='#FF9400', linewidth=1.5, markersize=6, alpha=0.7)
    ax.set_xlabel('Threshold (a.u.)'); ax.set_ylabel('Median volume (µm³)', color='#0279EE')
    ax2.set_ylabel('Number of objects', color='#FF9400')
    ax.set_title('Nuclear: Threshold sensitivity')
    ax.spines['top'].set_visible(False); ax2.spines['top'].set_visible(False)

    ax = axes[0, 1]
    if n_pairs > 0:
        n_bins = min(12, max(1, n_pairs))
        ax.hist(vols_nuc, bins=n_bins, color='#0279EE', edgecolor='black', alpha=0.8)
        ax.axvline(np.median(vols_nuc), color='#FF9400', linestyle='--', linewidth=2,
                   label=f'Median = {np.median(vols_nuc):.2f} µm³')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No matched pairs', ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
    ax.set_xlabel('Nuclear volume (µm³)'); ax.set_ylabel('Count')
    ax.set_title(f'Nuclear volume distribution (n={n_pairs})')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # Histone channel
    sens_med_his, sens_n_his = [], []
    for thr in thresholds:
        binary = smoothed_his > thr
        labeled_t, n_t = ndi.label(binary)
        sizes_t = ndi.sum(binary, labeled_t, range(1, n_t + 1)).astype(int)
        vols_t = sizes_t[sizes_t > 500] * voxel_vol
        sens_med_his.append(np.median(vols_t) if len(vols_t) > 0 else 0)
        sens_n_his.append(len(vols_t))

    ax = axes[1, 0]
    ax.plot(thresholds, sens_med_his, 'o-', c='#FD9BED', linewidth=2, markersize=8)
    ax2 = ax.twinx()
    ax2.plot(thresholds, sens_n_his, 's--', c='#FF9400', linewidth=1.5, markersize=6, alpha=0.7)
    ax.set_xlabel('Threshold (a.u.)'); ax.set_ylabel('Median volume (µm³)', color='#FD9BED')
    ax2.set_ylabel('Number of objects', color='#FF9400')
    ax.set_title('Histone: Threshold sensitivity')
    ax.spines['top'].set_visible(False); ax2.spines['top'].set_visible(False)

    ax = axes[1, 1]
    if n_pairs > 0:
        n_bins = min(12, max(1, n_pairs))
        ax.hist(vols_his, bins=n_bins, color='#FD9BED', edgecolor='black', alpha=0.8)
        ax.axvline(np.median(vols_his), color='#FF9400', linestyle='--', linewidth=2,
                   label=f'Median = {np.median(vols_his):.2f} µm³')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No matched pairs', ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
    ax.set_xlabel('Histone volume (µm³)'); ax.set_ylabel('Count')
    ax.set_title(f'Histone volume distribution (n={n_pairs})')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    safe_tight_layout(fig)
    safe_savefig(fig, f'{output_dir}/fig2_threshold_sensitivity_{image_id}.svg', format='svg')
    safe_savefig(fig, f'{output_dir}/fig2_threshold_sensitivity_{image_id}.png', format='png', dpi=150)
    plt.close()

    # Fig 3: Segmentation overlay for both channels
    mip_nuc = smoothed_nuc.max(axis=0)
    mip_his = smoothed_his.max(axis=0)
    mip_labels_nuc = labels_nuc.max(axis=0)
    mip_labels_his = labels_his.max(axis=0)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    axes[0, 0].imshow(mip_nuc, cmap='gray', vmin=0, vmax=np.percentile(mip_nuc, 99.5))
    axes[0, 0].set_title('Nuclear channel (MIP)'); axes[0, 0].axis('off')

    from matplotlib.colors import ListedColormap
    n_labels_nuc = int(mip_labels_nuc.max())
    cmap_nuc = ListedColormap(plt.cm.tab20(np.linspace(0, 1, max(n_labels_nuc, 20))))
    axes[0, 1].imshow(mip_nuc, cmap='gray', vmin=0, vmax=np.percentile(mip_nuc, 99.5))
    axes[0, 1].imshow(np.ma.masked_where(mip_labels_nuc == 0, mip_labels_nuc), cmap=cmap_nuc, alpha=0.4)
    axes[0, 1].set_title(f'Nuclear segmentation (n={n_labels_nuc})'); axes[0, 1].axis('off')

    axes[1, 0].imshow(mip_his, cmap='gray', vmin=0, vmax=np.percentile(mip_his, 99.5))
    axes[1, 0].set_title('Histone channel (MIP)'); axes[1, 0].axis('off')

    n_labels_his = int(mip_labels_his.max())
    cmap_his = ListedColormap(plt.cm.tab20(np.linspace(0, 1, max(n_labels_his, 20))))
    axes[1, 1].imshow(mip_his, cmap='gray', vmin=0, vmax=np.percentile(mip_his, 99.5))
    axes[1, 1].imshow(np.ma.masked_where(mip_labels_his == 0, mip_labels_his), cmap=cmap_his, alpha=0.4)
    axes[1, 1].set_title(f'Histone segmentation (n={n_labels_his})'); axes[1, 1].axis('off')

    safe_tight_layout(fig)
    safe_savefig(fig, f'{output_dir}/fig3_segmentation_overlay_{image_id}.svg', format='svg')
    safe_savefig(fig, f'{output_dir}/fig3_segmentation_overlay_{image_id}.png', format='png', dpi=150)
    plt.close()

    # Fig 4: Annotated map with pair IDs and ratios
    import matplotlib.patheffects

    fig, ax = plt.subplots(figsize=(14, 12))
    ax.imshow(mip_nuc, cmap='gray', vmin=0, vmax=np.percentile(mip_nuc, 99.5))
    ax.imshow(np.ma.masked_where(mip_labels_nuc == 0, mip_labels_nuc), cmap=cmap_nuc, alpha=0.35)

    ratio_lookup = {r['nuclear_object_id']: r['ratio_percent'] for r in results}
    pair_lookup = {r['nuclear_object_id']: r['pair_id'] for r in results}

    for nuc_id in range(1, n_labels_nuc + 1):
        mask_3d = labels_nuc == nuc_id
        mip_obj = mask_3d.any(axis=0)
        if mip_obj.sum() == 0:
            continue
        com_y, com_x = ndi.center_of_mass(mip_obj)
        pair_id = pair_lookup.get(nuc_id, 0)
        ratio = ratio_lookup.get(nuc_id, 0)
        if pair_id > 0:
            ax.text(com_x, com_y - 15, f'#{pair_id}', color='yellow', fontsize=7,
                    fontweight='bold', ha='center', va='center',
                    path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground='black')])
            ax.text(com_x, com_y + 15, f'{ratio:.1f}%', color='cyan', fontsize=6,
                    ha='center', va='center',
                    path_effects=[matplotlib.patheffects.withStroke(linewidth=1.5, foreground='black')])

    ax.set_title(f'{image_id}: Matched pairs (n={n_pairs})\n'
                 f'Yellow = pair ID, Cyan = histone/nuclear ratio (%)', fontsize=13)
    ax.axis('off')
    safe_tight_layout(fig)
    safe_savefig(fig, f'{output_dir}/fig4_annotated_map_{image_id}.svg', format='svg')
    safe_savefig(fig, f'{output_dir}/fig4_annotated_map_{image_id}.png', format='png', dpi=150)
    plt.close()

    # Fig 5: Ratio distribution + overlap quality
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    if n_pairs > 0:
        n_bins = min(12, max(1, n_pairs))
        ax.hist(ratios, bins=n_bins, color='#75A025', edgecolor='black', alpha=0.8)
        ax.axvline(np.median(ratios), color='#FF9400', linestyle='--', linewidth=2,
                   label=f'Median = {np.median(ratios):.1f}%')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No matched pairs', ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
    ax.set_xlabel('Histone/Nuclear volume ratio (%)'); ax.set_ylabel('Count')
    ax.set_title(f'Ratio distribution (n={n_pairs})')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    ax = axes[1]
    if n_pairs > 0:
        ax.scatter(overlaps, ratios, c='#0279EE', edgecolors='black', linewidth=0.5, s=60, alpha=0.8)
        ax.axhline(np.median(ratios), color='#FF9400', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(0.5, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label='Containment threshold')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No matched pairs', ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
    ax.set_xlabel('Containment fraction (histone in nucleus)'); ax.set_ylabel('Histone/Nuclear ratio (%)')
    ax.set_title('Containment quality vs ratio')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    safe_tight_layout(fig)
    safe_savefig(fig, f'{output_dir}/fig5_ratio_distribution_{image_id}.svg', format='svg')
    safe_savefig(fig, f'{output_dir}/fig5_ratio_distribution_{image_id}.png', format='png', dpi=150)
    plt.close()

    return {
        'intensity_vol_r_nuclear': r_p_nuc,
        'total_int_vol_r_nuclear': r_t_nuc,
        'intensity_vol_r_histone': r_p_his,
        'total_int_vol_r_histone': r_t_his,
    }


# ============================================================
# Main processing function
# ============================================================

def process_image(filepath, output_dir, channel_x_index, channel_y_index,
                  sample_ri=1.35, emission_wavelength_x=607, emission_wavelength_y=509,
                  deconv_iterations=40, strain=None):
    """Process a single image through the full dual-channel pipeline."""
    image_id = os.path.splitext(os.path.basename(filepath))[0]
    print(f"\n{'='*60}")
    print(f"Processing: {image_id}")
    print(f"{'='*60}")

    t_start = time.time()

    # Step 1: Read both channels
    if filepath.lower().endswith('.czi'):
        channel_x, channel_y, voxel_size, na, immersion_ri, immersion = read_czi_dual(
            filepath, channel_x_index, channel_y_index)
    else:
        raise ValueError(f"Only .czi files supported, got: {filepath}")

    print(f"  Channel X (nuclear) shape: {channel_x.shape}, voxel: {voxel_size}, NA={na}, immersion_ri={immersion_ri}")
    print(f"  Channel Y (histone) shape: {channel_y.shape}")

    # Step 2: Background subtract both channels
    channel_x_bs, bg_x = background_subtract(channel_x)
    channel_y_bs, bg_y = background_subtract(channel_y)
    print(f"  Background (median): nuclear={bg_x:.1f}, histone={bg_y:.1f}")

    # Step 3: PSF for both channels
    psf_x = theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength_x)
    psf_y = theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength_y)
    print(f"  PSF nuclear: {psf_x.shape}, histone: {psf_y.shape}")

    # Step 4: Deconvolve both channels
    print(f"  Deconvolving nuclear channel ({deconv_iterations} iterations)...")
    t0 = time.time()
    deconv_x = deconvolve(channel_x_bs, psf_x, iterations=deconv_iterations)
    print(f"  Nuclear deconvolution: {time.time()-t0:.1f}s")

    print(f"  Deconvolving histone channel ({deconv_iterations} iterations)...")
    t0 = time.time()
    deconv_y = deconvolve(channel_y_bs, psf_y, iterations=deconv_iterations)
    print(f"  Histone deconvolution: {time.time()-t0:.1f}s")

    # Steps 5-9: Segment both channels
    print(f"  Segmenting nuclear channel...")
    labels_nuc, smoothed_nuc, n_nuc, qc_nuc = segment_nuclei(deconv_x, voxel_size=voxel_size)
    n_removed_border_nuc = qc_nuc.pop('_n_removed_border', 0)
    if n_removed_border_nuc > 0:
        print(f"  Removed {n_removed_border_nuc} nuclear objects touching image border (kept {n_nuc})")
    print(f"  Found {n_nuc} nuclear objects")

    print(f"  Segmenting histone channel...")
    labels_his, smoothed_his, n_his, qc_his = segment_nuclei(deconv_y, voxel_size=voxel_size)
    n_removed_border_his = qc_his.pop('_n_removed_border', 0)
    if n_removed_border_his > 0:
        print(f"  Removed {n_removed_border_his} histone objects touching image border (kept {n_his})")
    print(f"  Found {n_his} histone objects")

    # Step 10: Match objects between channels (containment-based)
    print(f"  Matching objects between channels (containment criterion)...")
    pairs = match_objects(labels_nuc, labels_his)
    print(f"  Found {len(pairs)} matched pairs (from {n_nuc} nuclear + {n_his} histone objects)")

    # Step 11: Measure matched pairs
    raw_results = measure_matched_pairs(
        pairs, labels_nuc, labels_his, smoothed_nuc, smoothed_his,
        voxel_size, sample_ri, immersion_ri, qc_nuc=qc_nuc, qc_his=qc_his)

    # Reorder columns: put 'strain' and 'image_id' first
    results = []
    for r in raw_results:
        new_r = {'strain': strain if strain else '', 'image_id': image_id}
        for k, v in r.items():
            if k not in ('strain', 'image_id'):
                new_r[k] = v
        results.append(new_r)
    print(f"  Measured {len(results)} matched pairs")

    # Step 12: QC figures
    print(f"  Generating QC figures...")
    qc_stats = generate_qc_figures(results, labels_nuc, labels_his, smoothed_nuc, smoothed_his,
                                    image_id, output_dir, voxel_size, sample_ri, immersion_ri)

    # Build QC metrics dict
    ratios = np.array([r['ratio_percent'] for r in results])
    vols_nuc = np.array([r['nuclear_volume_um3'] for r in results])
    vols_his = np.array([r['histone_volume_um3'] for r in results])

    qc_metrics = {
        'image_id': image_id,
        'strain': strain or '',
        'n_matched_pairs': len(results),
        'n_nuclear_objects': n_nuc,
        'n_histone_objects': n_his,
        'ratio_median_percent': float(np.median(ratios)) if len(ratios) > 0 else 0,
        'ratio_mean_percent': float(ratios.mean()) if len(ratios) > 0 else 0,
        'ratio_std_percent': float(ratios.std()) if len(ratios) > 0 else 0,
        'nuclear_volume_median_um3': float(np.median(vols_nuc)) if len(vols_nuc) > 0 else 0,
        'histone_volume_median_um3': float(np.median(vols_his)) if len(vols_his) > 0 else 0,
        'intensity_vol_pearson_r_nuclear': qc_stats['intensity_vol_r_nuclear'],
        'total_int_vol_pearson_r_nuclear': qc_stats['total_int_vol_r_nuclear'],
        'intensity_vol_pearson_r_histone': qc_stats['intensity_vol_r_histone'],
        'total_int_vol_pearson_r_histone': qc_stats['total_int_vol_r_histone'],
        'axial_correction_factor': sample_ri / immersion_ri,
        'deconvolution': f'RL {deconv_iterations} iter, Gaussian PSF',
    }

    elapsed = time.time() - t_start
    print(f"  Done in {elapsed:.1f}s")
    if len(ratios) > 0:
        print(f"  Ratio: median={np.median(ratios):.1f}%, "
              f"nuclear_vol={np.median(vols_nuc):.3f} µm³, "
              f"histone_vol={np.median(vols_his):.3f} µm³")
    else:
        print(f"  WARNING: 0 matched pairs (no overlapping objects found)")

    return results, qc_metrics


# ============================================================
# CLI
# ============================================================

CHANNEL_NAMES = {1: 'mRFP1.2', 2: 'EGFP'}

# Emission wavelength is a fixed property of the fluorophore, so it follows
# the channel choice automatically. No need to specify it on the command line.
#   Channel 1 = mRFP1.2 -> 607 nm emission maximum
#   Channel 2 = EGFP    -> 509 nm emission maximum
CHANNEL_EMISSION_NM = {1: 607, 2: 509}


def resolve_channel(cli_value, channel_name):
    """Return channel index (1 or 2). Ask interactively if cli_value is None.

    - If cli_value is not None it was already validated by argparse choices=[1, 2].
    - If stdin is not a TTY and cli_value is None, abort with a clear message
      instead of silently falling back (avoids wrong-channel results in cron/batch runs).
    - On empty input (Enter) the default is 1 for nuclear, 2 for histone.
    """
    if cli_value is not None:
        return cli_value

    if not sys.stdin.isatty():
        raise SystemExit(
            f"Error: --channel-{channel_name} is required in non-interactive mode.\n"
            f"Use --channel-{channel_name} 1 for mRFP1.2 or --channel-{channel_name} 2 for EGFP."
        )

    default = 1 if channel_name == 'x' else 2
    while True:
        raw = input(f"Welcher Kanal ist der {channel_name}-Kanal ({'Kern' if channel_name == 'x' else 'Histon'})? "
                    f"[1] mRFP1.2 / [2] EGFP (Enter={default}): ").strip()
        if raw == "":
            return default
        if raw in ("1", "2"):
            return int(raw)
        print(f"  Bitte 1 oder 2 eingeben (oder Enter fuer {default}).")


def main():
    parser = argparse.ArgumentParser(
        description='Yeast nuclear/histone volume ratio quantification pipeline')
    parser.add_argument('--input', required=True,
                        help='Input .czi file or folder of .czi files')
    parser.add_argument('--output', required=True,
                        help='Output directory for results')
    parser.add_argument('--channel-x', type=int, default=None, choices=[1, 2],
                        help='Nuclear marker channel: 1=mRFP1.2, 2=EGFP. If omitted, script asks interactively.')
    parser.add_argument('--channel-y', type=int, default=None, choices=[1, 2],
                        help='Histone marker channel: 1=mRFP1.2, 2=EGFP. If omitted, script asks interactively.')
    parser.add_argument('--sample-ri', type=float, default=1.35,
                        help='Sample refractive index (default: 1.35 yeast)')
    parser.add_argument('--deconv-iterations', type=int, default=40,
                        help='RL deconvolution iterations (default: 40)')
    parser.add_argument('--strain-map',
                        help='CSV mapping filenames to strain labels (columns: filename, strain)')
    parser.add_argument('--no-deconvolution', action='store_true',
                        help='Skip deconvolution (faster, less accurate)')
    args = parser.parse_args()

    # Resolve channels BEFORE any file processing so the log is clear from the start.
    args.channel_x = resolve_channel(args.channel_x, 'x')
    args.channel_y = resolve_channel(args.channel_y, 'y')

    if args.channel_x == args.channel_y:
        raise SystemExit(
            "Error: --channel-x and --channel-y must be different channels.\n"
            f"Both were set to channel {args.channel_x}."
        )

    # Emission wavelengths follow automatically from the channel choice
    # (607 nm for mRFP1.2 = channel 1, 509 nm for EGFP = channel 2).
    emission_wavelength_x = CHANNEL_EMISSION_NM[args.channel_x]
    emission_wavelength_y = CHANNEL_EMISSION_NM[args.channel_y]

    print(f"Nuclear channel: {args.channel_x} ({CHANNEL_NAMES[args.channel_x]}, "
          f"emission {emission_wavelength_x} nm)")
    print(f"Histone channel: {args.channel_y} ({CHANNEL_NAMES[args.channel_y]}, "
          f"emission {emission_wavelength_y} nm)")

    os.makedirs(args.output, exist_ok=True)

    # Build file list
    if os.path.isdir(args.input):
        files = sorted([os.path.join(args.input, f) for f in os.listdir(args.input)
                       if f.lower().endswith(('.czi', '.tif', '.tiff'))])
    else:
        files = [args.input]

    if not files:
        print(f"No .czi or .tif files found in {args.input}")
        sys.exit(1)

    # Load strain mapping
    strain_map = {}
    if args.strain_map:
        with open(args.strain_map) as f:
            reader = csv.DictReader(f)
            for row in reader:
                strain_map[row['filename']] = row['strain']

    # Failed-image log: initialize with header. Batch continues past crashes;
    # each failure is logged with error type, message, and timestamp so the
    # user can inspect + rerun those specific images later.
    failed_path = f'{args.output}/failed_images.tsv'
    with open(failed_path, 'w') as f:
        f.write('filename\terror_type\terror_message\ttimestamp\n')

    # Process all files with per-image try/except so an overnight batch never
    # aborts on a single bad image.
    all_results = []
    all_qc = []
    n_success = 0
    n_failed = 0
    for filepath in files:
        fname = os.path.basename(filepath)
        strain = strain_map.get(fname, strain_map.get(os.path.splitext(fname)[0], None))
        try:
            results, qc = process_image(
                filepath, args.output,
                channel_x_index=args.channel_x,
                channel_y_index=args.channel_y,
                sample_ri=args.sample_ri,
                emission_wavelength_x=emission_wavelength_x,
                emission_wavelength_y=emission_wavelength_y,
                deconv_iterations=args.deconv_iterations,
                strain=strain,
            )
            all_results.extend(results)
            all_qc.append(qc)
            n_success += 1
        except KeyboardInterrupt:
            # Never swallow Ctrl+C -- the user must be able to abort the batch.
            print(f"\n  [ABORT] KeyboardInterrupt during {fname}. Stopping batch.")
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
            # Explicit: continue to the next image (implicit anyway, but clear)
            continue

    # Save combined matched pairs table (XLSX)
    if all_results:
        combined_path = f'{args.output}/all_matched_pairs.xlsx'
        save_xlsx(all_results, combined_path, sheet_name='all_matched_pairs')
        print(f"\nCombined results: {combined_path} ({len(all_results)} matched pairs total)")

        # ============================================================
        # Wide-format per-strain table with 3 sheets (matches R output).
        # ============================================================
        # One workbook, three sheets:
        #   Sheet 'nuclear_volume_um3' = nuclear volumes
        #   Sheet 'histone_volume_um3' = histone volumes
        #   Sheet 'ratio_percent'      = histone/nuclear ratio (%)
        #
        # Columns = strain names (alphabetically sorted). Pairs without
        # a strain assignment go into an 'unassigned' column, placed at
        # its natural alphabetical position.
        #
        # ROW PAIRING: row i in all three sheets belongs to the SAME cell.
        # This is guaranteed by collecting (nuclear, histone, ratio) tuples
        # per strain in processing order and writing them identically into
        # all three sheets. Shorter strain columns are padded with empty
        # cells so the table is rectangular (padding only ever differs
        # BETWEEN strains, never between sheets of the same strain).
        pairs_by_strain = {}
        for r in all_results:
            key = r.get('strain') or 'unassigned'
            if key == '':
                key = 'unassigned'
            pairs_by_strain.setdefault(key, []).append(
                (r['nuclear_volume_um3'], r['histone_volume_um3'], r['ratio_percent'])
            )
        strain_order = sorted(pairs_by_strain.keys())
        max_n = max(len(v) for v in pairs_by_strain.values())

        # value_index: 0 = nuclear volume, 1 = histone volume, 2 = ratio
        sheet_specs = [
            ('nuclear_volume_um3', 0),
            ('histone_volume_um3', 1),
            ('ratio_percent', 2),
        ]

        wide_wb = openpyxl.Workbook()
        wide_wb.remove(wide_wb.active)  # drop default sheet
        for sheet_name, value_index in sheet_specs:
            ws = wide_wb.create_sheet(title=sheet_name)
            ws.append(strain_order)
            for row_idx in range(max_n):
                ws.append([
                    pairs_by_strain[s][row_idx][value_index]
                    if row_idx < len(pairs_by_strain[s]) else None
                    for s in strain_order
                ])
            # Format data cells as float with 4 decimals
            for col_idx, _ in enumerate(strain_order, start=1):
                col_letter = get_column_letter(col_idx)
                for row in range(2, max_n + 2):
                    ws[f'{col_letter}{row}'].number_format = '0.0000'

        wide_path = f'{args.output}/Volumes_per_strain.xlsx'
        wide_wb.save(wide_path)
        print(f"Volumes per strain (3 sheets): {wide_path} "
              f"({max_n} rows, {len(strain_order)} strains)")

    # Save combined QC (XLSX)
    if all_qc:
        combined_qc_path = f'{args.output}/all_qc_metrics.xlsx'
        save_xlsx(all_qc, combined_qc_path, sheet_name='all_qc_metrics')
        print(f"Combined QC: {combined_qc_path}")

    # Final summary: success vs failure counts and pointer to failed_images.tsv
    print(f"\nPipeline complete. {n_success} of {len(files)} images successful, "
          f"{n_failed} failed, {len(all_results)} matched pairs total.")
    if n_failed > 0:
        print(f"See {failed_path} for details of the {n_failed} failed image(s).")


if __name__ == '__main__':
    main()
