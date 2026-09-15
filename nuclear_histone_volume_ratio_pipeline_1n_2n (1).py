"""
Yeast Nuclear/Histone Volume Ratio Quantification Pipeline
==========================================================
Processes fluorescence z-stacks to measure yeast nuclear and histone
volumes separately, match corresponding objects, and compute the
histone/nuclear volume ratio per cell.

Supported input formats:
  - .czi (Zeiss widefield, both channels in one file) -- legacy thesis data
  - OME-TIFF (VisiView/Visitron spinning-disk confocal): each channel is a
    separate file named '<prefix>_w<N><tag>.ome.tf2'. Files are paired by
    their shared prefix; the channel is selected by a fluorophore-tag
    substring (e.g. --channel-x GFP --channel-y RFP). DIC files
    (tag contains 'dic') are skipped automatically.

Pipeline:
  1. Read image(s) and extract both channels + metadata
  2. Background subtraction (median) per channel
  3. Theoretical PSF generation per channel (different emission wavelengths;
     widefield or confocal model depending on --modality)
  4. Richardson-Lucy deconvolution per channel (skimage, 40 iterations)
  5. Axial RI correction (oil -> aqueous sample)
  6. 3D Gaussian smoothing per channel
  7. Seed detection + size/z-span filtering per channel
  8. Watershed separation of touching nuclei per channel
  9. Per-object adaptive Otsu thresholding per channel
  10. Object matching between channels (histone centroid containment in nucleus)
  11. Volume + shape measurement per matched pair
  12. QC: intensity-volume correlation, threshold sensitivity, ratio distribution

All size/scale segmentation parameters are in PHYSICAL units (um / um^3)
and converted to pixels per image, so the same settings work for different
microscopes and z-step sizes. The defaults reproduce the original Zeiss
widefield settings exactly.

Usage:
  # .czi (Zeiss widefield): channels are 1-based indices as in ZEN
  python nuclear_histone_volume_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x 1 --channel-y 2

  # OME-TIFF (spinning disk): channels are fluorophore tags from the filename
  python nuclear_histone_volume_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x GFP --channel-y RFP --seed-threshold 150

  # With strain mapping CSV (columns: filename, strain).
  # For OME-TIFF the shared prefix (e.g. W21653_1) also works as filename key.
  python nuclear_histone_volume_ratio_pipeline_1n_2n.py --input ./images/ --output ./results --channel-x GFP --channel-y RFP --strain-map strains.csv

Channel selection:
  .czi:      --channel-x/-y = 1 (mRFP1.2) or 2 (EGFP); emission wavelength
             follows automatically (607 nm / 509 nm).
  OME-TIFF:  --channel-x/-y = tag substring, e.g. GFP or RFP; emission
             wavelength is derived from the tag (override with
             --emission-x-nm / --emission-y-nm if needed).

Other options:
  --modality {widefield,confocal}  PSF model (default: auto = widefield for
                                   .czi, confocal for OME-TIFF)
  --seed-threshold FLOAT           absolute seed threshold on the smoothed
                                   stack (default 8.0, calibrated for the
                                   Zeiss widefield camera; use ~150 for the
                                   Visitron sdc uint16 data)
  --no-deconvolution               skip RL deconvolution (faster)
  --sample-ri FLOAT                sample refractive index (default 1.35)

Requirements:
  pip install czifile tifffile scikit-image scipy numpy matplotlib openpyxl cairosvg
"""

import argparse
import gc
import os
import re
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
# Step 1b: OME-TIFF input (VisiView / Visitron spinning disk)
# ============================================================

# Emission wavelengths (nm) by fluorophore name; matched as substring of the
# filename tag, case-insensitive. NOTE: the OME 'Wavelength' attribute is the
# EXCITATION laser, not the emission -- hence this lookup.
FLUOROPHORE_EMISSION_NM = {
    'GFP': 509, 'EGFP': 509,
    'RFP': 607, 'mRFP': 607, 'MCHERRY': 610, 'YEMRFP': 607,
}

OME_TIFF_EXTS = ('.ome.tif', '.ome.tiff', '.ome.tf2', '.ome.tf8', '.tif', '.tiff')


def split_ome_prefix(filename):
    """
    Split a VisiView OME-TIFF name '<prefix>_w<N><tag>.ome.tf2' into
    (prefix, tag). Returns (None, None) if the pattern does not match.
    """
    base = os.path.basename(filename)
    m = re.match(r'^(.*)_w(\d+)(.*)$', base)
    if not m:
        return None, None
    prefix = m.group(1)
    tag = m.group(3)
    # strip the extension(s) from the tag
    for ext in OME_TIFF_EXTS:
        if tag.lower().endswith(ext):
            tag = tag[:-len(ext)]
            break
    return prefix, tag


def read_ome_tiff_pair(filepath_x, filepath_y):
    """
    Read two single-channel OME-TIFF files (nuclear + histone) + metadata.

    Voxel sizes come from PhysicalSizeX/Y/Z (assumed um, the OME default
    when no unit attribute is present). NA comes from Objective LensNA,
    immersion RI from ObjectiveSettings RefractiveIndex.
    """
    import tifffile

    def _read_one(path):
        with tifffile.TiffFile(path) as tf:
            data = np.squeeze(tf.asarray()).astype(np.float32)
            ome = tf.ome_metadata or ''
            # ImageJ fallback for plain TIFFs saved by Fiji (no OME XML):
            # X/YResolution hold pixels per unit; unit 'micron' -> um/px = 1/res.
            # 'spacing' holds the ImageJ voxel depth (z step), same unit.
            ij_dx = ij_dy = ij_dz = None
            if not ome and tf.is_imagej:
                ij = tf.imagej_metadata or {}
                xt = tf.pages[0].tags.get('XResolution')
                yt = tf.pages[0].tags.get('YResolution')
                if xt is not None and yt is not None and ij.get('unit') in ('micron', 'um', 'µm'):
                    ij_dx = xt.value[1] / xt.value[0]
                    ij_dy = yt.value[1] / yt.value[0]
                ij_dz = ij.get('spacing')
        if data.ndim != 3:
            raise ValueError(f"Expected a 3D (ZYX) stack in {path}, got shape {data.shape}")

        def _f(pat, default):
            m = re.search(pat, ome)
            return float(m.group(1)) if m else default

        dx = _f(r'PhysicalSizeX="([\d.eE+-]+)"', ij_dx if ij_dx else 0.1098)
        dy = _f(r'PhysicalSizeY="([\d.eE+-]+)"', ij_dy if ij_dy else 0.1098)
        dz = _f(r'PhysicalSizeZ="([\d.eE+-]+)"', ij_dz if ij_dz else 1.0)
        m = re.search(r'LensNA="([\d.]+)"', ome)
        na = float(m.group(1)) if m else 1.4
        m = re.search(r'RefractiveIndex="([\d.]+)"', ome)
        immersion_ri = float(m.group(1)) if m else 1.518
        return data, (dz, dy, dx), na, immersion_ri

    channel_x, voxel_size, na, immersion_ri = _read_one(filepath_x)
    channel_y, _, _, _ = _read_one(filepath_y)
    return channel_x, channel_y, voxel_size, na, immersion_ri, 'Oil'


# ============================================================
# Step 2: Background subtraction
# ============================================================

def background_subtract(stack):
    """Subtract stack-wide median background.

    In-place: the input is already a fresh float32 array produced by the
    reader, so modifying it avoids two extra full-stack temporaries.
    """
    bg = np.median(stack)
    stack -= bg
    np.clip(stack, 0, None, out=stack)
    return stack, bg


# ============================================================
# Step 3: Theoretical PSF generation
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


def segment_nuclei(
    stack,
    voxel_size=(0.2, 0.10238, 0.10238),
    smooth_sigma_um=(0.2, 0.20476190476190476, 0.20476190476190476),
    seed_threshold=8.0,
    seed_k=5.0,
    min_volume_um3=1.048,
    max_z_span_um=15.0,
    watershed_min_distance_um=1.02,
    otsu_dilation_xy_um=0.307,
    otsu_dilation_z_um=0.6,
):
    """
    Segment yeast nuclei: smoothing -> seeds -> watershed -> per-object Otsu.
    Returns labeled array, smoothed stack, n_objects, qc dict.

    seed_threshold may be a number (absolute intensity on the smoothed stack)
    or the string 'auto' (noise-adaptive: median + seed_k * sigma_robust via
    adaptive_seed_threshold). The effective numeric threshold is stored in the
    qc dict under the '_seed_threshold_effective' sentinel key.

    All size/scale parameters are in PHYSICAL units (um / um^3) and are
    converted to pixels using the image's voxel size, so the same settings
    work for different microscopes and z-step sizes. The defaults reproduce
    the original pixel parameters of the widefield CZI data exactly
    (voxel 0.2 x 0.10238 x 0.10238 um: sigma (1,2,2) px, 500 voxels,
    25 z-slices, 10 px watershed distance, 3 px isotropic Otsu dilation).
    """
    dz, dy, dx = voxel_size
    voxel_vol = dz * dy * dx

    # Convert physical parameters to pixel units for THIS image
    smooth_sigma = tuple(s / v for s, v in zip(smooth_sigma_um, voxel_size))
    min_voxels = int(round(min_volume_um3 / voxel_vol))
    max_z_span_slices = max(1, int(round(max_z_span_um / dz)))
    watershed_min_distance = max(1, int(round(watershed_min_distance_um / dy)))
    otsu_dilation_iters_xy = max(1, int(round(otsu_dilation_xy_um / dy)))
    otsu_dilation_iters_z = max(0, int(round(otsu_dilation_z_um / dz)))
    # Post-watershed and final size cutoffs (originals: 300 / 200 voxels)
    min_ws_voxels = int(round(0.629 / voxel_vol))
    min_final_voxels = int(round(0.419 / voxel_vol))

    nz = stack.shape[0]
    smoothed = gaussian_filter(stack.astype(np.float32), sigma=smooth_sigma)

    # Seed detection: absolute threshold, or noise-adaptive when 'auto'
    if isinstance(seed_threshold, str):
        if seed_threshold.lower() != 'auto':
            raise ValueError(
                f"seed_threshold must be a number or 'auto', got {seed_threshold!r}")
        seed_thr_eff = adaptive_seed_threshold(smoothed, k=seed_k)
    else:
        seed_thr_eff = float(seed_threshold)
    binary_seed = smoothed > seed_thr_eff
    labeled_seed, n_seed = ndi.label(binary_seed)
    sizes_seed = ndi.sum(binary_seed, labeled_seed, range(1, n_seed + 1)).astype(int)
    del binary_seed

    # Size + z-span filter
    keep_ids = []
    for obj_id in np.where(sizes_seed > min_voxels)[0] + 1:
        mask = labeled_seed == obj_id
        z_idx = np.where(mask.any(axis=(1, 2)))[0]
        z_span = int(z_idx.max() - z_idx.min() + 1)
        if z_span <= max_z_span_slices:
            keep_ids.append(obj_id)

    binary_filtered = np.isin(labeled_seed, keep_ids)
    del labeled_seed

    # Watershed
    distance = ndi.distance_transform_edt(binary_filtered)
    coords = peak_local_max(distance, min_distance=watershed_min_distance, labels=binary_filtered)
    mask_peaks = np.zeros(binary_filtered.shape, dtype=bool)
    mask_peaks[tuple(coords.T)] = True
    markers, _ = ndi.label(mask_peaks)
    del mask_peaks
    distance *= -1  # watershed needs minima at object centers; in-place
                    # avoids a second full-size float64 copy
    labeled_ws = watershed(distance, markers, mask=binary_filtered)
    del distance, markers

    # Per-object adaptive Otsu
    n_ws = labeled_ws.max()
    sizes_ws = ndi.sum(binary_filtered, labeled_ws, range(1, n_ws + 1)).astype(int)
    del binary_filtered

    final_labels = np.zeros_like(labeled_ws, dtype=np.int32)
    qc = {}
    new_id = 0

    for obj_id in range(1, n_ws + 1):
        if sizes_ws[obj_id - 1] < min_ws_voxels:
            continue
        mask_ws = labeled_ws == obj_id
        # Anisotropic dilation: separate iteration counts for XY and Z.
        # When both are equal (e.g. the CZI voxel), use the classic
        # isotropic dilation to reproduce the original behavior exactly.
        if otsu_dilation_iters_z == otsu_dilation_iters_xy:
            region = ndi.binary_dilation(mask_ws, iterations=otsu_dilation_iters_xy)
        else:
            struct_xy = np.zeros((1, 3, 3), dtype=bool)
            struct_xy[0] = ndi.generate_binary_structure(2, 1)
            struct_z = np.zeros((3, 1, 1), dtype=bool)
            struct_z[:, 0, 0] = True
            region = ndi.binary_dilation(mask_ws, structure=struct_xy,
                                         iterations=otsu_dilation_iters_xy)
            region = ndi.binary_dilation(region, structure=struct_z,
                                         iterations=otsu_dilation_iters_z)
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
        if final_size < min_final_voxels:
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
    del labeled_ws

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
    qc['_seed_threshold_effective'] = seed_thr_eff  # sentinel for logging/fig2
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

    # Memory-efficient containment: read the nuclear label directly at each
    # histone centroid and compute overlaps inside the histone bounding box.
    # Numerically identical to materializing one full-size boolean mask per
    # object (labels are disjoint, so at most one nuclear label is nonzero at
    # any voxel) but avoids n_objects x stack-size bool arrays, which OOMed
    # on large OME stacks (142 nuclei x 340 MB on an 81x2048x2048 stack).
    his_slices = ndi.find_objects(labels_his)
    for his_id in range(1, n_his + 1):
        sl = his_slices[his_id - 1]
        if sl is None:
            continue
        mask_his = labels_his[sl] == his_id
        n_his_voxels = mask_his.sum()
        if n_his_voxels == 0:
            continue

        com_local = ndi.center_of_mass(mask_his)
        com_his_int = tuple(int(round(c + s.start)) for c, s in zip(com_local, sl))

        # Find which nucleus contains the histone centroid
        containing_nuc = int(labels_nuc[com_his_int])
        if containing_nuc == 0:
            continue

        # Calculate overlap fraction (histone voxels inside nucleus / total histone voxels)
        overlap_voxels = (mask_his & (labels_nuc[sl] == containing_nuc)).sum()
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
                        image_id, output_dir, voxel_size, sample_ri, immersion_ri,
                        seed_threshold=8.0, min_volume_um3=1.048,
                        seed_threshold_nuc=None, seed_threshold_his=None):
    """
    Generate QC figures for an image.

    Robust to very low nucleus counts (n=0, 1, or many): correlations,
    polyfits, and histogram bin counts are all guarded. When n<2 the
    correlation-related plots show a placeholder and pearson r/p is NaN.

    The fig2 threshold-sensitivity sweep scales with seed_threshold and the
    object-size filter uses a physical volume (min_volume_um3), so the figure
    stays meaningful for cameras/voxels other than the Zeiss widefield.
    """
    dz_c = voxel_size[0] * (sample_ri / immersion_ri)
    voxel_vol = dz_c * voxel_size[1] * voxel_size[2]
    # Size cutoff for the sensitivity sweep, in voxels of THIS image.
    # Uses the raw (uncorrected) voxel volume like segment_nuclei, so the
    # default reproduces the original 500-voxel cutoff at the CZI voxel.
    raw_voxel_vol = voxel_size[0] * voxel_size[1] * voxel_size[2]
    min_voxels_sens = int(round(min_volume_um3 / raw_voxel_vol))

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

    # Fig 2: Threshold sensitivity + distribution for both channels.
    # Each sweep is centered on the EFFECTIVE seed threshold of its own
    # channel (relevant when --seed-threshold auto gives the two channels
    # different values). At the default 8.0 this reproduces the original
    # [6, 8, 10, 12, 15, 18, 20] exactly.
    base_nuc = seed_threshold_nuc if seed_threshold_nuc is not None else seed_threshold
    base_his = seed_threshold_his if seed_threshold_his is not None else seed_threshold
    if isinstance(base_nuc, str):
        base_nuc = 8.0  # fallback only when no effective value was passed
    if isinstance(base_his, str):
        base_his = 8.0
    thresholds_nuc = [base_nuc * f for f in (0.75, 1.0, 1.25, 1.5, 1.875, 2.25, 2.5)]
    thresholds_his = [base_his * f for f in (0.75, 1.0, 1.25, 1.5, 1.875, 2.25, 2.5)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Nuclear channel
    sens_med_nuc, sens_n_nuc = [], []
    for thr in thresholds_nuc:
        binary = smoothed_nuc > thr
        labeled_t, n_t = ndi.label(binary)
        sizes_t = ndi.sum(binary, labeled_t, range(1, n_t + 1)).astype(int)
        vols_t = sizes_t[sizes_t > min_voxels_sens] * voxel_vol
        sens_med_nuc.append(np.median(vols_t) if len(vols_t) > 0 else 0)
        sens_n_nuc.append(len(vols_t))

    ax = axes[0, 0]
    ax.plot(thresholds_nuc, sens_med_nuc, 'o-', c='#0279EE', linewidth=2, markersize=8)
    ax2 = ax.twinx()
    ax2.plot(thresholds_nuc, sens_n_nuc, 's--', c='#FF9400', linewidth=1.5, markersize=6, alpha=0.7)
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
    for thr in thresholds_his:
        binary = smoothed_his > thr
        labeled_t, n_t = ndi.label(binary)
        sizes_t = ndi.sum(binary, labeled_t, range(1, n_t + 1)).astype(int)
        vols_t = sizes_t[sizes_t > min_voxels_sens] * voxel_vol
        sens_med_his.append(np.median(vols_t) if len(vols_t) > 0 else 0)
        sens_n_his.append(len(vols_t))

    ax = axes[1, 0]
    ax.plot(thresholds_his, sens_med_his, 'o-', c='#FD9BED', linewidth=2, markersize=8)
    ax2 = ax.twinx()
    ax2.plot(thresholds_his, sens_n_his, 's--', c='#FF9400', linewidth=1.5, markersize=6, alpha=0.7)
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
                  deconv_iterations=40, strain=None, modality='widefield',
                  ome_pair=None, image_id_override=None, seed_threshold=8.0,
                  seed_k=5.0, max_z_span_um=15.0, dz_um=None):
    """
    Process a single image through the full dual-channel pipeline.

    For .czi, pass filepath and the 1-based channel indices. For OME-TIFF,
    pass ome_pair=(file_x, file_y); filepath is then only used for logging.
    """
    if image_id_override:
        image_id = image_id_override
    elif ome_pair:
        image_id = split_ome_prefix(ome_pair[0])[0] or \
            os.path.splitext(os.path.basename(ome_pair[0]))[0]
    else:
        image_id = os.path.splitext(os.path.basename(filepath))[0]
    print(f"\n{'='*60}")
    print(f"Processing: {image_id}")
    print(f"{'='*60}")

    t_start = time.time()

    # Step 1: Read both channels
    if ome_pair is not None:
        channel_x, channel_y, voxel_size, na, immersion_ri, immersion = read_ome_tiff_pair(
            ome_pair[0], ome_pair[1])
    elif filepath.lower().endswith('.czi'):
        channel_x, channel_y, voxel_size, na, immersion_ri, immersion = read_czi_dual(
            filepath, channel_x_index, channel_y_index)
    else:
        raise ValueError(f"Unsupported input (need .czi or ome_pair): {filepath}")

    if dz_um is not None:
        voxel_size = (dz_um, voxel_size[1], voxel_size[2])
        print(f"  [override] dz = {dz_um} um (metadata PhysicalSizeZ ignored)")

    print(f"  Channel X (nuclear) shape: {channel_x.shape}, voxel: {voxel_size}, NA={na}, immersion_ri={immersion_ri}")
    print(f"  Channel Y (histone) shape: {channel_y.shape}")

    # Step 2: Background subtract both channels (in-place: *_bs IS the
    # same array as the raw channel; drop the extra references so the
    # later del actually frees memory)
    channel_x_bs, bg_x = background_subtract(channel_x)
    channel_y_bs, bg_y = background_subtract(channel_y)
    del channel_x, channel_y
    print(f"  Background (median): nuclear={bg_x:.1f}, histone={bg_y:.1f}")

    # Step 3: PSF for both channels
    psf_x = theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength_x,
                            modality=modality)
    psf_y = theoretical_psf(voxel_size, na, immersion_ri, emission_wavelength_y,
                            modality=modality)
    print(f"  PSF nuclear: {psf_x.shape}, histone: {psf_y.shape} ({modality})")

    # Step 4: Deconvolve both channels (skipped when iterations=0)
    if deconv_iterations > 0:
        print(f"  Deconvolving nuclear channel ({deconv_iterations} iterations)...")
        t0 = time.time()
        deconv_x = deconvolve(channel_x_bs, psf_x, iterations=deconv_iterations)
        del channel_x_bs  # free before the second RL run lowers its peak
        print(f"  Nuclear deconvolution: {time.time()-t0:.1f}s")

        print(f"  Deconvolving histone channel ({deconv_iterations} iterations)...")
        t0 = time.time()
        deconv_y = deconvolve(channel_y_bs, psf_y, iterations=deconv_iterations)
        del channel_y_bs
        print(f"  Histone deconvolution: {time.time()-t0:.1f}s")
    else:
        print("  Deconvolution skipped (--no-deconvolution)")
        deconv_x = channel_x_bs
        deconv_y = channel_y_bs

    # Steps 5-9: Segment both channels
    print(f"  Segmenting nuclear channel...")
    labels_nuc, smoothed_nuc, n_nuc, qc_nuc = segment_nuclei(
        deconv_x, voxel_size=voxel_size,
        seed_threshold=seed_threshold, seed_k=seed_k, max_z_span_um=max_z_span_um)
    del deconv_x  # only smoothed_nuc/labels_nuc are needed from here on
    n_removed_border_nuc = qc_nuc.pop('_n_removed_border', 0)
    seed_thr_nuc = qc_nuc.pop('_seed_threshold_effective', None)
    if n_removed_border_nuc > 0:
        print(f"  Removed {n_removed_border_nuc} nuclear objects touching image border (kept {n_nuc})")
    print(f"  Found {n_nuc} nuclear objects")

    print(f"  Segmenting histone channel...")
    labels_his, smoothed_his, n_his, qc_his = segment_nuclei(
        deconv_y, voxel_size=voxel_size,
        seed_threshold=seed_threshold, seed_k=seed_k, max_z_span_um=max_z_span_um)
    del deconv_y
    n_removed_border_his = qc_his.pop('_n_removed_border', 0)
    seed_thr_his = qc_his.pop('_seed_threshold_effective', None)
    if n_removed_border_his > 0:
        print(f"  Removed {n_removed_border_his} histone objects touching image border (kept {n_his})")
    print(f"  Found {n_his} histone objects")

    if isinstance(seed_threshold, str):
        print(f"  Adaptive seed ('{seed_threshold}', k={seed_k}): "
              f"nuclear={seed_thr_nuc:.2f}, histone={seed_thr_his:.2f}")

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
                                    image_id, output_dir, voxel_size, sample_ri, immersion_ri,
                                    seed_threshold=seed_threshold,
                                    seed_threshold_nuc=seed_thr_nuc,
                                    seed_threshold_his=seed_thr_his)

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
        'deconvolution': (f'RL {deconv_iterations} iter, Gaussian PSF ({modality})'
                          if deconv_iterations > 0 else 'none'),
    }

    elapsed = time.time() - t_start
    print(f"  Done in {elapsed:.1f}s")
    if len(ratios) > 0:
        print(f"  Ratio: median={np.median(ratios):.1f}%, "
              f"nuclear_vol={np.median(vols_nuc):.3f} µm³, "
              f"histone_vol={np.median(vols_his):.3f} µm³")
    else:
        print(f"  WARNING: 0 matched pairs (no overlapping objects found)")

    # Release per-image stacks (smoothed_*, labels_*) before the next image
    # so peak RSS does not grow across a batch.
    gc.collect()

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

    - If cli_value is not None it was already parsed as int by main().
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
                        help='Input .czi file, OME-TIFF file, or folder of image files')
    parser.add_argument('--output', required=True,
                        help='Output directory for results')
    parser.add_argument('--channel-x', type=str, default=None,
                        help='Nuclear marker channel. .czi: 1=mRFP1.2, 2=EGFP. '
                             'OME-TIFF: fluorophore tag substring from the filename, e.g. GFP or RFP. '
                             'If omitted for .czi, script asks interactively.')
    parser.add_argument('--channel-y', type=str, default=None,
                        help='Histone marker channel. Same rules as --channel-x.')
    parser.add_argument('--sample-ri', type=float, default=1.35,
                        help='Sample refractive index (default: 1.35 yeast)')
    parser.add_argument('--deconv-iterations', type=int, default=40,
                        help='RL deconvolution iterations (default: 40)')
    parser.add_argument('--strain-map',
                        help='CSV mapping filenames to strain labels (columns: filename, strain). '
                             'For OME-TIFF, the shared prefix (e.g. W21653_1) also works.')
    parser.add_argument('--no-deconvolution', action='store_true',
                        help='Skip deconvolution (faster, less accurate)')
    parser.add_argument('--modality', choices=['widefield', 'confocal'], default=None,
                        help='PSF model for deconvolution. Default: auto '
                             '(widefield for .czi, confocal for spinning-disk OME-TIFF).')
    parser.add_argument('--emission-x-nm', type=float, default=None,
                        help='Override emission wavelength (nm) of the nuclear channel.')
    parser.add_argument('--emission-y-nm', type=float, default=None,
                        help='Override emission wavelength (nm) of the histone channel.')
    parser.add_argument('--seed-threshold', default=8.0, metavar='FLOAT|auto',
                        help="Seed detection threshold on the smoothed stack: an absolute "
                             "intensity value (default: 8.0, calibrated for the Zeiss "
                             "widefield data) or 'auto' for a noise-adaptive threshold "
                             "(median + k*sigma, robust one-sided MAD per channel and "
                             "image; see --seed-k).")
    parser.add_argument('--seed-k', type=float, default=5.0,
                        help="Factor k for --seed-threshold auto: threshold = median + "
                             "k * sigma_robust (default: 5.0). Ignored otherwise.")
    parser.add_argument('--max-z-span-um', type=float, default=15.0,
                        help='Maximum z-extent of a seed object in um (default: 15.0). '
                             '15 um keeps an artifact guard for the Visitron stacks at '
                             'dz=0.2 um (longest observed real seed: 10.4 um); use 20+ '
                             'to effectively disable the filter on 81-slice stacks.')
    parser.add_argument('--dz-um', type=float, default=None,
                        help='Override the z voxel size in um when the metadata value '
                             'is wrong (e.g. VisiView writes PhysicalSizeZ=1 for a '
                             '0.2 um z-step). xy is still taken from the metadata.')
    args = parser.parse_args()

    # --seed-threshold accepts 'auto' or a number
    if isinstance(args.seed_threshold, str):
        if args.seed_threshold.lower() == 'auto':
            args.seed_threshold = 'auto'
        else:
            try:
                args.seed_threshold = float(args.seed_threshold)
            except ValueError:
                raise SystemExit(
                    "Error: --seed-threshold must be a number or 'auto', "
                    f"got '{args.seed_threshold}'.")

    os.makedirs(args.output, exist_ok=True)

    # Build file list
    if os.path.isdir(args.input):
        files = sorted([os.path.join(args.input, f) for f in os.listdir(args.input)
                       if f.lower().endswith(('.czi',) + OME_TIFF_EXTS)])
    else:
        files = [args.input]

    if not files:
        print(f"No .czi or OME-TIFF files found in {args.input}")
        sys.exit(1)

    czi_files = [f for f in files if f.lower().endswith('.czi')]
    ome_files = [f for f in files if not f.lower().endswith('.czi')]

    if czi_files and ome_files:
        raise SystemExit(
            "Error: mixed .czi and OME-TIFF inputs in one run are not supported.\n"
            "Please run the two formats in separate pipeline calls."
        )

    # ------------------------------------------------------------------
    # Mode-specific setup: channels, emission wavelengths, PSF modality,
    # and the job list (one job per image = one file for .czi, one
    # (nuclear file, histone file) pair for OME-TIFF).
    # ------------------------------------------------------------------
    jobs = []  # each: {'filepath' (czi) or 'file_x'+'file_y'+'prefix' (ome)}

    if ome_files:
        # ---- OME-TIFF mode: channels are filename tags ----
        if args.channel_x is None or args.channel_y is None:
            tags = sorted({split_ome_prefix(f)[1] for f in ome_files
                           if split_ome_prefix(f)[1] is not None})
            raise SystemExit(
                "Error: --channel-x and --channel-y are required for OME-TIFF input.\n"
                f"Available channel tags in this input: {tags}\n"
                "Example: --channel-x GFP --channel-y RFP"
            )
        tag_x, tag_y = args.channel_x, args.channel_y
        if tag_x.lower() == tag_y.lower():
            raise SystemExit(
                "Error: --channel-x and --channel-y must be different channels.\n"
                f"Both were set to '{tag_x}'."
            )

        def _emission_for_tag(tag):
            for name, nm in FLUOROPHORE_EMISSION_NM.items():
                if name.lower() in tag.lower():
                    return nm
            return None

        emission_wavelength_x = args.emission_x_nm or _emission_for_tag(tag_x)
        emission_wavelength_y = args.emission_y_nm or _emission_for_tag(tag_y)
        if emission_wavelength_x is None or emission_wavelength_y is None:
            raise SystemExit(
                "Error: could not derive emission wavelength from the channel tag.\n"
                f"Known fluorophores: {sorted(FLUOROPHORE_EMISSION_NM)}.\n"
                "Use --emission-x-nm / --emission-y-nm to set them explicitly."
            )

        modality = args.modality or 'confocal'  # OME-TIFF here = spinning disk
        print(f"Nuclear channel tag: '{tag_x}' (emission {emission_wavelength_x:.0f} nm)")
        print(f"Histone channel tag: '{tag_y}' (emission {emission_wavelength_y:.0f} nm)")
        print(f"PSF modality: {modality}")

        # Group files by shared prefix and pair the two channels
        groups = {}
        for f in ome_files:
            prefix, tag = split_ome_prefix(f)
            if prefix is None:
                print(f"  [WARN] Skipping file with unrecognized name: {os.path.basename(f)}")
                continue
            if 'dic' in tag.lower():
                continue  # transmitted-light channel, not used
            groups.setdefault(prefix, {})[tag] = f

        for prefix in sorted(groups):
            chans = groups[prefix]
            fx = next((f for t, f in chans.items() if tag_x.lower() in t.lower()), None)
            fy = next((f for t, f in chans.items() if tag_y.lower() in t.lower()), None)
            if fx is None or fy is None:
                print(f"  [WARN] Skipping '{prefix}': channel tag(s) '{tag_x}'/'{tag_y}' "
                      f"not found (available: {sorted(chans)})")
                continue
            jobs.append({'prefix': prefix, 'file_x': fx, 'file_y': fy})

        if not jobs:
            print("No complete channel pairs found. Check --channel-x/--channel-y.")
            sys.exit(1)
        print(f"Found {len(jobs)} image pair(s) to process.")

    else:
        # ---- CZI mode: channels are 1-based indices ----
        try:
            cx = int(args.channel_x) if args.channel_x is not None else None
            cy = int(args.channel_y) if args.channel_y is not None else None
        except ValueError:
            raise SystemExit(
                "Error: for .czi input, --channel-x/--channel-y must be channel "
                "numbers (1 or 2).\nFluorophore tags like 'GFP' only apply to "
                "OME-TIFF input."
            )
        args.channel_x = resolve_channel(cx, 'x')
        args.channel_y = resolve_channel(cy, 'y')

        if args.channel_x == args.channel_y:
            raise SystemExit(
                "Error: --channel-x and --channel-y must be different channels.\n"
                f"Both were set to channel {args.channel_x}."
            )

        # Emission wavelengths follow automatically from the channel choice
        # (607 nm for mRFP1.2 = channel 1, 509 nm for EGFP = channel 2).
        emission_wavelength_x = args.emission_x_nm or CHANNEL_EMISSION_NM[args.channel_x]
        emission_wavelength_y = args.emission_y_nm or CHANNEL_EMISSION_NM[args.channel_y]
        modality = args.modality or 'widefield'

        print(f"Nuclear channel: {args.channel_x} ({CHANNEL_NAMES[args.channel_x]}, "
              f"emission {emission_wavelength_x:.0f} nm)")
        print(f"Histone channel: {args.channel_y} ({CHANNEL_NAMES[args.channel_y]}, "
              f"emission {emission_wavelength_y:.0f} nm)")
        print(f"PSF modality: {modality}")

        for f in czi_files:
            jobs.append({'filepath': f})

    # --no-deconvolution: the flag now actually takes effect by setting
    # iterations to 0 (previously it was parsed but ignored).
    deconv_iterations = 0 if args.no_deconvolution else args.deconv_iterations

    # Load strain mapping (robust: strip whitespace, register extension-less
    # alias; lookup tries full filename, then extension-less name, then the
    # OME prefix).
    strain_map = {}
    if args.strain_map:
        with open(args.strain_map, newline='') as f:
            # German Excel saves "CSV" with semicolons -- sniff the delimiter
            # from the header so both comma- and semicolon-separated files work.
            header = f.readline()
            f.seek(0)
            try:
                dialect = csv.Sniffer().sniff(header, delimiters=',;\t')
            except csv.Error:
                dialect = 'excel'  # comma default
            reader = csv.DictReader(f, dialect=dialect)
            fields = {h.strip() for h in (reader.fieldnames or [])}
            if not {'filename', 'strain'} <= fields:
                raise SystemExit(
                    "Error: the strain map needs the columns 'filename' and "
                    f"'strain' (found: {reader.fieldnames}).")
            for row in reader:
                fname_key = row['filename'].strip()
                strain_val = row['strain'].strip()
                strain_map[fname_key] = strain_val
                # extension-less alias so both 'img.czi' and 'img' work
                strain_map.setdefault(os.path.splitext(fname_key)[0], strain_val)
        print(f"Loaded strain map: {len(strain_map)} entries from {args.strain_map}")
    else:
        print("NOTE: no --strain-map given; all pairs get strain column 'unassigned'.")

    def lookup_strain(*keys):
        for k in keys:
            if k is None:
                continue
            k = k.strip()
            if k in strain_map:
                return strain_map[k]
        return None

    # Failed-image log: initialize with header. Batch continues past crashes;
    # each failure is logged with error type, message, and timestamp so the
    # user can inspect + rerun those specific images later.
    failed_path = f'{args.output}/failed_images.tsv'
    with open(failed_path, 'w') as f:
        f.write('filename\terror_type\terror_message\ttimestamp\n')

    # Process all jobs with per-image try/except so an overnight batch never
    # aborts on a single bad image.
    all_results = []
    all_qc = []
    n_success = 0
    n_failed = 0
    n_unassigned = 0
    for job in jobs:
        if 'prefix' in job:  # OME-TIFF pair
            fname = job['prefix']
            strain = lookup_strain(os.path.basename(job['file_x']),
                                   os.path.splitext(os.path.basename(job['file_x']))[0],
                                   job['prefix'])
        else:  # CZI
            fname = os.path.basename(job['filepath'])
            strain = lookup_strain(fname, os.path.splitext(fname)[0])
        if strain is None:
            n_unassigned += 1
            print(f"  [WARN] No strain assignment for '{fname}' -> column 'unassigned'.")
        try:
            if 'prefix' in job:
                results, qc = process_image(
                    job['file_x'], args.output,
                    channel_x_index=None, channel_y_index=None,
                    sample_ri=args.sample_ri,
                    emission_wavelength_x=emission_wavelength_x,
                    emission_wavelength_y=emission_wavelength_y,
                    deconv_iterations=deconv_iterations,
                    strain=strain, modality=modality,
                    ome_pair=(job['file_x'], job['file_y']),
                    image_id_override=job['prefix'],
                    seed_threshold=args.seed_threshold,
                    seed_k=args.seed_k,
                    max_z_span_um=args.max_z_span_um,
                    dz_um=args.dz_um,
                )
            else:
                results, qc = process_image(
                    job['filepath'], args.output,
                    channel_x_index=args.channel_x,
                    channel_y_index=args.channel_y,
                    sample_ri=args.sample_ri,
                    emission_wavelength_x=emission_wavelength_x,
                    emission_wavelength_y=emission_wavelength_y,
                    deconv_iterations=deconv_iterations,
                    strain=strain, modality=modality,
                    seed_threshold=args.seed_threshold,
                    seed_k=args.seed_k,
                    max_z_span_um=args.max_z_span_um,
                    dz_um=args.dz_um,
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
    print(f"\nPipeline complete. {n_success} of {len(jobs)} images successful, "
          f"{n_failed} failed, {len(all_results)} matched pairs total.")
    if n_unassigned > 0:
        print(f"Note: {n_unassigned} image(s) had no strain assignment "
              f"(column 'unassigned'). Check --strain-map.")
    if n_failed > 0:
        print(f"See {failed_path} for details of the {n_failed} failed image(s).")


if __name__ == '__main__':
    main()
