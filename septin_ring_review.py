#!/usr/bin/env python3
"""
Interactive septin ring selection tool (companion to septin_ring_pipeline.py).

After a pipeline run, start this tool on the results folder:

    python septin_ring_review.py --results ./results

A gallery grid opens (24 rings per page, same panels as the QC gallery):
  - LEFT-CLICK a panel  : toggle keep (green frame = kept)
  - SCROLL over a panel : zoom in/out, centered on the cursor
  - 'r'                 : reset zoom on the current page
  - LEFT/RIGHT arrows   : previous/next page (buttons work too)
  - 's' / Save button   : write ring_diameter_per_strain_final.xlsx
  - closing the window  : also writes the final table

Nothing is selected at start -- click every ring you want to keep.
The selection is stored in selection_state.json after every click, so the
review can be interrupted and resumed at any time.

Inputs (all in the results folder):
  all_rings.xlsx           written by the pipeline (measurements)
  crops_<image>/ring_*.png per-ring panels written by the pipeline
                           (a placeholder panel is shown for rings without
                           a crop file, e.g. from older runs)

Output:
  ring_diameter_per_strain_final.xlsx
      One column per strain, the kept rings' primary diameters
      (xy_diameter_um, object width) in measurement order below.
      Shorter columns stay empty at the bottom.

NOTE: standalone script -- no imports from other pipeline files.
Requires a desktop Python with a GUI backend (TkAgg is the default on
Windows/macOS standard Python installs).
"""

import argparse
import json
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.widgets import Button
from matplotlib.image import imread

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

def _pick_font_family():
    """Use Jost when it is actually installed, else DejaVu Sans.

    Requesting a missing family makes matplotlib log a 'findfont' warning
    for EVERY text element -- checking once up front keeps the console
    clean. To use Jost: install the static TTFs and delete the matplotlib
    font cache (%USERPROFILE%/.matplotlib/fontlist-*.json)."""
    from matplotlib import font_manager
    try:
        font_manager.findfont('Jost', fallback_to_default=False)
        return ['Jost', 'DejaVu Sans']
    except ValueError:
        return ['DejaVu Sans']


rcParams['font.family'] = _pick_font_family()
rcParams['svg.fonttype'] = 'none'

REVIEW_VERSION = '1.4.0'

GRID_COLS = 6
GRID_ROWS = 4
PANELS_PER_PAGE = GRID_COLS * GRID_ROWS

STATE_FILENAME = 'selection_state.json'
FINAL_TABLE_FILENAME = 'ring_diameter_per_strain_final.xlsx'


# ============================================================
# Excel output helper (same formatting as the pipeline)
# ============================================================

FLOAT_FORMAT_DEFAULT = '0.0000'
HEADER_FONT = Font(bold=True)
HEADER_FILL = PatternFill(start_color='D9D9D9', end_color='D9D9D9',
                          fill_type='solid')
HEADER_ALIGN = Alignment(horizontal='center', vertical='center')


def save_xlsx(records, path, sheet_name):
    """Save a list of dict records as a formatted .xlsx file.

    - Bold header row with light-grey background, centered
    - Frozen top row
    - Column widths auto-fit to content (capped at 50 chars)
    - Float columns get 4 decimals

    Empty records: writes an empty file with no header, so downstream tools
    still see a file exists.
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
                f"Please close it and re-run the review tool."
            )
        return

    keys = list(records[0].keys())

    for col_idx, key in enumerate(keys, start=1):
        cell = ws.cell(row=1, column=col_idx, value=key)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = HEADER_ALIGN

    for row_idx, record in enumerate(records, start=2):
        for col_idx, key in enumerate(keys, start=1):
            value = record.get(key)
            if isinstance(value, (bool, np.bool_)):
                cell = ws.cell(row=row_idx, column=col_idx,
                               value=str(bool(value)))
            elif isinstance(value, (np.integer,)):
                cell = ws.cell(row=row_idx, column=col_idx, value=int(value))
            elif isinstance(value, (np.floating,)):
                cell = ws.cell(row=row_idx, column=col_idx, value=float(value))
            else:
                cell = ws.cell(row=row_idx, column=col_idx, value=value)

            if isinstance(value, (bool, np.bool_)):
                pass
            elif isinstance(value, (int, float, np.integer, np.floating)):
                cell.number_format = ('0' if isinstance(value, (int, np.integer))
                                      else FLOAT_FORMAT_DEFAULT)

    ws.freeze_panes = 'A2'

    for col_idx, key in enumerate(keys, start=1):
        max_len = len(str(key))
        for record in records:
            value = record.get(key)
            if value is None:
                continue
            if isinstance(value, (bool, np.bool_)):
                s = str(bool(value))
            elif isinstance(value, (float, np.floating)):
                s = f"{float(value):.4f}"
            else:
                s = str(value)
            if len(s) > max_len:
                max_len = len(s)
        width = min(max(max_len + 2, 8), 50)
        ws.column_dimensions[get_column_letter(col_idx)].width = width

    try:
        wb.save(path)
    except PermissionError:
        raise PermissionError(
            f"Cannot write '{path}'. File is likely open in Excel. "
            f"Please close it and re-run the review tool."
        )


# ============================================================
# Safe-text mode (FreeType raster overflow guard)
# ============================================================
#
# The interactive GUI rasterizes text through FreeType (unlike the
# pipeline, which saves SVG and converts via cairosvg). Certain FreeType
# builds crash with 'FT_Render_Glyph ... raster overflow' (error 0x62).
# The trigger is NOT the font family -- it is specific glyphs (em dash,
# micro sign, ...) and/or the antialiased smooth rasterizer, so the
# crash happens with DejaVu Sans too. The tool therefore escalates
# through increasingly safe text modes until rendering works:
#   level 0: normal (configured font, antialiased, Unicode)
#   level 1: antialiasing off  -> FreeType MONO rasterizer
#   level 2: + ASCII-only text (em dash -> '-', 'µ' -> 'u', ...)
#   level 3: + DejaVu Sans     (in case the configured font is the trigger)
# The real fix is updating FreeType:  conda update freetype matplotlib

_SAFE_LEVEL = 0
_SAFE_LEVEL_MAX = 3

_ASCII_REPL = {
    '—': '-', '–': '-', '−': '-', 'µ': 'u', '✓': 'OK',
    'ä': 'ae', 'ö': 'oe', 'ü': 'ue', 'Ä': 'Ae', 'Ö': 'Oe', 'Ü': 'Ue',
    'ß': 'ss',
}


def _is_raster_overflow(exc):
    return 'raster overflow' in str(exc) or 'FT_Render' in str(exc)


def _safe_text(s):
    """ASCII transliteration, active from safe-text level 2 on."""
    if _SAFE_LEVEL < 2:
        return s
    for k, v in _ASCII_REPL.items():
        s = s.replace(k, v)
    return s.encode('ascii', 'replace').decode('ascii')


def _text_kw():
    """Extra kwargs for text creation at the current safe level."""
    return {'antialiased': False} if _SAFE_LEVEL >= 1 else {}


def _escalate_safe_text():
    """Advance one safe-text level. Returns False when already at max."""
    global _SAFE_LEVEL
    if _SAFE_LEVEL >= _SAFE_LEVEL_MAX:
        return False
    _SAFE_LEVEL += 1
    if _SAFE_LEVEL >= 3:
        rcParams['font.family'] = ['DejaVu Sans']
    what = {1: 'Text-Antialiasing deaktiviert',
            2: 'nur ASCII-Zeichen',
            3: "Schriftart 'DejaVu Sans'"}[_SAFE_LEVEL]
    print(f"  [warn] FreeType raster overflow — sicherer Text-Modus "
          f"{_SAFE_LEVEL}/{_SAFE_LEVEL_MAX}: {what}.\n"
          f"         Tipp: 'conda update freetype matplotlib' behebt "
          f"die Ursache dauerhaft.")
    return True


def probe_font_or_fallback():
    """Probe-render the GUI's special characters at several sizes and
    escalate the safe-text mode until rendering works."""
    while True:
        try:
            fig = plt.figure(figsize=(3, 2))
            for size in (7, 8, 11, 14):
                fig.text(0.1, 0.2 * size / 14,
                         _safe_text('probe äöü — ✓ µm 0123'),
                         fontsize=size, **_text_kw())
            fig.canvas.draw()
            plt.close(fig)
            return
        except RuntimeError as e:
            plt.close('all')
            if not _is_raster_overflow(e):
                raise
            if not _escalate_safe_text():
                raise SystemExit(
                    "\nABBRUCH: matplotlib kann auf diesem System keinen "
                    "Text rastern — der FreeType raster overflow tritt bei "
                    "JEDER Glyphe auf (jede Schriftart, jeder Render-Modus, "
                    "selbst reines ASCII mit der mitgelieferten DejaVu "
                    "Sans).\n"
                    "Das ist eine defekte matplotlib/FreeType-Installation "
                    "(z.B. DLL-Konflikt), kein Problem dieses Skripts.\n\n"
                    "Reparatur im Anaconda Prompt:\n"
                    "  1) conda install --force-reinstall -c conda-forge "
                    "matplotlib freetype\n"
                    "  2) falls das nicht hilft: python -m pip install "
                    "--force-reinstall --no-cache-dir matplotlib\n"
                    "     (das pip-Paket bringt sein eigenes FreeType mit)\n"
                    "Danach dieses Tool erneut starten.")


# ============================================================
# Selection model (GUI-independent, headless-testable)
# ============================================================

class RingSelection:
    """Loads all_rings.xlsx, manages the keep set, writes the final table."""

    def __init__(self, results_dir):
        self.results_dir = results_dir
        self.records = self._load_records(
            os.path.join(results_dir, 'all_rings.xlsx'))
        self.kept = set()  # keys "image_id|ring_id"
        self.load_state()

    @staticmethod
    def _load_records(path):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"all_rings.xlsx not found in '{os.path.dirname(path)}'. "
                f"Run septin_ring_pipeline.py first.")
        wb = openpyxl.load_workbook(path, read_only=True)
        ws = wb.active
        rows = list(ws.iter_rows(values_only=True))
        wb.close()
        if not rows:
            return []
        keys = [str(k) for k in rows[0]]
        records = []
        for row in rows[1:]:
            rec = {k: v for k, v in zip(keys, row)}
            if rec.get('image_id') is None or rec.get('ring_id') is None:
                continue
            records.append(rec)
        return records

    @staticmethod
    def key(rec):
        return f"{rec['image_id']}|{rec['ring_id']}"

    def __len__(self):
        return len(self.records)

    def is_kept(self, idx):
        return self.key(self.records[idx]) in self.kept

    def toggle(self, idx):
        k = self.key(self.records[idx])
        if k in self.kept:
            self.kept.discard(k)
        else:
            self.kept.add(k)
        self.save_state()

    @property
    def n_kept(self):
        return len(self.kept)

    def crop_path(self, rec):
        return os.path.join(self.results_dir, f"crops_{rec['image_id']}",
                            f"ring_{int(rec['ring_id']):03d}.png")

    # --- state persistence ---
    @property
    def state_path(self):
        return os.path.join(self.results_dir, STATE_FILENAME)

    def save_state(self):
        with open(self.state_path, 'w') as f:
            json.dump({'kept': sorted(self.kept)}, f)

    def load_state(self):
        if not os.path.exists(self.state_path):
            return
        try:
            with open(self.state_path) as f:
                state = json.load(f)
            valid = {self.key(r) for r in self.records}
            self.kept = {k for k in state.get('kept', []) if k in valid}
        except (json.JSONDecodeError, OSError) as e:
            print(f"  [warn] could not read {STATE_FILENAME} ({e}), "
                  f"starting with empty selection")

    # --- final table ---
    def write_final_table(self):
        """ring_diameter_per_strain_final.xlsx: one column per strain,
        kept rings' diameters in measurement order, shorter columns
        padded empty. All strains appear as columns even if a strain
        currently has no kept rings."""
        strains = sorted({str(r.get('strain') or 'unassigned')
                          for r in self.records})
        kept_recs = [r for r in self.records if self.key(r) in self.kept]
        diam_by_strain = {s: [] for s in strains}
        for r in kept_recs:
            d = r.get('xy_diameter_um')
            if d is None:
                continue
            d = float(d)
            if np.isfinite(d):
                diam_by_strain[str(r.get('strain') or 'unassigned')].append(d)
        n_max = max((len(v) for v in diam_by_strain.values()), default=0)
        rows = [{s: (diam_by_strain[s][i] if i < len(diam_by_strain[s])
                     else None)
                 for s in strains}
                for i in range(n_max)]
        path = os.path.join(self.results_dir, FINAL_TABLE_FILENAME)
        save_xlsx(rows, path, sheet_name='ring_diameter_per_strain')
        return path, len(kept_recs)


# ============================================================
# GUI
# ============================================================

class ReviewGUI:
    """Gallery grid with click-to-keep selection."""

    def __init__(self, model):
        self.model = model
        self.page = 0
        self.n_pages = max(1, int(np.ceil(len(model) / PANELS_PER_PAGE)))

        self.fig = plt.figure(figsize=(18, 9))
        self.axes = []
        for i in range(PANELS_PER_PAGE):
            ax = self.fig.add_subplot(GRID_ROWS, GRID_COLS, i + 1)
            ax.set_xticks([]); ax.set_yticks([])
            self.axes.append(ax)

        # navigation / save buttons at the bottom
        ax_prev = self.fig.add_axes([0.30, 0.005, 0.10, 0.035])
        ax_save = self.fig.add_axes([0.45, 0.005, 0.10, 0.035])
        ax_next = self.fig.add_axes([0.60, 0.005, 0.10, 0.035])
        self.btn_prev = Button(ax_prev, '< Prev')
        self.btn_save = Button(ax_save, 'Save (s)')
        self.btn_next = Button(ax_next, 'Next >')
        self.btn_prev.on_clicked(lambda _e: self.prev_page())
        self.btn_next.on_clicked(lambda _e: self.next_page())
        self.btn_save.on_clicked(lambda _e: self.save())

        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.95,
                                 bottom=0.05, hspace=0.35, wspace=0.05)
        self.render_page()

    # --- rendering ---
    def render_page(self):
        start = self.page * PANELS_PER_PAGE
        for pi, ax in enumerate(self.axes):
            ax.clear()
            ax.set_xticks([]); ax.set_yticks([])
            idx = start + pi
            if idx >= len(self.model):
                ax.axis('off')
                continue
            rec = self.model.records[idx]
            crop_file = self.model.crop_path(rec)
            if os.path.exists(crop_file):
                ax.imshow(imread(crop_file))
            else:
                # placeholder for rings without a saved crop (old runs)
                ax.text(0.5, 0.5,
                        _safe_text(f"Ring {rec['ring_id']}\n"
                                   f"d={self._fmt_diam(rec)}\n"
                                   f"[{rec.get('qc_flags', '')}]"),
                        ha='center', va='center', fontsize=9,
                        transform=ax.transAxes, **_text_kw())
                ax.set_facecolor('0.95')
            kept = self.model.is_kept(idx)
            title = f"{rec['image_id']}  #{rec['ring_id']}"
            if kept:
                title = u'\u2713 ' + title  # checkmark
            ax.set_title(_safe_text(title), fontsize=8,
                         color='green' if kept else 'black', **_text_kw())
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color('green' if kept else '0.75')
                spine.set_linewidth(3 if kept else 1)
        self.fig.suptitle(
            _safe_text(f"Page {self.page + 1}/{self.n_pages}   —   "
                       f"{self.model.n_kept} of {len(self.model)} rings kept"
                       f"   —   click = keep/unkeep, scroll = zoom, "
                       f"r = reset zoom, arrows = page, s = save"),
            fontsize=11, **_text_kw())
        # covers artists created without _text_kw (e.g. button labels)
        self._apply_safe_text_everywhere()
        while True:
            try:
                self.fig.canvas.draw_idle()
                break
            except RuntimeError as e:
                # the startup probe may pass while a real render still
                # trips the overflow (other glyphs/sizes) -- escalate the
                # safe-text mode and update the existing artists
                if not _is_raster_overflow(e) or not _escalate_safe_text():
                    raise
                self._apply_safe_text_everywhere()

    def _apply_safe_text_everywhere(self):
        """Apply the current safe-text mode to every existing artist."""
        if _SAFE_LEVEL == 0:
            return
        texts = []
        if self.fig._suptitle is not None:
            texts.append(self.fig._suptitle)
        for ax in self.fig.axes:
            texts.append(ax.title)
            texts.extend(ax.texts)
        for t in texts:
            if _SAFE_LEVEL >= 1:
                t.set_antialiased(False)
            if _SAFE_LEVEL >= 2:
                t.set_text(_safe_text(t.get_text()))
            if _SAFE_LEVEL >= 3:
                t.set_fontfamily('DejaVu Sans')

    @staticmethod
    def _fmt_diam(rec):
        d = rec.get('xy_diameter_um')
        if d is None or not np.isfinite(float(d)):
            return 'n/a'
        return f"{float(d):.3f} um"

    # --- events ---
    ZOOM_STEP = 0.75   # limit scale per scroll notch when zooming in
    ZOOM_MIN_PX = 5    # stop zooming in at this many pixels across

    def on_click(self, event):
        if event.inaxes not in self.axes:
            return
        idx = self.page * PANELS_PER_PAGE + self.axes.index(event.inaxes)
        if idx >= len(self.model):
            return
        self.model.toggle(idx)
        self.render_page()

    def on_scroll(self, event):
        """Scroll-wheel zoom on the panel under the cursor, anchored at the
        cursor position. Limits stay clamped to the image; changing the
        page or pressing 'r' resets the zoom (render_page clears axes)."""
        ax = event.inaxes
        if ax not in self.axes or not ax.images:
            return
        h, w = ax.images[0].get_array().shape[:2]
        factor = self.ZOOM_STEP if event.button == 'up' else 1.0 / self.ZOOM_STEP
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()  # inverted (imshow origin='upper')
        x = event.xdata if event.xdata is not None else (x0 + x1) / 2
        y = event.ydata if event.ydata is not None else (y0 + y1) / 2
        nx = self._clamp_limits((x - (x - x0) * factor,
                                 x + (x1 - x) * factor), w)
        ny = self._clamp_limits((y - (y - y0) * factor,
                                 y + (y1 - y) * factor), h)
        if nx is None or ny is None:
            return  # max zoom reached
        ax.set_xlim(*nx)
        ax.set_ylim(ny[1], ny[0])  # keep y inverted
        self.fig.canvas.draw_idle()

    def _clamp_limits(self, lim, n_px):
        """Clamp a (lo, hi) limit pair to the image extent [-0.5, n_px-0.5],
        preserving width. Returns None when already at the zoom-in limit."""
        lo, hi = min(lim), max(lim)
        width = hi - lo
        if width < self.ZOOM_MIN_PX:
            return None
        if width >= n_px:
            return (-0.5, n_px - 0.5)
        if lo < -0.5:
            hi -= lo + 0.5
            lo = -0.5
        if hi > n_px - 0.5:
            lo -= hi - (n_px - 0.5)
            hi = n_px - 0.5
        return (lo, hi)

    def on_key(self, event):
        if event.key == 'right':
            self.next_page()
        elif event.key == 'left':
            self.prev_page()
        elif event.key == 's':
            self.save()
        elif event.key == 'r':
            self.render_page()  # axes are cleared -> zoom reset

    def on_close(self, _event):
        self.save()

    # --- actions ---
    def next_page(self):
        if self.page < self.n_pages - 1:
            self.page += 1
            self.render_page()

    def prev_page(self):
        if self.page > 0:
            self.page -= 1
            self.render_page()

    def save(self):
        path, n = self.model.write_final_table()
        print(f"Saved {n} kept rings -> {path}")

    def show(self):
        plt.show()


def build_gui_with_fallback(model):
    """Build the ReviewGUI; on FreeType raster overflow escalate the
    safe-text mode and rebuild. Text artists capture font settings at
    creation time, so the whole figure must be rebuilt after a mode
    change."""
    while True:
        try:
            return ReviewGUI(model)
        except RuntimeError as e:
            plt.close('all')
            if not _is_raster_overflow(e) or not _escalate_safe_text():
                raise
            print("  [warn] Fenster wird mit sicherem Text-Modus neu "
                  "aufgebaut.")


# ============================================================
# CLI
# ============================================================

_INTERACTIVE_BACKENDS = {'tkagg', 'qtagg', 'qt5agg', 'qt6agg', 'wxagg',
                         'macosx', 'gtk3agg', 'gtk4agg', 'webagg'}


def ensure_interactive_backend():
    """Make sure a window-capable backend is active.

    Started from Jupyter (!python / %run) matplotlib often runs on the
    non-interactive 'agg'/'inline' backend -- plt.show() then opens NO
    window. Try to switch to a GUI backend; if none is available, exit
    with instructions instead of silently showing nothing.
    """
    if matplotlib.get_backend().lower() in _INTERACTIVE_BACKENDS:
        return
    for candidate in ('TkAgg', 'QtAgg', 'Qt5Agg'):
        try:
            matplotlib.use(candidate)
        except Exception:
            continue
        if matplotlib.get_backend().lower() in _INTERACTIVE_BACKENDS:
            print(f"  Backend auf {candidate} umgeschaltet.")
            return
    raise SystemExit(
        f"\nABBRUCH: Kein Fenster-Backend verfügbar "
        f"(aktiv: '{matplotlib.get_backend()}').\n"
        "Das interaktive Fenster kann so nicht angezeigt werden.\n\n"
        "Bitte im ANACONDA PROMPT starten (nicht in Jupyter):\n"
        "  cd /d \"<Ordner mit septin_ring_review.py>\"\n"
        "  python septin_ring_review.py --results ./results\n\n"
        "Falls dann dasselbe kommt, fehlt tkinter:\n"
        "  conda install tk\n"
        "und erneut versuchen.")


def main():
    parser = argparse.ArgumentParser(
        description='Interactive septin ring selection: click the rings '
                    'to keep, then save ring_diameter_per_strain_final.xlsx.')
    parser.add_argument('--results', required=True,
                        help='Results folder from septin_ring_pipeline.py '
                             '(contains all_rings.xlsx and crops_* folders)')
    args = parser.parse_args()

    print(f"septin_ring_review.py v{REVIEW_VERSION}")
    try:
        import matplotlib.ft2font as _ft
        print(f"  matplotlib {matplotlib.__version__}, "
              f"FreeType {_ft.__freetype_version__}")
    except Exception:
        pass
    ensure_interactive_backend()
    print(f"  Backend: {matplotlib.get_backend()}")
    probe_font_or_fallback()

    model = RingSelection(args.results)
    if len(model) == 0:
        print(f"No rings in {args.results}/all_rings.xlsx — nothing to "
              f"review.")
        return
    print(f"{len(model)} rings loaded, {model.n_kept} already kept "
          f"(state: {model.state_path})")
    gui = build_gui_with_fallback(model)
    gui.show()


if __name__ == '__main__':
    main()
