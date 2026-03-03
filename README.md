# BarDiameter2 MVP (Rhino 8)

Minimal Rhino 8 Python script to estimate rebar diameter from noisy merged STL scans.

## What it does
- Select mesh once.
- Repeatedly pick 2 points on the same bar.
- Script computes a shortest path on mesh between points (or straight fallback if disconnected).
- Every N mm (default 50 mm), it slices orthogonal to local path direction.
- For each valid slice, it computes area and equivalent diameter:
  - `d = sqrt(4A/pi)`
- Returns per measurement:
  - per-slice diameters
  - mean, median, std dev
  - 95% CI of mean

## New workflow improvements
- **Continuous measuring loop**: after each measurement, it asks if you want to measure another bar.
- **Bar naming + auto numbering**: enter bar ID like `405A`, results are saved as `405A_1`, `405A_2`, etc.
- **CSV export**: after finishing, save all collected measurements to CSV.
- **Visual confirmation geometry**:
  - sample points are added
  - chosen slice polylines are added
  - both are added to `BarDiameterMVP_Slices` sublayer under current layer.

## Run in Rhino 8
1. Open Rhino 8.
2. Run `EditPythonScript`.
3. Open `rhino8_bar_diameter_mvp.py` and run.
4. Follow prompts:
   - Select mesh
   - Slice spacing (mm)
   - Repeatedly pick start/end points, name bar, continue or stop
   - Export CSV at end

## Notes / limitations
- If bars are fused, pick points on a segment where one bar dominates.
- Works best when points are 100–500 mm apart (up to ~1000 mm supported).
- Loop-center distance is automatic (derived from spacing).
- If topology between picks is disconnected, script uses straight fallback between picks.
