# BarDiameter2 MVP (Rhino 8)

Minimal Rhino 8 Python script to estimate rebar diameter from noisy merged STL scans.

## What it does
- Select mesh once (it is automatically unselected right after pick so shading stays readable).
- Repeatedly pick as many points as needed on the same bar path; press Enter to finish each bar path.
- Script computes a shortest path on mesh between points (or straight fallback if disconnected).
- Every N mm (default 50 mm), it slices orthogonal to local path direction.
- For each valid slice, it computes area and equivalent diameter:
  - `d = sqrt(4A/pi)`
- Returns per measurement:
  - per-slice diameters
  - mean, median, std dev
  - 95% CI of mean

## New workflow improvements
- **Continuous measuring loop**: after naming a bar, it immediately starts next measurement. Press **Esc** at the first point prompt to finish session.
- **Bar naming + auto numbering**: enter bar ID like `405A`, results are saved as `405A_1`, `405A_2`, etc.
- **CSV export**: after finishing, save all collected measurements to CSV.
- **Visual confirmation geometry per bar**:
  - chosen slice polylines are added
  - sample points are added
  - each measurement gets its own sublayer under the active layer, named after saved bar ID (`TopBar_1`, `TopBar_2`, etc.).

## Run in Rhino 8
1. Open Rhino 8.
2. Run `EditPythonScript`.
3. Open `rhino8_bar_diameter_mvp.py` and run.
4. Follow prompts:
   - Select mesh
   - Slice spacing (mm)
   - For each bar: pick points along bar path (as many as needed), press Enter, then name bar
   - Press Esc at first point prompt when done with all bars
   - Export CSV at end

## Notes / limitations
- If bars are fused, pick points on a segment where one bar dominates.
- Works best when points are 100–500 mm apart (up to ~1000 mm supported).
- Loop-center distance is automatic (derived from spacing).
- If topology between picks is disconnected, script uses straight fallback between picks.
