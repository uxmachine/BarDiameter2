# BarDiameter2 MVP (Rhino 8)

Minimal Rhino 8 Python script to estimate rebar diameter from noisy merged STL scans.

## What it does
- You pick one mesh.
- You pick 2 points on the **same bar**.
- Script computes a shortest path on mesh between those points.
- Every N mm (default 50 mm), it slices orthogonal to the local path direction.
- For each valid slice, it computes area and converts to equivalent diameter:
  - `d = sqrt(4A/pi)`
- Returns:
  - per-slice diameters
  - mean, median, std dev
  - 95% CI of mean

## Why this is MVP
- No install outside Rhino 8.
- One script (`rhino8_bar_diameter_mvp.py`).
- Manual point-picking stays (robust and simple for internal use).

## Run in Rhino 8
1. Open Rhino 8.
2. Run `EditPythonScript`.
3. Open `rhino8_bar_diameter_mvp.py` and run.
4. Follow prompts:
   - Select mesh
   - Pick start/end on same bar
   - Slice spacing (mm)
   - Max loop-center distance

## Notes / limitations
- If bars are fused, you must pick points on a segment where one bar dominates.
- Works best when points are 100–500 mm apart.
- Bent bars are handled better than straight-axis slicing because tangent comes from mesh path.
- If too few valid slices are found, increase `Max loop-center distance` or repick cleaner points.
