"""
Rhino 8 MVP script: estimate rebar diameter from a noisy merged STL mesh.

Usage in Rhino 8:
1) Open PythonScript editor.
2) Load this file and run `BarDiameterMVP()`.
3) Select mesh once, then keep measuring bars in a loop.
"""

import csv
import math
from datetime import datetime

import Rhino
import System
import rhinoscriptsyntax as rs
import scriptcontext as sc


def _dist2(a, b):
    dx = a.X - b.X
    dy = a.Y - b.Y
    dz = a.Z - b.Z
    return dx * dx + dy * dy + dz * dz


def _lerp_point(a, b, t):
    return Rhino.Geometry.Point3d(
        a.X + (b.X - a.X) * t,
        a.Y + (b.Y - a.Y) * t,
        a.Z + (b.Z - a.Z) * t,
    )


def _nearest_mesh_vertex_index(mesh, point3d):
    best_i = None
    best_d2 = float("inf")
    for i, v in enumerate(mesh.Vertices):
        d2 = _dist2(v, point3d)
        if d2 < best_d2:
            best_d2 = d2
            best_i = i
    return best_i


def _closest_topology_vertex_index(mesh, point3d):
    topo = mesh.TopologyVertices
    mp = mesh.ClosestMeshPoint(point3d, 0.0)

    if mp and hasattr(mp, "FaceIndex") and mp.FaceIndex is not None and mp.FaceIndex >= 0:
        fi = mp.FaceIndex
        if fi < mesh.Faces.Count:
            face = mesh.Faces[fi]
            mesh_vis = [face.A, face.B, face.C]
            if face.IsQuad:
                mesh_vis.append(face.D)
            best_topo = None
            best_d2 = float("inf")
            for mvi in mesh_vis:
                tvi = topo.TopologyVertexIndex(mvi)
                if tvi < 0:
                    continue
                d2 = _dist2(topo[tvi], point3d)
                if d2 < best_d2:
                    best_d2 = d2
                    best_topo = tvi
            if best_topo is not None:
                return best_topo

    near_mvi = _nearest_mesh_vertex_index(mesh, point3d)
    if near_mvi is None:
        return None
    tvi = topo.TopologyVertexIndex(near_mvi)
    return tvi if tvi >= 0 else None


def _build_mesh_graph(mesh):
    topo = mesh.TopologyVertices
    edge_count = mesh.TopologyEdges.Count
    graph = {i: set() for i in range(topo.Count)}
    for ei in range(edge_count):
        tv = mesh.TopologyEdges.GetTopologyVertices(ei)
        a = tv.I
        b = tv.J
        if a >= 0 and b >= 0:
            graph[a].add(b)
            graph[b].add(a)
    return graph


def _dijkstra_path(mesh, graph, start, end):
    import heapq

    topo = mesh.TopologyVertices
    dist = {start: 0.0}
    prev = {}
    visited = set()
    heap = [(0.0, start)]

    while heap:
        cur_d, u = heapq.heappop(heap)
        if u in visited:
            continue
        visited.add(u)
        if u == end:
            break

        pu = topo[u]
        for v in graph.get(u, []):
            if v in visited:
                continue
            pv = topo[v]
            w = pu.DistanceTo(pv)
            nd = cur_d + w
            if nd < dist.get(v, float("inf")):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(heap, (nd, v))

    if end not in dist:
        return None

    path = [end]
    cur = end
    while cur != start:
        cur = prev[cur]
        path.append(cur)
    path.reverse()
    return path


def _resample_polyline(points, spacing):
    if len(points) < 2:
        return points

    cum = [0.0]
    for i in range(1, len(points)):
        cum.append(cum[-1] + points[i - 1].DistanceTo(points[i]))

    total = cum[-1]
    if total <= spacing:
        return points

    targets = []
    d = 0.0
    while d <= total:
        targets.append(d)
        d += spacing
    if total - targets[-1] > 0.25 * spacing:
        targets.append(total)

    out = []
    seg = 0
    for t in targets:
        while seg < len(cum) - 2 and cum[seg + 1] < t:
            seg += 1
        a = points[seg]
        b = points[seg + 1]
        seg_len = cum[seg + 1] - cum[seg]
        if seg_len <= 1e-9:
            out.append(a)
            continue
        alpha = (t - cum[seg]) / seg_len
        out.append(_lerp_point(a, b, alpha))
    return out


def _local_tangent(points, i):
    n = len(points)
    if n < 2:
        return Rhino.Geometry.Vector3d.Unset
    if i == 0:
        v = points[1] - points[0]
    elif i == n - 1:
        v = points[-1] - points[-2]
    else:
        v = points[i + 1] - points[i - 1]
    if not v.Unitize():
        return Rhino.Geometry.Vector3d.Unset
    return v


def _line_samples(p1, p2, spacing):
    crv = Rhino.Geometry.LineCurve(p1, p2)
    seg_len = p1.DistanceTo(p2)
    if seg_len <= spacing:
        return [p1, p2]
    tvals = crv.DivideByLength(spacing, True)
    if not tvals:
        return [p1, p2]
    return [crv.PointAt(t) for t in tvals]


def _equivalent_diameter_from_area(area):
    if area <= 0:
        return None
    return math.sqrt((4.0 * area) / math.pi)


def _curve_area(curve):
    amp = Rhino.Geometry.AreaMassProperties.Compute(curve)
    if not amp:
        return None
    return abs(amp.Area)


def _auto_max_center_dist(spacing_mm):
    return max(12.0, min(35.0, 0.5 * spacing_mm))


def _pick_relevant_loop(section_curves, sample_point, max_center_dist):
    best_curve = None
    best_dist = float("inf")

    for c in section_curves:
        if not c.IsClosed:
            continue
        ok, _ = c.TryGetPlane()
        if not ok:
            continue
        amp = Rhino.Geometry.AreaMassProperties.Compute(c)
        if not amp:
            continue
        center = amp.Centroid
        d = center.DistanceTo(sample_point)
        if d < best_dist and d <= max_center_dist:
            best_dist = d
            best_curve = c

    return best_curve


def _ensure_child_layer(parent_index, child_name, color=None):
    parent_layer = sc.doc.Layers[parent_index]
    full = "{}::{}".format(parent_layer.FullPath, child_name)
    existing = sc.doc.Layers.FindByFullPath(full, -1)
    if existing >= 0:
        return existing

    layer = Rhino.DocObjects.Layer()
    layer.Name = child_name
    layer.ParentLayerId = parent_layer.Id
    if color is not None:
        layer.Color = color
    return sc.doc.Layers.Add(layer)


def _add_curve_on_layer(curve, layer_index):
    if layer_index is None or layer_index < 0:
        sc.doc.Objects.AddCurve(curve)
        return
    attr = Rhino.DocObjects.ObjectAttributes()
    attr.LayerIndex = layer_index
    sc.doc.Objects.AddCurve(curve, attr)


def _add_point_on_layer(pt, layer_index):
    if layer_index is None or layer_index < 0:
        sc.doc.Objects.AddPoint(pt)
        return
    attr = Rhino.DocObjects.ObjectAttributes()
    attr.LayerIndex = layer_index
    sc.doc.Objects.AddPoint(pt, attr)




def _samples_for_segment(mesh, graph, p_start, p_end, spacing):
    s_topo = _closest_topology_vertex_index(mesh, p_start)
    e_topo = _closest_topology_vertex_index(mesh, p_end)
    if s_topo is None or e_topo is None:
        return None

    topo = mesh.TopologyVertices
    path_ids = _dijkstra_path(mesh, graph, s_topo, e_topo)
    if path_ids and len(path_ids) >= 2:
        raw_pts = [topo[i] for i in path_ids]
        Rhino.RhinoApp.WriteLine("Path mode: mesh topology path used (segment).")
        return _resample_polyline(raw_pts, spacing)

    Rhino.RhinoApp.WriteLine("Path mode: straight fallback used (segment path not connected).")
    return _line_samples(p_start, p_end, spacing)


def _build_samples_from_waypoints(mesh, graph, waypoints, spacing):
    all_samples = []
    for i in range(len(waypoints) - 1):
        seg = _samples_for_segment(mesh, graph, waypoints[i], waypoints[i + 1], spacing)
        if not seg or len(seg) < 2:
            continue
        if all_samples:
            all_samples.extend(seg[1:])
        else:
            all_samples.extend(seg)
    return all_samples


def _measurement_stats(diameters):
    n = len(diameters)
    mean_d = sum(diameters) / n
    sorted_d = sorted(diameters)
    median_d = sorted_d[n // 2] if n % 2 == 1 else 0.5 * (sorted_d[n // 2 - 1] + sorted_d[n // 2])
    var = sum((x - mean_d) ** 2 for x in diameters) / max(1, n - 1)
    std = math.sqrt(var)
    half_ci = 1.96 * std / math.sqrt(n)
    return {
        "n": n,
        "mean": mean_d,
        "median": median_d,
        "std": std,
        "ci_low": mean_d - half_ci,
        "ci_high": mean_d + half_ci,
    }


def _measure_once(mesh_obj, mesh, graph, spacing):
    point_count = rs.GetInteger("Number of points for this bar path (2-10)", 3, 2, 10)
    if point_count is None:
        return None

    waypoints = []
    for i in range(point_count):
        if i == 0:
            msg = "Pick point 1/{} on bar (Esc to finish session)".format(point_count)
        else:
            msg = "Pick point {}/{} on SAME bar path".format(i + 1, point_count)
        p = rs.GetPointOnMesh(mesh_obj.ObjectId, msg)
        if not p:
            if i == 0:
                return None
            Rhino.RhinoApp.WriteLine("Point picking cancelled; measurement skipped.")
            return False
        waypoints.append(p)

    max_center_dist = _auto_max_center_dist(spacing)
    Rhino.RhinoApp.WriteLine("Using auto loop-center distance: {:.1f} mm".format(max_center_dist))

    samples = _build_samples_from_waypoints(mesh, graph, waypoints, spacing)
    if not samples or len(samples) < 2:
        Rhino.RhinoApp.WriteLine("Could not build path samples from selected points.")
        return False

    diameters = []
    sample_points = []
    sample_curves = []

    for i, pt in enumerate(samples):
        tan = _local_tangent(samples, i)
        if not tan.IsValid or tan.IsTiny():
            continue

        plane = Rhino.Geometry.Plane(pt, tan)
        polylines = Rhino.Geometry.Intersect.Intersection.MeshPlane(mesh, plane)
        if not polylines:
            continue

        section_curves = []
        for pl in polylines:
            if pl is None or pl.Count < 4:
                continue
            first_pt = pl[0]
            last_pt = pl[pl.Count - 1]
            if not first_pt.EpsilonEquals(last_pt, 1e-6):
                pl.Add(first_pt)
            c = Rhino.Geometry.PolylineCurve(pl)
            if c and c.IsClosed:
                section_curves.append(c)

        if not section_curves:
            continue

        chosen = _pick_relevant_loop(section_curves, pt, max_center_dist)
        if not chosen:
            continue

        area = _curve_area(chosen)
        if area is None or area <= 0:
            continue

        d = _equivalent_diameter_from_area(area)
        if d is None:
            continue

        diameters.append(d)
        sample_points.append(pt)
        sample_curves.append(chosen.DuplicateCurve())

    n = len(diameters)
    if n < 3:
        Rhino.RhinoApp.WriteLine("Not enough valid slices. Repick cleaner points or increase slice spacing.")
        return False

    stats = _measurement_stats(diameters)
    Rhino.RhinoApp.WriteLine("--- BarDiameterMVP ---")
    Rhino.RhinoApp.WriteLine("Valid slices: {}".format(stats["n"]))
    Rhino.RhinoApp.WriteLine("Diameters (mm): {}".format(", ".join("{:.2f}".format(x) for x in diameters)))
    Rhino.RhinoApp.WriteLine("Mean diameter (mm): {:.2f}".format(stats["mean"]))
    Rhino.RhinoApp.WriteLine("Median diameter (mm): {:.2f}".format(stats["median"]))
    Rhino.RhinoApp.WriteLine("Std dev (mm): {:.2f}".format(stats["std"]))
    Rhino.RhinoApp.WriteLine("95% CI for mean (mm): [{:.2f}, {:.2f}]".format(stats["ci_low"], stats["ci_high"]))

    sc.doc.Views.Redraw()
    stats["diameters"] = diameters
    stats["sample_points"] = sample_points
    stats["sample_curves"] = sample_curves
    return stats


def _unique_bar_name(base_name, name_counts):
    base = (base_name or "BAR").strip() or "BAR"
    count = name_counts.get(base, 0) + 1
    name_counts[base] = count
    return "{}_{}".format(base, count)


def _csv_safe(value):
    try:
        unicode_type = unicode  # type: ignore[name-defined]
    except NameError:
        unicode_type = str

    if isinstance(value, unicode_type):
        try:
            return value.encode("utf-8")
        except Exception:
            return str(value)
    return value


def _export_csv(rows):
    path = rs.SaveFileName(
        "Save bar diameter results as CSV",
        "CSV Files (*.csv)|*.csv||",
        None,
        "bar_diameter_results.csv",
        "csv",
    )
    if not path:
        Rhino.RhinoApp.WriteLine("CSV export skipped.")
        return

    headers = [
        "measurement_name",
        "input_bar_name",
        "timestamp",
        "slice_spacing_mm",
        "valid_slices",
        "mean_mm",
        "median_mm",
        "std_mm",
        "ci95_low_mm",
        "ci95_high_mm",
        "slice_diameters_mm",
    ]

    with open(path, "wb") as f:
        w = csv.writer(f)
        w.writerow([_csv_safe(h) for h in headers])
        for r in rows:
            w.writerow([
                _csv_safe(r["name"]),
                _csv_safe(r["base_name"]),
                _csv_safe(r["timestamp"]),
                _csv_safe("{:.2f}".format(r["spacing"])),
                _csv_safe(r["n"]),
                _csv_safe("{:.3f}".format(r["mean"])),
                _csv_safe("{:.3f}".format(r["median"])),
                _csv_safe("{:.3f}".format(r["std"])),
                _csv_safe("{:.3f}".format(r["ci_low"])),
                _csv_safe("{:.3f}".format(r["ci_high"])),
                _csv_safe(";".join("{:.3f}".format(x) for x in r["diameters"])),
            ])

    Rhino.RhinoApp.WriteLine("CSV exported: {}".format(path))


def BarDiameterMVP():
    go = Rhino.Input.Custom.GetObject()
    go.SetCommandPrompt("Select merged rebar mesh")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Mesh
    go.SubObjectSelect = False
    go.EnableUnselectObjectsOnExit(True)
    go.Get()
    if go.CommandResult() != Rhino.Commands.Result.Success:
        return go.CommandResult()

    mesh_obj = go.Object(0)
    mesh = mesh_obj.Mesh()

    # Keep mesh visible during workflow: remove persistent selection highlight.
    rs.UnselectAllObjects()
    sc.doc.Views.Redraw()
    if not mesh:
        Rhino.RhinoApp.WriteLine("No mesh selected.")
        return Rhino.Commands.Result.Failure

    spacing = rs.GetReal("Slice spacing (mm)", 50.0, 10.0, 1000.0)
    if spacing is None:
        return Rhino.Commands.Result.Cancel

    mesh.Normals.ComputeNormals()
    mesh.Compact()
    graph = _build_mesh_graph(mesh)

    parent_layer_index = sc.doc.Layers.CurrentLayerIndex

    rows = []
    name_counts = {}

    while True:
        result = _measure_once(mesh_obj, mesh, graph, spacing)
        if result is None:
            break
        if result is False:
            continue

        base_name = rs.GetString("Bar ID for this measurement", "BAR")
        if base_name is None:
            base_name = "BAR"
        unique_name = _unique_bar_name(base_name, name_counts)

        result["name"] = unique_name
        result["base_name"] = base_name
        result["spacing"] = spacing
        result["timestamp"] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

        layer_index = _ensure_child_layer(parent_layer_index, unique_name, System.Drawing.Color.OrangeRed)
        if layer_index < 0:
            Rhino.RhinoApp.WriteLine("Could not create sublayer {}; using active layer.".format(unique_name))
        for c in result.get("sample_curves", []):
            _add_curve_on_layer(c, layer_index)
        for pt in result.get("sample_points", []):
            _add_point_on_layer(pt, layer_index)
        sc.doc.Views.Redraw()

        rows.append(result)
        Rhino.RhinoApp.WriteLine("Saved measurement as: {}".format(unique_name))


    if rows:
        _export_csv(rows)

    return Rhino.Commands.Result.Success


if __name__ == "__main__":
    BarDiameterMVP()
