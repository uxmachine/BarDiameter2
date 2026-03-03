"""
Rhino 8 MVP script: estimate rebar diameter from a noisy merged STL mesh.

Usage in Rhino 8:
1) Open PythonScript editor.
2) Load this file and run `BarDiameterMVP()`.
3) Pick one mesh, then two points on the same bar segment.

Output:
- Per-slice equivalent diameters.
- Mean / median diameter.
- 95% confidence interval (normal approximation).
"""

import math

import Rhino
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
    """Map picked point to a topology vertex robustly across Rhino versions."""
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


def BarDiameterMVP():
    go = Rhino.Input.Custom.GetObject()
    go.SetCommandPrompt("Select merged rebar mesh")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Mesh
    go.SubObjectSelect = False
    go.Get()
    if go.CommandResult() != Rhino.Commands.Result.Success:
        return go.CommandResult()

    mesh_obj = go.Object(0)
    mesh = mesh_obj.Mesh()
    if not mesh:
        Rhino.RhinoApp.WriteLine("No mesh selected.")
        return Rhino.Commands.Result.Failure

    p1 = rs.GetPointOnMesh(mesh_obj.ObjectId, "Pick start point on one bar")
    if not p1:
        return Rhino.Commands.Result.Cancel
    p2 = rs.GetPointOnMesh(mesh_obj.ObjectId, "Pick end point on the SAME bar (100-1000mm away)")
    if not p2:
        return Rhino.Commands.Result.Cancel

    spacing = rs.GetReal("Slice spacing (mm)", 50.0, 10.0, 1000.0)
    if spacing is None:
        return Rhino.Commands.Result.Cancel

    max_center_dist = _auto_max_center_dist(spacing)
    Rhino.RhinoApp.WriteLine("Using auto loop-center distance: {:.1f} mm".format(max_center_dist))

    mesh.Normals.ComputeNormals()
    mesh.Compact()

    s_topo = _closest_topology_vertex_index(mesh, p1)
    e_topo = _closest_topology_vertex_index(mesh, p2)
    if s_topo is None or e_topo is None:
        Rhino.RhinoApp.WriteLine("Could not map picked points to mesh topology vertices.")
        return Rhino.Commands.Result.Failure

    topo = mesh.TopologyVertices
    graph = _build_mesh_graph(mesh)
    path_ids = _dijkstra_path(mesh, graph, s_topo, e_topo)

    if path_ids and len(path_ids) >= 2:
        raw_pts = [topo[i] for i in path_ids]
        samples = _resample_polyline(raw_pts, spacing)
        Rhino.RhinoApp.WriteLine("Path mode: mesh topology path used.")
    else:
        samples = _line_samples(p1, p2, spacing)
        Rhino.RhinoApp.WriteLine("Path mode: straight fallback used (mesh path not connected).")

    diameters = []
    good_points = []

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
            if not pl[0].EpsilonEquals(pl[-1], 1e-6):
                pl.Add(pl[0])
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
        good_points.append(pt)

    n = len(diameters)
    if n < 3:
        Rhino.RhinoApp.WriteLine("Not enough valid slices. Repick cleaner points or increase slice spacing.")
        return Rhino.Commands.Result.Failure

    mean_d = sum(diameters) / n
    sorted_d = sorted(diameters)
    median_d = sorted_d[n // 2] if n % 2 == 1 else 0.5 * (sorted_d[n // 2 - 1] + sorted_d[n // 2])
    var = sum((x - mean_d) ** 2 for x in diameters) / max(1, n - 1)
    std = math.sqrt(var)
    half_ci = 1.96 * std / math.sqrt(n)

    Rhino.RhinoApp.WriteLine("--- BarDiameterMVP ---")
    Rhino.RhinoApp.WriteLine("Valid slices: {}".format(n))
    Rhino.RhinoApp.WriteLine("Diameters (mm): {}".format(", ".join("{:.2f}".format(x) for x in diameters)))
    Rhino.RhinoApp.WriteLine("Mean diameter (mm): {:.2f}".format(mean_d))
    Rhino.RhinoApp.WriteLine("Median diameter (mm): {:.2f}".format(median_d))
    Rhino.RhinoApp.WriteLine("Std dev (mm): {:.2f}".format(std))
    Rhino.RhinoApp.WriteLine("95% CI for mean (mm): [{:.2f}, {:.2f}]".format(mean_d - half_ci, mean_d + half_ci))

    for pt in good_points:
        sc.doc.Objects.AddPoint(pt)
    sc.doc.Views.Redraw()

    return Rhino.Commands.Result.Success


if __name__ == "__main__":
    BarDiameterMVP()
