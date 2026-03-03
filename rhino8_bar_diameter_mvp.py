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


def _closest_mesh_vertex_index(mesh, point3d):
    mp = mesh.ClosestMeshPoint(point3d, 0.0)
    if not mp:
        return None
    vi = mp.VertexIndex
    if vi is None or vi < 0:
        # Fallback: nearest vertex by distance
        best_i = None
        best_d2 = float("inf")
        for i, v in enumerate(mesh.Vertices):
            d2 = v.DistanceToSquared(point3d)
            if d2 < best_d2:
                best_d2 = d2
                best_i = i
        return best_i
    return vi


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

    # cumulative arc length
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
        out.append(a + alpha * (b - a))
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


def _equivalent_diameter_from_area(area):
    if area <= 0:
        return None
    return math.sqrt((4.0 * area) / math.pi)


def _curve_area(curve):
    amp = Rhino.Geometry.AreaMassProperties.Compute(curve)
    if not amp:
        return None
    return abs(amp.Area)


def _pick_relevant_loop(section_curves, sample_point, max_center_dist):
    best_curve = None
    best_dist = float("inf")

    for c in section_curves:
        if not c.IsClosed:
            continue
        ok, plane = c.TryGetPlane()
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
    p2 = rs.GetPointOnMesh(mesh_obj.ObjectId, "Pick end point on the SAME bar (100-500mm away)")
    if not p2:
        return Rhino.Commands.Result.Cancel

    spacing = rs.GetReal("Slice spacing (mm)", 50.0, 10.0, 200.0)
    if spacing is None:
        return Rhino.Commands.Result.Cancel
    max_center_dist = rs.GetReal("Max loop-center distance to path point (mm)", 30.0, 5.0, 100.0)
    if max_center_dist is None:
        return Rhino.Commands.Result.Cancel

    mesh.Normals.ComputeNormals()
    mesh.Compact()

    s = _closest_mesh_vertex_index(mesh, p1)
    e = _closest_mesh_vertex_index(mesh, p2)
    if s is None or e is None:
        Rhino.RhinoApp.WriteLine("Could not map picked points to mesh vertices.")
        return Rhino.Commands.Result.Failure

    topo = mesh.TopologyVertices
    # Convert mesh vertex indices to topology vertex indices when needed
    s_topo = mesh.TopologyVertices.TopologyVertexIndex(s) if s < mesh.Vertices.Count else s
    e_topo = mesh.TopologyVertices.TopologyVertexIndex(e) if e < mesh.Vertices.Count else e

    graph = _build_mesh_graph(mesh)
    path_ids = _dijkstra_path(mesh, graph, s_topo, e_topo)
    if not path_ids or len(path_ids) < 2:
        Rhino.RhinoApp.WriteLine("Failed to find path between picked points. Pick cleaner points on one bar.")
        return Rhino.Commands.Result.Failure

    raw_pts = [topo[i] for i in path_ids]
    samples = _resample_polyline(raw_pts, spacing)

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
        Rhino.RhinoApp.WriteLine("Not enough valid slices. Try larger max loop distance or pick cleaner points.")
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

    # Optional visualization: add points where valid slices were taken.
    for pt in good_points:
        sc.doc.Objects.AddPoint(pt)
    sc.doc.Views.Redraw()

    return Rhino.Commands.Result.Success


if __name__ == "__main__":
    BarDiameterMVP()
