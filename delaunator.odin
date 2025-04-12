package delaunator


import "core:math"
import "core:math/linalg"
import "core:slice"

Idx :: u32
// Configere this to f64 for double precision
Float :: f32
Point :: [2]Float

EMPTY :: max(Idx)

FloatAndIdx :: struct {
	float: Float,
	idx:   Idx,
}

// https://docs.rs/robust/latest/src/robust/lib.rs.html#1-2650

// Returns a positive value if the coordinates `pa`, `pb`, and `pc` occur in counterclockwise order
// (`pc` lies to the **left** of the directed line defined by coordinates `pa` and `pb`).  
// Returns a negative value if they occur in clockwise order (`pc` lies to the **right** of the directed line `pa, pb`).  
// Returns `0` if they are **collinear**. 
orient2d :: proc(a: Point, b: Point, c: Point) -> Float {
	detleft := (a.x - c.x) * (b.y - c.y)
	detright := (a.y - c.y) * (b.x - c.x)
	det := detleft - detright
	return det
}

circumdelta :: proc(a: Point, b: Point, c: Point) -> (offset_from_a: Point) {
	dx := b.x - a.x
	dy := b.y - a.y
	ex := c.x - a.x
	ey := c.y - a.y

	bl := dx * dx + dy * dy
	cl := ex * ex + ey * ey
	d := 0.5 / (dx * ey - dy * ex)

	x := (ey * bl - dy * cl) * d
	y := (dx * cl - ex * bl) * d
	return {x, y}
}
circumradius2 :: proc(a: Point, b: Point, c: Point) -> Float {
	off := circumdelta(a, b, c)
	return off.x * off.x + off.y * off.y
}
circumcenter :: proc(a: Point, b: Point, c: Point) -> Point {
	offset_from_a := circumdelta(a, b, c)
	return a + offset_from_a
}
nearly_equals :: proc(a: Point, b: Point) -> bool {
	return abs(a.x - b.x) <= linalg.F32_EPSILON * 2 && abs(a.y - b.y) <= linalg.F32_EPSILON * 2
}

in_circle :: proc(a: Point, b: Point, c: Point, p: Point) -> bool {
	dx := a.x - p.x
	dy := a.y - p.y
	ex := b.x - p.x
	ey := b.y - p.y
	fx := c.x - p.x
	fy := c.y - p.y

	ap := dx * dx + dy * dy
	bp := ex * ex + ey * ey
	cp := fx * fx + fy * fy

	return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0.0
}


/// Next halfedge in a triangle.
next_halfedge :: proc(i: Idx) -> Idx {
	if i % 3 == 2 {
		return i - 2
	} else {
		return i + 1
	}
}

/// Previous halfedge in a triangle.
prev_halfedge :: proc(i: Idx) -> Idx {
	if i % 3 == 0 {
		return i + 2
	} else {
		return i - 1
	}
}


Triangulation :: struct {
	//An array of point indices where each triple represents a Delaunay triangle.
	//  All triangles are directed counter-clockwise.
	triangles: [dynamic]Idx,
	// A vector of adjacent halfedge indices that allows traversing the triangulation graph.
	//
	// `i`-th half-edge in the array corresponds to vertex `triangles[i]`
	// the half-edge is coming from. `halfedges[i]` is the index of a twin half-edge
	// in an adjacent triangle (or `EMPTY` for outer half-edges on the convex hull).
	halfedges: [dynamic]Idx,
	// An array of indices that reference points on the convex hull of the triangulation, counter-clockwise.
	hull:      [dynamic]Idx,
}
// data structure for tracking the edges of the advancing convex hull
Hull :: struct {
	prev:   [dynamic]Idx,
	next:   [dynamic]Idx,
	tri:    [dynamic]Idx,
	hash:   [dynamic]Idx,
	start:  Idx,
	center: Point,
}


triangulation_make :: proc(n: int) -> Triangulation {
	max_triangles := 2 * n - 5 if n > 2 else 0
	return Triangulation {
		triangles = make([dynamic]Idx, 0, max_triangles * 3),
		halfedges = make([dynamic]Idx, 0, max_triangles * 3),
		hull = make([dynamic]Idx),
	}
}
triangulation_drop :: proc(this: ^Triangulation) {
	delete(this.triangles)
	delete(this.halfedges)
	delete(this.hull)
}

triangulation_n_triangles :: proc(this: Triangulation) -> int {
	return len(this.triangles) * 3
}

triangulation_add_triangle :: proc(
	this: ^Triangulation,
	i0: Idx,
	i1: Idx,
	i2: Idx,
	a: Idx,
	b: Idx,
	c: Idx,
) -> Idx {


	t := Idx(len(this.triangles))

	append(&this.triangles, i0)
	append(&this.triangles, i1)
	append(&this.triangles, i2)

	append(&this.halfedges, a)
	append(&this.halfedges, b)
	append(&this.halfedges, c)

	if a != EMPTY {
		this.halfedges[a] = t
	}
	if b != EMPTY {
		this.halfedges[b] = t + 1
	}
	if c != EMPTY {
		this.halfedges[c] = t + 2
	}

	return t
}


triangulation_legalize :: proc(this: ^Triangulation, a: Idx, points: []Point, hull: ^Hull) -> Idx {
	b := this.halfedges[a]

	// if the pair of triangles doesn't satisfy the Delaunay condition
	// (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
	// then do the same check/flip recursively for the new pair of triangles
	//
	//           pl                    pl
	//          /||\                  /  \
	//       al/ || \bl            al/    \a
	//        /  ||  \              /      \
	//       /  a||b  \    flip    /___ar___\
	//     p0\   ||   /p1   =>   p0\---bl---/p1
	//        \  ||  /              \      /
	//       ar\ || /br             b\    /br
	//          \||/                  \  /
	//           pr                    pr
	//

	ar := prev_halfedge(a)

	if b == EMPTY {
		return ar
	}

	al := next_halfedge(a)
	bl := prev_halfedge(b)

	p0 := this.triangles[ar]
	pr := this.triangles[a]
	pl := this.triangles[al]
	p1 := this.triangles[bl]

	illegal := in_circle(points[p0], points[pr], points[pl], points[p1])
	if illegal {
		this.triangles[a] = p1
		this.triangles[b] = p0

		hbl := this.halfedges[bl]
		har := this.halfedges[ar]

		// edge swapped on the other side of the hull (rare); fix the halfedge reference
		if hbl == EMPTY {
			e := hull.start
			for {
				if hull.tri[e] == bl {
					hull.tri[e] = a
					break
				}
				e = hull.prev[e]
				if e == hull.start {
					break
				}
			}
		}

		this.halfedges[a] = hbl
		this.halfedges[b] = har
		this.halfedges[ar] = bl

		if hbl != EMPTY {
			this.halfedges[hbl] = a
		}
		if har != EMPTY {
			this.halfedges[har] = b
		}
		if bl != EMPTY {
			this.halfedges[bl] = ar
		}

		br := next_halfedge(b)
		triangulation_legalize(this, a, points, hull)
		return triangulation_legalize(this, br, points, hull)
	}
	return ar
}


hull_make :: proc(n: int, center: Point, i0: Idx, i1: Idx, i2: Idx, points: []Point) -> Hull {
	hash_len := Idx(math.sqrt(Float(n)))

	hull := Hull {
		prev   = make([dynamic]Idx, n), // edge to prev edge
		next   = make([dynamic]Idx, n), // edge to next edge
		tri    = make([dynamic]Idx, n), // edge to adjacent halfedge
		hash   = make([dynamic]Idx, hash_len), // angular edge hash
		start  = i0,
		center = center,
	}
	for &h in hull.hash {
		h = EMPTY
	}

	hull.next[i0] = i1
	hull.prev[i2] = i1
	hull.next[i1] = i2
	hull.prev[i0] = i2
	hull.next[i2] = i0
	hull.prev[i1] = i0

	hull.tri[i0] = 0
	hull.tri[i1] = 1
	hull.tri[i2] = 2

	hull_hash_edge(&hull, points[i0], i0)
	hull_hash_edge(&hull, points[i1], i1)
	hull_hash_edge(&hull, points[i2], i2)


	return hull
}

hull_drop :: proc(this: ^Hull) {
	delete(this.prev)
	delete(this.next)
	delete(this.tri)
	delete(this.hash)
}

hull_hash_key :: proc(this: Hull, p: Point) -> u64 {
	dx := p.x - this.center.x
	dy := p.y - this.center.y

	p := dx / (abs(dx) + abs(dy))
	a: Float = (3.0 - p if dy > 0.0 else 1.0 + p) / 4.0 // [0..1]

	len := u64(len(this.hash))
	return u64(math.floor(Float(len) * a)) % len
}


hull_hash_edge :: proc(this: ^Hull, p: Point, i: Idx) {
	key := hull_hash_key(this^, p)
	this.hash[key] = i

}

hull_find_visible_edge :: proc(this: Hull, p: Point, points: []Point) -> (Idx, bool) {
	start: Idx = 0
	key := hull_hash_key(this, p)
	len := u64(len(this.hash))
	for j in 0 ..< len {
		start = this.hash[(key + j) % len]
		if start != EMPTY && this.next[start] != EMPTY {
			break
		}
	}
	start = this.prev[start]
	e := start

	for orient2d(p, points[e], points[this.next[e]]) <= 0 {
		e = this.next[e]
		if e == start {
			return EMPTY, false
		}
	}
	return e, e == start
}

calc_bbox_center :: proc(points: []Point) -> Point {
	min_x := max(Float)
	min_y := max(Float)
	max_x := min(Float)
	max_y := min(Float)
	for p in points {
		min_x = min(min_x, p.x)
		min_y = min(min_y, p.y)
		max_x = max(max_x, p.x)
		max_y = max(max_y, p.y)
	}
	return Point{(min_x + max_x) / 2.0, (min_y + max_y) / 2.0}
}

find_closest_point :: proc(points: []Point, p0: Point) -> (idx: Idx, ok: bool) {
	min_dist := max(Float)
	k: Idx = 0
	for p, i in points {
		d := linalg.length2(p0 - p)
		if d > 0.0 && d < min_dist {
			k = Idx(i)
			min_dist = d
		}
	}
	if min_dist == max(Float) {
		return max(Idx), false
	} else {
		return k, true
	}
}


find_seed_triangle :: proc(points: []Point) -> (idx_1: Idx, idx_2: Idx, idx_3: Idx, ok: bool) {
	// pick a seed point close to the center
	bbox_center := calc_bbox_center(points)
	i0 := find_closest_point(points, bbox_center) or_return
	p0 := points[i0]

	// find the point closest to the seed
	i1 := find_closest_point(points, p0) or_return
	p1 := points[i1]

	// find the third point which forms the smallest circumcircle with the first two
	min_radius := max(Float)
	i2: Idx = 0
	for p, i in points {
		i := Idx(i)
		if i == i0 || i == i1 {
			continue
		}
		r := circumradius2(p0, p1, p)
		if r < min_radius {
			i2 = i
			min_radius = r
		}
	}

	if min_radius == max(Float) {
		return {}, {}, {}, false
	} else {
		// swap the order of the seed points for counter-clockwise orientation
		if orient2d(p0, p1, points[i2]) > 0. {
			return i0, i2, i1, true
		} else {
			return i0, i1, i2, true
		}
	}
}


sortf :: proc(f: []FloatAndIdx) {
	slice.sort_by(f, proc(a: FloatAndIdx, b: FloatAndIdx) -> bool {
		return a.float < b.float
	})
}


/// Order collinear points by dx (or dy if all x are identical) and return the list as a hull
handle_collinear_points :: proc(points: []Point) -> Triangulation {
	if len(points) == 0 {
		return {}
	}

	first_pt := points[0]


	dists := make([]FloatAndIdx, len(points))
	defer delete(dists)
	for p, i in points {
		d := p.x - first_pt.x
		if d == 0 {
			d = p.y - first_pt.y
		}
		dists[i] = FloatAndIdx{d, Idx(i)}
	}
	sortf(dists)

	triangulation := triangulation_make(0)
	d0 := min(Float)
	for dist_and_i in dists {
		if dist_and_i.float > d0 {
			append(&triangulation.hull, dist_and_i.idx)
			d0 = dist_and_i.float
		}
	}

	return triangulation
}


// Triangulate a set of 2D points.
// Returns the triangulation for the input points.
// For the degenerated case when all points are collinear, returns an empty triangulation where all points are in the hull.
triangulate :: proc(points: []Point) -> Triangulation {
	i0, i1, i2, seed_triangle_ok := find_seed_triangle(points)
	if !seed_triangle_ok {
		return handle_collinear_points(points)
	}

	n := len(points)
	center := circumcenter(points[i0], points[i1], points[i2])

	triangulation := triangulation_make(n)
	triangulation_add_triangle(&triangulation, i0, i1, i2, EMPTY, EMPTY, EMPTY)

	// sort the points by distance from the seed triangle circumcenter
	dists := make([]FloatAndIdx, len(points))
	defer delete(dists)
	for p, i in points {
		dists[i] = FloatAndIdx{linalg.length2(center - p), Idx(i)}
	}
	sortf(dists)

	hull := hull_make(n, center, i0, i1, i2, points)
	defer hull_drop(&hull)

	for el, k in dists {

		i := el.idx
		p := points[i]

		// skip near-duplicates
		if k > 0 && nearly_equals(p, points[dists[k - 1].idx]) {
			continue
		}
		// skip seed triangle points
		if i == i0 || i == i1 || i == i2 {
			continue
		}

		// find a visible edge on the convex hull using edge hash
		e, walk_back := hull_find_visible_edge(hull, p, points)
		if e == EMPTY {
			continue // likely a near-duplicate point; skip it
		}

		// add the first triangle from the point
		t := triangulation_add_triangle(
			&triangulation,
			e,
			i,
			hull.next[e],
			EMPTY,
			EMPTY,
			hull.tri[e],
		)

		// recursively flip triangles from the point until they satisfy the Delaunay condition
		hull.tri[i] = triangulation_legalize(&triangulation, t + 2, points, &hull)
		hull.tri[e] = t // keep track of boundary triangles on the hull


		// walk forward through the hull, adding more triangles and flipping recursively
		n := hull.next[e]
		for {
			q := hull.next[n]
			if orient2d(p, points[n], points[q]) <= 0 {
				break
			}
			t := triangulation_add_triangle(
				&triangulation,
				n,
				i,
				q,
				hull.tri[i],
				EMPTY,
				hull.tri[n],
			)
			hull.tri[i] = triangulation_legalize(&triangulation, t + 2, points, &hull)
			hull.next[n] = EMPTY // mark as removed
			n = q
		}


		// walk backward from the other side, adding more triangles and flipping
		if walk_back {
			for {
				q := hull.prev[e]
				if orient2d(p, points[q], points[e]) <= 0 {
					break
				}
				t := triangulation_add_triangle(
					&triangulation,
					q,
					i,
					e,
					EMPTY,
					hull.tri[e],
					hull.tri[q],
				)
				triangulation_legalize(&triangulation, t + 2, points, &hull)
				hull.tri[q] = t
				hull.next[e] = EMPTY // mark as removed
				e = q
			}
		}

		// update the hull indices
		hull.prev[i] = e
		hull.next[i] = n
		hull.prev[n] = i
		hull.next[e] = i
		hull.start = e

		// save the two new edges in the hash table
		hull_hash_edge(&hull, p, i)
		hull_hash_edge(&hull, points[e], e)
	}

	// expose hull as a vector of point indices
	e := hull.start
	for {
		append(&triangulation.hull, e)
		e = hull.next[e]
		if e == hull.start {
			break
		}
	}

	shrink(&triangulation.triangles)
	shrink(&triangulation.halfedges)

	return triangulation
}
