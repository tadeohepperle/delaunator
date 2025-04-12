package example

import delaunay "../"
import "core:math"
import "core:math/linalg"
import "core:math/noise"
import "core:math/rand"
import "core:slice"
import q "shared:quat/quat"
import engine "shared:quat/quat/engine"


Vec2 :: [2]f32

Triangle :: [3]u32

main :: proc() {
	settings := engine.DEFAULT_ENGINE_SETTINGS
	settings.hot_reload_shaders = false
	settings.shaders_dir_path = ""
	engine.init(settings)
	defer engine.deinit()
	cam := engine.camera_controller_create()


	n_pts := 1000
	n_ptrs_sqrt := int(math.sqrt(f32(n_pts)))

	pts := make([]Vec2, n_pts)
	for &pt, i in pts {
		pt = Vec2{f32(i / n_ptrs_sqrt), f32(i % n_ptrs_sqrt)}
		offset := linalg.normalize(Vec2{rand.float32(), rand.float32()}) * 0.3
		pt += offset
	}

	triangulation := delaunay.triangulate(pts)
	triangles := slice.from_ptr(
		cast(^Triangle)raw_data(triangulation.triangles),
		len(triangulation.triangles) / 3,
	)

	for engine.next_frame() {
		engine.camera_controller_update(&cam)
		engine.draw_gizmos_coords()
		// engine.draw_gizmos_circle({3, 3}, engine.get_osc(7, 1, 2))

		for pt in pts {
			engine.draw_gizmos_circle(pt, 0.1)
		}

		for tri in triangles {
			a := pts[tri[0]]
			b := pts[tri[1]]
			c := pts[tri[2]]
			engine.draw_gizmos_triangle(a, b, c, q.ColorSoftOrange)

		}


	}


}
