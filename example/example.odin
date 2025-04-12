package example

import delaunator "../"
import "core:math"
import "core:math/linalg"
import "core:math/noise"
import "core:math/rand"
import "core:slice"
import rl "vendor:raylib"

Vec2 :: [2]f32
Triangle :: [3]u32

main :: proc() {
	screen_width: i32 = 1000
	screen_height: i32 = 800
	rl.InitWindow(screen_width, screen_height, "Delaunay Triangulation")
	defer rl.CloseWindow()

	rl.SetTargetFPS(60)
	rl.SetConfigFlags(rl.ConfigFlags{.VSYNC_HINT})

	// Camera setup
	camera := rl.Camera2D {
		offset   = rl.Vector2{f32(screen_width) / 2.0, f32(screen_height) / 2.0},
		target   = rl.Vector2{0, 0},
		rotation = 0,
		zoom     = 25,
	}

	N_PTS :: 1000
	n_ptrs_sqrt := int(math.sqrt(f32(N_PTS)))

	pts := make([]Vec2, N_PTS)
	for &pt, i in pts {
		pt = Vec2{f32(i / n_ptrs_sqrt), f32(i % n_ptrs_sqrt)} - f32(n_ptrs_sqrt / 2)
		offset := linalg.normalize(Vec2{rand.float32(), rand.float32()}) * 3.0
		pt += offset
	}

	triangulation := delaunator.triangulate(pts)
	triangles := slice.from_ptr(
		cast(^Triangle)raw_data(triangulation.triangles),
		len(triangulation.triangles) / 3,
	)

	for !rl.WindowShouldClose() {
		// Input
		ZOOM_SPEED :: 1.02
		if rl.IsKeyDown(rl.KeyboardKey.R) {
			camera.zoom *= ZOOM_SPEED
		} else if rl.IsKeyDown(rl.KeyboardKey.F) {
			camera.zoom /= ZOOM_SPEED
		}
		SPEED :: 0.3
		if rl.IsKeyDown(rl.KeyboardKey.RIGHT) || rl.IsKeyDown(rl.KeyboardKey.D) {
			camera.target[0] += SPEED
		} else if rl.IsKeyDown(rl.KeyboardKey.LEFT) || rl.IsKeyDown(rl.KeyboardKey.A) {
			camera.target[0] -= SPEED
		}
		if rl.IsKeyDown(rl.KeyboardKey.DOWN) || rl.IsKeyDown(rl.KeyboardKey.S) {
			camera.target[1] += SPEED
		} else if rl.IsKeyDown(rl.KeyboardKey.UP) || rl.IsKeyDown(rl.KeyboardKey.W) {
			camera.target[1] -= SPEED
		}

		rl.BeginDrawing()
		rl.ClearBackground(rl.BLACK)

		rl.BeginMode2D(camera)
		for pt in pts {
			rl.DrawCircleV(pt, 0.2, rl.PURPLE)
		}
		for tri in triangles {
			a := pts[tri[0]]
			b := pts[tri[1]]
			c := pts[tri[2]]
			rl.DrawTriangleLines(a, b, c, rl.ORANGE)
		}
		rl.EndMode2D()

		rl.DrawText("Use arrow keys or WASD to pan, R/F to zoom", 10, 10, 20, rl.RAYWHITE)

		rl.EndDrawing()
	}
}
