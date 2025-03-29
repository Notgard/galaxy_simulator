#pragma once

//SDL
#define WINDOW_WIDTH 1020
#define WINDOW_HEIGHT 1000

#define NORTH_EAST 0
#define SOUTH_EAST 1
#define SOUTH_WEST 2
#define NORTH_WEST 3

#define TOP_LEFT 0
#define TOP_RIGHT 1
#define BOTTOM_RIGHT 2
#define BOTTOM_LEFT 3

#define QUADTREE_MAX_VALUES 4
#define QUADTREE_THRESHOLD 16
#define QUADTREE_MAX_DEPTH 50

#define BOX_WIDTH 1000.0
#define BOX_HEIGHT 1000.0
#define BOX_LEFT 10.0
#define BOX_TOP 10.0

//simulation constants
#define delta_time 0.1
#define G 6.67430e-11
#define t_end (delta_time * 365 * 10)
#define ETA 0.01 //softening factor
#define DISTANCE_MULTIPLE 1e9 //meter
#define GFACTOR 2

#define PARTICLE_RADIUS 1.0f
#define PARTICLE_MASS 1.0
#define PARTICLE_POSITION_HISTORY_SIZE 10

#define SUN_RADIUS 10.0f
#define SUN_MASS 30000000/* 1.989e30 */

#define PLANET_RADIUS 5.0f

#define THETA 0.8

#define SMALL_GALAXY_THRESHOLD 2

#define CELESTIAL_BODY_COUNT 1

#define SUB_TREE_PROCESSING 1
#define CHUNKED_PROCESSING 2

#define DEBUG 0
