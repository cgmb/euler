#pragma once

/*
   These debug functions are used to print out the data structures used in the
   pressure solve. The output format is designed so output can be copy/pasted
   directly into Matlab / Octave for further analysis.

   The printing functions also handy just to understand how the data structures
   in memory relate to the matrixes and vectors in the math. The indexing for a
   lot of matrixes and vectors expect you to ignore fluid cells. I do that by
   using matrixes large enough to hold all cells and then just skipping over
   the non-fluid cells. I suspect that compacting and using some kind of
   mapping between indexes would be more efficient.

   This file doesn't really stand on its own. It expects to be embedded in
   main.cxx. Just include it right before whatever function is being debugged.
*/

// Fills q with the dense fluid matrix A.
void get_dense_a(int q[Y*X][Y*X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        q[X*y+x][X*y+x] = g_a[y][x].a_diag;
        if (x > 0)
          q[X*y+x][X*y+x-1] = get_a_minus_i(y,x);
        if (x+1 < X)
          q[X*y+x][X*y+x+1] = g_a[y][x].a_plus_i;
        if (y > 0)
          q[X*y+x][X*(y-1)+x] = get_a_minus_j(y,x);
        if (y+1 < X)
          q[X*y+x][X*(y+1)+x] = g_a[y][x].a_plus_j;
      }
    }
  }
}

// Prints a dense fluid matrix. e.g. A
void print_fluid_matrix(FILE* f, const char* name, int q[Y*X][Y*X]) {
  fprintf(f, "%s = [", name);
  for (size_t y0 = 0; y0 < Y; ++y0) {
    for (size_t x0 = 0; x0 < X; ++x0) {
      for (size_t y1 = 0; y1 < Y; ++y1) {
        for (size_t x1 = 0; x1 < X; ++x1) {
          if (is_water(y0,x0) && is_water(y1,x1)) {
            fprintf(f, "% 2d ", q[y0*X+x0][y1*X+x1]);
          }
        }
      }
      if (is_water(y0,x0)) {
        fprintf(f, "; ");
      }
    }
  }
  fprintf(f, "]\n");
}

// Prints a fluid vector
// Fluid vectors have one value per fluid cell, e.g. pressure, color, etc...
// Despite its appearance, q is a 1D vector with an index based on y and x.
// It's not quite as simple as 'idx = y*X + x', because cells that don't
// contain fluid are not part of the vector.
void print_fluid_vector(FILE* f, const char* name, double q[Y][X]) {
  fprintf(f, "%s = [", name);
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        fprintf(f, "%f ", q[y][x]);
      }
    }
  }
  fprintf(f, "].'\n");
}

// Prints a matrix
// This has the function same signature as print_fluid_vector, but in this case
// q really is a 2D matrix. All cells are part of the matrix, fluid or not.
void print_matrix(FILE* f, const char* name, float q[Y][X]) {
  fprintf(f, "%s = [", name);
  for (size_t y = Y; y--;) {
    for (size_t x = 0; x < X; ++x) {
      fprintf(f, "%f ", q[y][x]);
    }
    fprintf(f, ";\n");
  }
  fprintf(f, "]\n");
}
