#include <stdio.h>
#include <stdlib.h>

#include "world_grid.h"
#include "misc/file.h"

WorldGrid WorldGrid::from_file(const char* filename, size_t X, size_t Y) {
  int length;
  char* contents = load_file(filename, &length);
  if (!contents) {
    fprintf(stderr, "Could not load %s!\n", filename);
    exit(1);
  }

  // find the width and height of the scenario
  int width = 0;
  int height = 0;
  for (int i = 0; i < length; ++i) {
    int row_start = i;
    while (i < length && contents[i] != '\n') {
      ++i;
    }
    ++height;
    int row_width = i - row_start;
    if (width < row_width) {
      width = row_width;
    }
  }

/*
  // allocate a grid large enough for the scenario
  // plus an extra ring of sinks around the outside
  WorldGrid grid(size_t(width)+2, size_t(height)+2);
*/

  if (size_t(width) > X-2 || size_t(height) > Y-2) {
    fprintf(stderr, "Failed to load %s.\n"
      "Scenario dimensions [%d,%d] exceed maximums [%zu,%zu]!\n",
      filename, width, height, X-2, Y-2);
    exit(1);
  }
  WorldGrid grid(X,Y);

  // initialize the grid
  // note that text files start at the top of the page and go down,
  // so our grid is filled from the top to bottom
  size_t x = 1;
  size_t y = grid.height() - 2;
  for (int i = 0; i < length; ++i, ++x) {
    char c = contents[i];
    if (c == '\n') {
      x = 0;
      --y;
    } else if (c == 'X') {
      grid.flag_as_solid({x,y});
    } else if (c == '0') {
      grid.flag_as_fluid({x,y});
    } else if (c == '?') {
      grid.flag_as_fluid({x,y});
      grid.flag_as_source({x,y});
    } else if (c == '=') {
      grid.flag_as_sink({x,y});
    }
  }

  release_file(contents);

  // add sinks around the outside so we have fewer edge cases
  for (size_t y = 0; y < grid.height(); ++y) {
    grid.flag_as_sink({0, y});
    grid.flag_as_sink({grid.width()-1, y});
  }
  for (size_t x = 0; x < grid.width(); ++x) {
    grid.flag_as_sink({x, 0});
    grid.flag_as_sink({x, grid.height()-1});
  }
  return grid;
}
