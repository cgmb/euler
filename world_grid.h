#pragma once

#include <stddef.h>
#include <vector>

#include "index_calculator.h"

class WorldGrid {
public:
  WorldGrid();
  WorldGrid(size_t x, size_t y);

  size_t width() const;
  size_t height() const;
  size_t size() const;

  void set_fluid(vec2zu index, bool value);

  void flag_as_fluid(vec2zu index);
  void flag_as_solid(vec2zu index);
  void flag_as_source(vec2zu index);
  void flag_as_sink(vec2zu index);

  bool is_fluid(vec2zu index) const;
  bool is_solid(vec2zu index) const;
  bool is_source(vec2zu index) const;
  bool is_sink(vec2zu index) const;

  static WorldGrid from_file(const char* filename, size_t X, size_t Y);

private:
  IndexCalculator ic_;
  std::vector<bool> fluid_;
  std::vector<bool> solid_;
  std::vector<bool> source_;
  std::vector<bool> sink_;
};

inline WorldGrid::WorldGrid()
  : ic_(0,0) 
{}

inline WorldGrid::WorldGrid(size_t x, size_t y)
  : ic_(x, y)
  , fluid_(ic_.size())
  , solid_(ic_.size())
  , source_(ic_.size())
  , sink_(ic_.size())
{}

inline size_t WorldGrid::width() const {
  return ic_.width();
}

inline size_t WorldGrid::height() const {
  return ic_.height();
}

inline size_t WorldGrid::size() const {
  return ic_.size();
}

inline void WorldGrid::set_fluid(vec2zu index, bool value) {
  fluid_[ic_.to1d(index)] = value;
}

inline void WorldGrid::flag_as_fluid(vec2zu index) {
  fluid_[ic_.to1d(index)] = true;
}

inline void WorldGrid::flag_as_solid(vec2zu index) {
  solid_[ic_.to1d(index)] = true;
}

inline void WorldGrid::flag_as_source(vec2zu index) {
  source_[ic_.to1d(index)] = true;
}

inline void WorldGrid::flag_as_sink(vec2zu index) {
  sink_[ic_.to1d(index)] = true;
}

inline bool WorldGrid::is_fluid(vec2zu index) const {
  return fluid_[ic_.to1d(index)];
}

inline bool WorldGrid::is_solid(vec2zu index) const {
  return solid_[ic_.to1d(index)];
}

inline bool WorldGrid::is_source(vec2zu index) const {
  return source_[ic_.to1d(index)];
}

inline bool WorldGrid::is_sink(vec2zu index) const {
  return sink_[ic_.to1d(index)];
}
