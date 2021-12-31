// Copyright 2021 Global Phasing Ltd.
//
// Finding maxima or "blobs" in a Grid (map).
// Similar to CCP4 PEAKMAX and COOT's "Unmodelled blobs".

#ifndef GEMMI_BLOB_HPP_
#define GEMMI_BLOB_HPP_

#include "grid.hpp"   // for Grid

namespace gemmi {

struct Blob {
  double volume = 0.0;
  double score = 0.0;
  double peak_value = 0.0;
  gemmi::Position centroid;
  gemmi::Position peak_pos;
  explicit operator bool() const { return volume != 0. ; }
};

struct BlobCriteria {
  double min_volume = 10.0;
  double min_score = 15.0;
  double min_peak = 0.0;
  double cutoff;
};

namespace impl {

struct GridConstPoint {
  int u, v, w;
  float value;
};

inline Blob make_blob_of_points(const std::vector<GridConstPoint>& points,
                                const gemmi::Grid<float>& grid,
                                const BlobCriteria& criteria) {
  Blob blob;
  if (points.size() < 3)
    return blob;
  double volume_per_point = grid.unit_cell.volume / grid.point_count();
  double volume = points.size() * volume_per_point;
  if (volume < criteria.min_volume)
    return blob;
  double sum[3] = {0., 0., 0.};
  const GridConstPoint* peak_point = &points[0];
  blob.peak_value = points[0].value;
  double score = 0.;
  for (const GridConstPoint& point : points) {
    score += point.value;
    if (point.value > blob.peak_value) {
      blob.peak_value = point.value;
      peak_point = &point;
    }
    sum[0] += double(point.u) * point.value;
    sum[1] += double(point.v) * point.value;
    sum[2] += double(point.w) * point.value;
  }
  if (blob.peak_value < criteria.min_peak)
    return blob;
  blob.score = score * volume_per_point;
  if (blob.score < criteria.min_score)
    return blob;
  gemmi::Fractional fract(sum[0] / (score * grid.nu),
                          sum[1] / (score * grid.nv),
                          sum[2] / (score * grid.nw));
  blob.centroid = grid.unit_cell.orthogonalize(fract);
  blob.peak_pos = grid.get_position(peak_point->u, peak_point->v, peak_point->w);
  blob.volume = volume;
  return blob;
}

} // namespace impl

inline std::vector<Blob> find_blobs_by_flood_fill(const gemmi::Grid<float>& grid,
                                                  const BlobCriteria& criteria) {
  std::vector<Blob> blobs;
  std::array<std::array<int, 3>, 6> moves = {{{{-1, 0, 0}}, {{1, 0, 0}},
                                              {{0 ,-1, 0}}, {{0, 1, 0}},
                                              {{0, 0, -1}}, {{0, 0, 1}}}};
  // the mask will be used as follows:
  // -1=in blob,  0=in asu, not in blob (so far),  1=in neither
  std::vector<std::int8_t> mask = grid.get_asu_mask<std::int8_t>();
  std::vector<gemmi::GridOp> ops = grid.get_scaled_ops_except_id();
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        assert(idx == grid.index_q(u, v, w));
        if (mask[idx] != 0)
          continue;
        float value = grid.data[idx];
        if (value < criteria.cutoff)
          continue;
        std::vector<impl::GridConstPoint> points;
        points.push_back({u, v, w, value});
        mask[idx] = -1;
        for (size_t j = 0; j < points.size()/*increasing!*/; ++j)
          for (const std::array<int, 3>& mv : moves) {
            int nabe_u = points[j].u + mv[0];
            int nabe_v = points[j].v + mv[1];
            int nabe_w = points[j].w + mv[2];
            size_t nabe_idx = grid.index_s(nabe_u, nabe_v, nabe_w);
            if (mask[nabe_idx] == -1)
              continue;
            float nabe_value = grid.data[nabe_idx];
            if (nabe_value > criteria.cutoff) {
              if (mask[nabe_idx] != 0)
                for (const gemmi::GridOp& op : ops) {
                  auto t = op.apply(nabe_u, nabe_v, nabe_w);
                  size_t mate_idx = grid.index_s(t[0], t[1], t[2]);
                  if (mask[mate_idx] == 0)
                    mask[mate_idx] = 1;
                }
              mask[nabe_idx] = -1;
              points.push_back({nabe_u, nabe_v, nabe_w, nabe_value});
            }
          }
        if (Blob blob = impl::make_blob_of_points(points, grid, criteria))
          blobs.push_back(blob);
      }
  return blobs;
}

} // namespace gemmi
#endif
