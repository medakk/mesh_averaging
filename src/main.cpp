#include "openvdb/openvdb.h"
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/VolumeToMesh.h"
#include "happly.h"
#include <iostream>
#include <string>
#include <vector>

openvdb::Vec3s centroid(const std::vector<openvdb::Vec3s>& points) {
  openvdb::Vec3s c(0.0, 0.0, 0.0);
  for(const auto& pt : points) {
    c += pt;
  }
  return c / points.size();
}

openvdb::Vec3s bounds(const std::vector<openvdb::Vec3s>& points) {
  openvdb::Vec3s min(1e9, 1e9, 1e9);
  openvdb::Vec3s max(-1e9, -1e9, -1e9);
  for(const auto& pt : points) {
    for(int i=0; i<3; i++) {
      min[i] = std::min(min[i], pt[i]);
      max[i] = std::max(max[i], pt[i]);
    }
  }
  return max - min;
}

openvdb::Vec3s min(const std::vector<openvdb::Vec3s>& points) {
  openvdb::Vec3s min(1e9, 1e9, 1e9);
  for(const auto& pt : points) {
    for(int i=0; i<3; i++) {
      min[i] = std::min(min[i], pt[i]);
    }
  }
  return min;
}
openvdb::Vec3s max(const std::vector<openvdb::Vec3s>& points) {
  openvdb::Vec3s max(-1e9, -1e9, -1e9);
  for(const auto& pt : points) {
    for(int i=0; i<3; i++) {
      max[i] = std::max(max[i], pt[i]);
    }
  }
  return max;
}

int main(int argc, char* argv[]) {
  openvdb::initialize();

  std::vector<std::string> filenames;
  for(int i=1; i<argc; i++) {
    filenames.emplace_back(argv[i]);
  }

  const auto test_filename = filenames[0];
  happly::PLYData plyIn(test_filename);

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> tris;
  std::vector<openvdb::Vec4I> quads;
  for(const auto pt : plyIn.getVertexPositions()) {
    points.push_back({(float)pt[0], (float)pt[1], (float)pt[2]});
  }
  for(const auto tri : plyIn.getFaceIndices()) {
    tris.push_back({(unsigned int)tri[0], (unsigned int)tri[1], (unsigned int)tri[2]});
  }

  const auto pts_min = min(points);
  const auto pts_max = max(points);
  std::cout << "min: " << pts_min << "\n";
  std::cout << "max: " << pts_max << "\n";

  openvdb::math::Transform xform; // ?
  xform.postScale(2);
  xform.postTranslate(min(points));
  const float narrow_band = 5;
  auto vdb_grid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(xform, points, tris, quads, narrow_band, narrow_band);
  std::cout << "number of voxels: " << vdb_grid->activeVoxelCount() << "\n";
  openvdb::io::File vdb_file("out.vdb");
  openvdb::GridPtrVec grids; grids.push_back(vdb_grid);
  vdb_file.write(grids);
  vdb_file.close();

  points.clear();
  tris.clear();
  quads.clear();
  openvdb::tools::volumeToMesh(*vdb_grid, points, tris, quads);

  std::cout << "number of output points: " << points.size() << "\n";
  std::cout << "number of output tris: " << tris.size() << "\n";
  std::cout << "number of output quads: " << quads.size() << "\n";

  happly::PLYData plyOut;
  std::vector<std::array<double, 3>> meshVertexPositions;
  std::vector<std::vector<size_t>> meshFaceIndices;
  for(const auto& pt : points) {
    meshVertexPositions.push_back({pt[0], pt[1], pt[2]});
  }
  for(const auto& tri : tris) {
    meshFaceIndices.push_back({tri[0], tri[1], tri[2]});
  }
  for(const auto& quad : quads) {
    meshFaceIndices.push_back({quad[0], quad[1], quad[2], quad[3]});
  }
  plyOut.addVertexPositions(meshVertexPositions);
  plyOut.addFaceIndices(meshFaceIndices);
  plyOut.write("out.ply", happly::DataFormat::ASCII);

  return 0;
}
