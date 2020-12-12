#include "openvdb/openvdb.h"
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/VolumeToMesh.h"
#include "happly.h"
#include <iostream>
#include <string>
#include <vector>

struct VDBMesh {
  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> tris;
  std::vector<openvdb::Vec4I> quads;
};

VDBMesh load_ply(const std::string& filename) {

  happly::PLYData plyIn(filename);

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> tris;
  std::vector<openvdb::Vec4I> quads;
  for(const auto pt : plyIn.getVertexPositions()) {
    points.push_back({(float)pt[0], (float)pt[1], (float)pt[2]});
  }
  for(const auto tri : plyIn.getFaceIndices()) {
    tris.push_back({(unsigned int)tri[0], (unsigned int)tri[1], (unsigned int)tri[2]});
  }

  return {
    std::move(points),
    std::move(tris),
    std::move(quads),
  };
}

void write_ply(const std::string& filename, const VDBMesh& vdb_mesh, bool flip_normals=true) {
  happly::PLYData plyOut;
  std::vector<std::array<double, 3>> meshVertexPositions;
  std::vector<std::vector<size_t>> meshFaceIndices;
  for(const auto& pt : vdb_mesh.points) {
    meshVertexPositions.push_back({pt[0], pt[1], pt[2]});
  }
  for(const auto& tri : vdb_mesh.tris) {
    if(flip_normals) {
      meshFaceIndices.push_back({tri[1], tri[0], tri[2]});
    } else {
      meshFaceIndices.push_back({tri[0], tri[1], tri[2]});
    }
  }
  for(const auto& quad : vdb_mesh.quads) {
    if(flip_normals) {
      meshFaceIndices.push_back({quad[0], quad[3], quad[2], quad[1]});
    } else {
      meshFaceIndices.push_back({quad[0], quad[1], quad[2], quad[3]});
    }
  }
  plyOut.addVertexPositions(meshVertexPositions);
  plyOut.addFaceIndices(meshFaceIndices);
  plyOut.write(filename, happly::DataFormat::ASCII);
}

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

  if(argc <= 1) {
    std::cerr << "usage: mesh_averaging [file1.ply] [file2.ply] [file3.ply]....\n";
    std::cerr << "Output is written to out.ply in the current folder.\n";
    return 1;
  }

  std::vector<std::string> filenames;
  for(int i=1; i<argc; i++) {
    filenames.emplace_back(argv[i]);
  }

  auto vdb_mesh0 = load_ply(filenames[0]);

  const auto pts_min = min(vdb_mesh0.points);
  const auto pts_max = max(vdb_mesh0.points);
  std::cout << "min: " << pts_min << "\n";
  std::cout << "max: " << pts_max << "\n";

  openvdb::math::Transform xform;
  xform.postScale(1); // arbitrary. reduce it for higher resolution
  xform.postTranslate(min(vdb_mesh0.points) - 20.0); // arbitrary subtract 20 to account for other samples
  const float narrow_band = 20;
  auto vdb_grid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
          xform, vdb_mesh0.points, vdb_mesh0.tris, vdb_mesh0.quads, narrow_band, narrow_band);
  std::cout << "number of voxels: " << vdb_grid->activeVoxelCount() << "\n";
  /*
  openvdb::io::File vdb_file("out.vdb");
  openvdb::GridPtrVec grids; grids.push_back(vdb_grid);
  vdb_file.write(grids);
  vdb_file.close();
  */

  // add other grids
  for(int i=1; i<filenames.size(); i++) {
    auto vdb_mesh_i = load_ply(filenames[i]);
    auto vdb_grid_i = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
            xform, vdb_mesh_i.points, vdb_mesh_i.tris, vdb_mesh_i.quads, narrow_band, narrow_band);
    vdb_grid->tree().combine(vdb_grid_i->tree(), [](const float& a, const float& b, float& result) {
      result = a + b;

      // uncomment this and remove the divide after the loop to do union instead.
      // result = std::min(a, b);
    });
  }
  // divide to get average
  openvdb::tools::foreach(vdb_grid->beginValueAll(), [&](const auto& iter) {
    iter.setValue(*iter / filenames.size());
  });






  // back to mesh
  vdb_mesh0.points.clear();
  vdb_mesh0.tris.clear();
  vdb_mesh0.quads.clear();
  openvdb::tools::volumeToMesh(*vdb_grid, vdb_mesh0.points, vdb_mesh0.tris, vdb_mesh0.quads);

  std::cout << "number of output points: " << vdb_mesh0.points.size() << "\n";
  std::cout << "number of output tris: " << vdb_mesh0.tris.size() << "\n";
  std::cout << "number of output quads: " << vdb_mesh0.quads.size() << "\n";

  write_ply("out.ply", vdb_mesh0);

  return 0;
}
