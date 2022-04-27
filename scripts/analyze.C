using namespace ROOT::VecOps;
using tracksterRVec = ROOT::VecOps::RVec<ticl::Trackster>;
using simVertexRVec = ROOT::VecOps::RVec<SimVertex>;

auto doSumEnergy = [](const tracksterRVec & v) { 
  auto energies = Map(v, [](auto const & t) {return t.raw_energy();});
  return Sum(energies);
};
auto stdEnergy = [](const tracksterRVec & v) { 
  auto energies = Map(v, [](auto const & t) {return t.raw_energy();});
  return StdDev(energies);
};

auto doCount = [](const tracksterRVec & v) {
  return (int)v.size();
};

auto doDistance = [](const simVertexRVec & v) {
  auto distanceV = v[0].position()-v[1].position();
  float distance = ceilf(distanceV.P() * 100)/100;
  return distance;
};

void analyze() {
  ROOT::RDataFrame R("Events",
    "/Users/bainbrid/Desktop/HGCAL/step3ticl_pt50_eta21_run0.root"
    );
  auto Re = R
    .Define("trackstersSumEnergy", doSumEnergy, {"ticlTracksters_ticlTrackstersCLUE3D3__TICL.obj"})
    .Define("numTracksters", doCount, {"ticlTracksters_ticlTrackstersCLUE3D3__TICL.obj"})
    .Define("distance", doDistance, {"SimVertexs_g4SimHits__SIM.obj"});
  Re.Snapshot("Tracksters", "tracksters.root", 
      {
      "trackstersSumEnergy", "numTracksters",
      "distance"
      });
}
