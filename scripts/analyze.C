// Run interactively within root using: .x analyze.C

// Data types
using namespace ROOT::VecOps;
using tracksterRVec = ROOT::VecOps::RVec<ticl::Trackster>;
using simVertexRVec = ROOT::VecOps::RVec<SimVertex>;

// Define methods to produce new variables
auto doDistance = [](const simVertexRVec & v) {
  auto distanceV = v[0].position()-v[1].position();
  float distance = ceilf(distanceV.P() * 100)/100;
  return distance;
};
auto doCount = [](const tracksterRVec & v) {
  return (int)v.size();
};
auto doSumEnergy = [](const tracksterRVec & v) { 
  auto energies = Map(v, [](auto const & t) {return t.raw_energy();});
  return Sum(energies);
};
auto stdEnergy = [](const tracksterRVec & v) { 
  auto energies = Map(v, [](auto const & t) {return t.raw_energy();});
  return StdDev(energies);
};

void analyze() {

  std::string path = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/bainbrid/";
  path += "Photons/samples/FineCalo/CloseByPhotons/";
  std::vector<std::string> files;
  for ( int i = 0; i<100; ++i ) {
    std::stringstream ss;
    ss << path << "step3ticl_pt50_eta21_run" << i << ".root";
    files.push_back(ss.str());
  }
  
  // Read input file
  ROOT::RDataFrame R("Events",files);

  // Define which variables to produce
  auto Re = R
    .Define("dist", doDistance, {"SimVertexs_g4SimHits__SIM.obj"})
    .Define("nTotTS_EM", doCount, {"ticlTracksters_ticlTrackstersEM__RECO.obj"})
    .Define("eTotTS_EM", doSumEnergy, {"ticlTracksters_ticlTrackstersEM__RECO.obj"})
    .Define("nTotTS_Merge", doCount, {"ticlTracksters_ticlTrackstersMerge__RECO.obj"})
    .Define("eTotTS_Merge", doSumEnergy, {"ticlTracksters_ticlTrackstersMerge__RECO.obj"})
    .Define("nTotTS_EM3", doCount, {"ticlTracksters_ticlTrackstersEM3__TICL.obj"})
    .Define("eTotTS_EM3", doSumEnergy, {"ticlTracksters_ticlTrackstersEM3__TICL.obj"})
    .Define("nTotTS_CLUE3D", doCount, {"ticlTracksters_ticlTrackstersCLUE3D3__TICL.obj"})
    .Define("eTotTS_CLUE3D", doSumEnergy, {"ticlTracksters_ticlTrackstersCLUE3D3__TICL.obj"});

  // Create output file
  Re.Snapshot("Tracksters",
	      "tracksters.root",
	      {
		"dist",
		"nTotTS_EM",
		"eTotTS_EM",
		"nTotTS_Merge",
		"eTotTS_Merge",
		"nTotTS_EM3",
		"eTotTS_EM3",
		"nTotTS_CLUE3D",
		"eTotTS_CLUE3D",
		}
	      );
}
