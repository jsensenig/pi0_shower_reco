////////////////////////////////////////////////////////////////////////
// Class:       Pi0Shower
// Plugin Type: analyzer (art v2_07_03)
// File:        Pi0Shower_module.cc
//
// Author: Leigh Whitehead using cetskelgen
//
// Example module that will access PFParticle objects tagged as beam
// particles by Pandora. A lot of useful information about these
// objects is stored in art via associations. These complex links
// are encapsulated by the ProtoDUNEPFParticleUtils class used here
// to simplify the process 
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"


#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNESliceUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

#include "lardata/ArtDataHelper/MVAReader.h"
// Cluster stuff
#include "larreco/RecoAlg/DBScanAlg.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

#include "art_root_io/TFileService.h"
#include "TFile.h"

// ROOT includes
#include "TVirtualFitter.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TROOT.h"

namespace protoana {
  class Pi0Shower;

  // Functions to fit a line to a cluster of points
  double Dist2D( double x1, double y1, double x2, double y2 );
  double distance2(double x, double y, double z, double * p);
  void line(double t, double * p, double & x, double & y, double & z);
  void SumDistance2(int &, double *, double & sum, double * par, int);
  TVector3 FitLine(const std::vector<TVector3> & input);

  struct ClusterProps {

    explicit ClusterProps(double r, double x, double y, double a, double aspan, double l, double q) : 
                          dist(r), dirX(x), dirY(y), angle(a), angle_span(aspan), len(l), charge(q) { }
    double dist;
    double dirX;
    double dirY;
    double angle;
    double angle_span;
    double len;
    double charge;
    TVector3 dirXY;
    std::vector<const recob::Hit*> hits;
  };

}

//protoana::ClusterProps::ClusterProps() : dist(0.), dirX(0.), dirY(0.), angle(0.), angle_span(0.), len(0.), charge(0.) { }

class protoana::Pi0Shower : public art::EDAnalyzer {
public:

  explicit Pi0Shower(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Shower(Pi0Shower const &) = delete;
  Pi0Shower(Pi0Shower &&) = delete;
  Pi0Shower & operator = (Pi0Shower const &) = delete;
  Pi0Shower & operator = (Pi0Shower &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;
  std::vector<art::Ptr<recob::Hit>> ClassifyHits( const art::Event &evt );
  std::map<size_t, std::map<size_t, std::vector<const recob::Hit*>>> ClusterHits( const art::Event &evt, std::vector<art::Ptr<recob::Hit>> &hitvec );
  void PolarClusterMerging( const art::Event &evt, std::map<size_t, std::map<size_t, std::vector<const recob::Hit*>>> &plane_cluster_map, std::pair<size_t, float> &beam_vertex );
  void MergeCluster( ClusterProps &up_cluster, ClusterProps &down_cluster );
  protoana::ClusterProps CharacterizeCluster( std::vector<const recob::Hit*> &clusterHits, std::pair<size_t, float> &beam_vertex );
  std::vector<int> GetHitPdg( const art::Event &evt, std::vector<const recob::Hit*> &hitvec );
  void TransformPoint( TVector3& point, const TVector3& shower_start, const TVector3& shower_dir );
  void reset();

private:

  // object that implements the DB scan algorithm
  cluster::DBScanAlg fDBScan;

  // fcl parameters
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fSliceTag;
  std::string fGeneratorTag;
  std::string fHitTag;

  // CNN Hit classifier cuts
  double fCollectionCnnCut;
  double fInductionCnnCut;

  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNESliceUtils sliceUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  art::ServiceHandle<geo::Geometry> geom;

  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  art::ServiceHandle< cheat::BackTrackerService > bt_serv;

  TTree *fTree;

  int run;
  int subrun;
  int event;

  /// Beam slice hits
  double beam_vertex_hit_channel, beam_vertex_hit_time;
  std::vector<double> selected_hits_channel, selected_hits_time;

  /// Clustered hits
  
  // Pre-cut/cluster
  std::vector<double> coll_hits_channel_initial, coll_hits_time_initial, coll_hits_charge_initial;
  std::vector<double> ind0_hits_channel_initial, ind0_hits_time_initial, ind0_hits_charge_initial;
  std::vector<double> ind1_hits_channel_initial, ind1_hits_time_initial, ind1_hits_charge_initial;

  // DBScan clusters
  std::vector<std::vector<int>> coll_hit_pdg_dbscan, ind0_hit_pdg_dbscan,  ind1_hit_pdg_dbscan;
  std::vector<std::vector<double>> coll_hits_channel_dbscan, coll_hits_time_dbscan, coll_hits_charge_dbscan;
  std::vector<std::vector<double>> ind0_hits_channel_dbscan, ind0_hits_time_dbscan, ind0_hits_charge_dbscan;
  std::vector<std::vector<double>> ind1_hits_channel_dbscan, ind1_hits_time_dbscan, ind1_hits_charge_dbscan;

  // Polar Merged clusters
  std::vector<std::vector<int>> coll_hit_pdg_polar_merge, ind0_hit_pdg_polar_merge,  ind1_hit_pdg_polar_merge;
  std::vector<std::vector<double>> coll_hits_channel_polar_merge, coll_hits_time_polar_merge, coll_hits_charge_polar_merge;
  std::vector<std::vector<double>> ind0_hits_channel_polar_merge, ind0_hits_time_polar_merge, ind0_hits_charge_polar_merge;
  std::vector<std::vector<double>> ind1_hits_channel_polar_merge, ind1_hits_time_polar_merge, ind1_hits_charge_polar_merge;

};


protoana::Pi0Shower::Pi0Shower(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fDBScan(p.get<fhicl::ParameterSet>("DBScanAlg")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fSliceTag(p.get<std::string>("SliceTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fCollectionCnnCut(p.get<double>("CollectionCNNCut")),
  fInductionCnnCut(p.get<double>("InductionCNNCut")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils"))
{

}

void protoana::Pi0Shower::analyze(art::Event const & evt)
{

  //reset containers
  reset();

  run    = evt.run();
  subrun = evt.subRun();
  event  = evt.id().event();

  std::cout << "************************************" << std::endl;
  std::cout << "Run: " << run << " Subrun: " << subrun << " Event: " << event << std::endl;

  // If this event is MC then we can check what the true beam particle is
  if ( !evt.isRealData() ) {
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    const simb::MCParticle* beam_particle = truthUtil.GetGeantGoodParticle( (*mcTruths)[0], evt );

    if( beam_particle == 0x0 ) { 
      std::cout << "Found no Geant particle!" << std::endl; 
      return; 
    }
    std::cout << "Beam particle PDG: " << beam_particle->PdgCode() << std::endl;
  }

  // Start by getting the reco PF beam particle
  // We can use its interaction vertex as a reference frame for the shower analysis
  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice( evt, fPFParticleTag );
  if( beamParticles.empty() ) { std::cout << "No beam particle found!" << std::endl; return; }
  std::vector<const recob::Hit*> beam_hits = pfpUtil.GetPFParticleHitsFromPlane( *beamParticles.at(0), evt, fPFParticleTag, 2); 

  // Get the beam particle interaction vertex, this translates to the max channel on the collection plane
  std::pair<size_t, float> max_channel_time(0, 0.0);
  for( auto hit : beam_hits ) if( hit->Channel() > max_channel_time.first ) max_channel_time = std::make_pair( hit->Channel(), hit->PeakTime() );
  beam_vertex_hit_channel = max_channel_time.first;
  beam_vertex_hit_time = max_channel_time.second;

  // Step 1. Classify hit with CNN
  auto selected_hits = ClassifyHits( evt );

  // 2. Cluster selected hits into shower segments using DBScan
  auto clustered_hit_map = ClusterHits( evt, selected_hits );

  // Step 3. First stage merging of clusters
  PolarClusterMerging( evt, clustered_hit_map, max_channel_time );

  fTree -> Fill(); 

}


// Step 1. Select only the Hits which are EM-like as classified by the CNN
std::vector<art::Ptr<recob::Hit>> protoana::Pi0Shower::ClassifyHits( const art::Event &evt ) {

  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit, 4> * CNN_results = 0x0;
  CNN_results = new anab::MVAReader<recob::Hit, 4>(evt, "emtrkmichelid:emtrkmichel" );

  // Get handle to the hits in the event
  //auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(fHitTag);

  // Get all hits from the beam slice
  // TODO I assume the cosmics are removed from the slice but need to check
  unsigned short beam_slice = pfpUtil.GetBeamSlice(evt, fPFParticleTag);
  const std::vector<art::Ptr<recob::Hit>> sliceHits = sliceUtil.GetRecoSliceHits_Ptrs(beam_slice, evt, fSliceTag);

  // The hits which pass our cuts will be stored here
  std::vector<art::Ptr<recob::Hit>> selected_hits;

  // Note: 
  // collection plane = 2
  // two induction planes = 0,1

  // Loop over all the Hits and make a cut on the EM-likeness of the Hit
  for ( size_t h = 0; h < sliceHits.size(); ++h ) {

    art::Ptr<recob::Hit> hit = sliceHits[h];
    std::array<float,4> cnn_out = CNN_results->getOutput( hit );
    double em_score     = cnn_out[ CNN_results->getIndex("em") ];
    //double michel_score = cnn_out[ CNN_results->getIndex("michel") ];

    // Make a seperate cut on collection and induction planes
    // the induction planes are noisier than collection so make a tighter cut on them e.g. >0.9 ind. and >0.5 coll.
    size_t thePlane = hit.get()->WireID().asPlaneID().Plane;
    
    // Save all the intial hits, pre-cut by plane
    switch( thePlane ) {
      case 0: {
                ind0_hits_channel_initial.emplace_back( hit.get()->Channel() );
                ind0_hits_time_initial.emplace_back( hit.get()->PeakTime() );
                ind0_hits_charge_initial.emplace_back( hit.get()->Integral() );
                break;
              }  
      case 1: {
                ind1_hits_channel_initial.emplace_back( hit.get()->Channel() );
                ind1_hits_time_initial.emplace_back( hit.get()->PeakTime() );
                ind1_hits_charge_initial.emplace_back( hit.get()->Integral() );
                break;
              }
      case 2: {
                coll_hits_channel_initial.emplace_back( hit.get()->Channel() );
                coll_hits_time_initial.emplace_back( hit.get()->PeakTime() );
                coll_hits_charge_initial.emplace_back( hit.get()->Integral() );
                break;
              }
      default: 
               std::cout << "Unknown plane! " << thePlane << std::endl;
               break;
    } //switch

    if( (thePlane == 2 && em_score < fCollectionCnnCut) || (thePlane != 2 && em_score < fInductionCnnCut) ) continue;

    selected_hits.emplace_back( hit );

    if( thePlane == 2 ) {
      selected_hits_channel.emplace_back( hit.get()->Channel() );
      selected_hits_time.emplace_back( hit.get()->PeakTime() );
    }

  }

  return selected_hits;

}

// Step 2. Cluster the selected hits in time/channel for each plane
std::map<size_t, std::map<size_t, std::vector<const recob::Hit*>>> protoana::Pi0Shower::ClusterHits( const art::Event &evt, std::vector<art::Ptr<recob::Hit>> &hitvec ) {

  // Try DBScan to cluster hits in channel/time 

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
  lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();
  lariov::ChannelStatusProvider::ChannelSet_t const BadChannels = channelStatus.BadChannels();

  // the vector to fill
  std::vector<art::Ptr<recob::Hit>> allhits;

  // make a map of the geo::PlaneID (collection or induction) to vectors of art::Ptr<recob::Hit>
  std::map<size_t, std::map<size_t, std::vector<const recob::Hit*>>> plane_cluster_map;
  std::map<size_t, std::vector<art::Ptr<recob::Hit>>> planeIDToHits;
  for (size_t i = 0; i < hitvec.size(); ++i) planeIDToHits[hitvec.at(i).get()->WireID().asPlaneID().Plane].push_back(hitvec[i]);

  // Loop over each wire plane and cluster
  for (auto& itr : planeIDToHits) {

    allhits.resize(itr.second.size());
    allhits.swap(itr.second);
    // Now initialize and run the clustering on the plane
    fDBScan.InitScan(clock_data, det_prop, allhits, BadChannels);
    fDBScan.run_cluster();

    std::cout << "Found " << fDBScan.fclusters.size() << " clusters!" << std::endl;

    // Now fill a map with the clusters of hits for each cluster and plane: the map = <planeID, <clusterID, Hit_vector>>
    for (size_t i = 0; i < fDBScan.fclusters.size(); ++i) {
      for (size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j) {
        if (fDBScan.fpointId_to_clusterId[j] == i) {
          plane_cluster_map[itr.first][i].emplace_back( allhits[j].get() );
        }
      }
    }
  }

  // Label hits and fill the histograms
  for( auto &plane : plane_cluster_map ) { // wire plane
    for( auto &cluster : plane.second ) { // cluster
      // Select the plane
      switch( plane.first ) {
        case 0: {
                  ind0_hit_pdg_dbscan.push_back( GetHitPdg( evt, cluster.second ) );
                  ind0_hits_channel_dbscan.push_back( std::vector<double>() );
                  ind0_hits_time_dbscan.push_back( std::vector<double>() );
                  ind0_hits_charge_dbscan.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.second ) {
                    ind0_hits_channel_dbscan.back().emplace_back( hit->Channel() );
                    ind0_hits_time_dbscan.back().emplace_back( hit->PeakTime() );
                    ind0_hits_charge_dbscan.back().emplace_back( hit->Integral() );
                  }
                  break;
                }  
        case 1: {
                  ind1_hit_pdg_dbscan.push_back( GetHitPdg( evt, cluster.second ) );
                  ind1_hits_channel_dbscan.push_back( std::vector<double>() );
                  ind1_hits_time_dbscan.push_back( std::vector<double>() );
                  ind1_hits_charge_dbscan.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.second ) {
                    ind1_hits_channel_dbscan.back().emplace_back( hit->Channel() );
                    ind1_hits_time_dbscan.back().emplace_back( hit->PeakTime() );
                    ind1_hits_charge_dbscan.back().emplace_back( hit->Integral() );
                  }
                  break;
                }
        case 2: {
                  coll_hit_pdg_dbscan.push_back( GetHitPdg( evt, cluster.second ) );
                  coll_hits_channel_dbscan.push_back( std::vector<double>() );
                  coll_hits_time_dbscan.push_back( std::vector<double>() );
                  coll_hits_charge_dbscan.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.second ) {                               
                    coll_hits_channel_dbscan.back().emplace_back( hit->Channel() );
                    coll_hits_time_dbscan.back().emplace_back( hit->PeakTime() );
                    coll_hits_charge_dbscan.back().emplace_back( hit->Integral() );
                  }
                  break;
                }
        default: 
                 std::cout << "Unknown plane! " << plane.first << std::endl;
                 break;
      } //switch
    } //cluster
  } //plane

  return plane_cluster_map;

}


// Step 3. Cluster shower segments using polar coordinate merging (map is <planeID, <clusterID, Hit_vector>> )
void protoana::Pi0Shower::PolarClusterMerging( const art::Event &evt, std::map<size_t, std::map<size_t, std::vector<const recob::Hit*>>> &plane_cluster_map, std::pair<size_t, float> &beam_vertex ) {

  std::map<size_t, std::vector<ClusterProps>> plane_cluster;

  // Loop over each (3) wire plane
  for( auto &plane : plane_cluster_map ) {
    //if( plane.first != 2 ) continue; // only collection plane for now
    for( auto &cluster : plane.second ) { // Loop over each cluster in the plane
      plane_cluster[plane.first].push_back( CharacterizeCluster(cluster.second, beam_vertex) );
    }
  }

  // Sort the clusters from closest (upstream) to farthest (downstream) from the beam vertex for each of the wire planes
  for( auto plane_cls : plane_cluster ) {
    std::sort(plane_cls.second.begin(), plane_cls.second.end(), [] ( ClusterProps& a, ClusterProps& b ) { return a.dist < b.dist; });
  }

  // Now try to merge shower segments under 3 conditions (upstream means closest to beam vertex)
  // 1. Upstream cluster Q > downstream cluster
  // 2. Downstream cluster angle within angle span of upstream cluster
  // 3. Distance between up/downstream clusters < length of upstream cluster

  // Merge clusters for each wire plane
  for( auto &clusters : plane_cluster ) {
    // Loop over the merging until no more showers are merged  
    bool still_merging = true;
    while( still_merging ) {
      still_merging = false;
                                                                                                                                          
      // Try to merge the "jth" cluster into the "ith" cluster
      // The "ith" cluster is always upstream of the "jth" cluster
      for( size_t up = 0; up < clusters.second.size(); up++ ) {
        for( size_t down = up+1; down < clusters.second.size(); down++ ) { 
                                                                                                                                          
          double halfspan = 0.5 * clusters.second.at(up).angle_span;
          if( clusters.second.at(up).charge < clusters.second.at(down).charge ) continue; // (1.)
          if( ((clusters.second.at(up).angle - halfspan) > clusters.second.at(down).angle) &&
              ((clusters.second.at(up).angle + halfspan) < clusters.second.at(down).angle) ) continue; // (2.)
          //////////////////////
          if( clusters.second.at(up).dirXY.Dot(clusters.second.at(down).dirXY) < 0.9) continue;
          //////////////////////
          if( (clusters.second.at(down).dist - clusters.second.at(up).dist) > clusters.second.at(up).len ) continue; // (3.)
                                                                                                                                         
          // Merge the downstream into upstream cluster, MergeCluster() reference arguemnts avoids a copy of the merged structure
          MergeCluster( clusters.second.at(up), clusters.second.at(down) );
          clusters.second.erase(clusters.second.begin() + down); // now remove the merged cluster
                                                                                                                                          
          std::cout << "Cluster merged on plane " << clusters.first << " Remaining clusters: " << clusters.second.size() << std::endl;
          still_merging = true;
                                                                                                                                          
        } // jth loop
      } // ith loop
    } // while loop
  } // plane loop

  // Label hits and fill the histograms
  for( auto &plane : plane_cluster ) { // wire plane
    for( auto &cluster : plane.second ) { // cluster
      // Select the plane
      switch( plane.first ) {
        case 0: {
                  ind0_hit_pdg_polar_merge.push_back( GetHitPdg( evt, cluster.hits ) );
                  ind0_hits_channel_polar_merge.push_back( std::vector<double>() );
                  ind0_hits_time_polar_merge.push_back( std::vector<double>() );
                  ind0_hits_charge_polar_merge.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.hits ) {
                    ind0_hits_channel_polar_merge.back().emplace_back( hit->Channel() );
                    ind0_hits_time_polar_merge.back().emplace_back( hit->PeakTime() );
                    ind0_hits_charge_polar_merge.back().emplace_back( hit->Integral() );
                  }
                  break;
                }  
        case 1: {
                  ind1_hit_pdg_polar_merge.push_back( GetHitPdg( evt, cluster.hits ) );
                  ind1_hits_channel_polar_merge.push_back( std::vector<double>() );
                  ind1_hits_time_polar_merge.push_back( std::vector<double>() );
                  ind1_hits_charge_polar_merge.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.hits ) {
                    ind1_hits_channel_polar_merge.back().emplace_back( hit->Channel() );
                    ind1_hits_time_polar_merge.back().emplace_back( hit->PeakTime() );
                    ind1_hits_charge_polar_merge.back().emplace_back( hit->Integral() );
                  }
                  break;
                }
        case 2: {
                  coll_hit_pdg_polar_merge.push_back( GetHitPdg( evt, cluster.hits ) );
                  coll_hits_channel_polar_merge.push_back( std::vector<double>() );
                  coll_hits_time_polar_merge.push_back( std::vector<double>() );
                  coll_hits_charge_polar_merge.push_back( std::vector<double>() ); 
                  for( auto &hit : cluster.hits ) {                               
                    coll_hits_channel_polar_merge.back().emplace_back( hit->Channel() );
                    coll_hits_time_polar_merge.back().emplace_back( hit->PeakTime() );
                    coll_hits_charge_polar_merge.back().emplace_back( hit->Integral() );
                  }
                  break;
                }
        default: 
                 std::cout << "Unknown plane! " << plane.first << std::endl;
                 break;
      } //switch
    } //cluster
  } //plane
 
}

// Essentially the += operator for merging the ClusterProps structure
void protoana::Pi0Shower::MergeCluster( ClusterProps &up_cluster, ClusterProps &down_cluster ) {

  // No need to re-fill these variables

  //up_cluster.dist 
  //up_cluster.angle // re-fit the cluster Hits
 
  // TODO should recalculate angle_span
  //up_cluster.angle_span
  
  up_cluster.len += down_cluster.len;
  up_cluster.charge += down_cluster.charge;
  up_cluster.hits.insert(up_cluster.hits.begin(), down_cluster.hits.begin(), down_cluster.hits.end());

  // Unfortunately we need to redo the charge-weighted Hit fit`
  std::vector<TVector3> tvec;
  for( auto hit : up_cluster.hits ) tvec.emplace_back( TVector3(hit->Channel(), hit->PeakTime(), hit->Integral()) );
  TVector3 newdir = FitLine( tvec );

  up_cluster.dirX = newdir.X();
  up_cluster.dirY = newdir.Y();
  up_cluster.dirXY = newdir;

}

// Characterize the clusters, beam_vertex = <PeakTime, Channel>
protoana::ClusterProps protoana::Pi0Shower::CharacterizeCluster( std::vector<const recob::Hit*> &clusterHits, std::pair<size_t, float> &beam_vertex ) {

  double clusterQ = 0.;
  double rmin = 1000., rmax = 0.;
  double anglemin = 1000., anglemax = 0.;
  std::vector<TVector3> fit_hits_vec(clusterHits.size());

  for( auto &hit : clusterHits ) {
    double ptime   = static_cast<double>(hit->PeakTime()); // from float
    double channel = static_cast<double>(hit->Channel());  // from size_t
    double charge  = static_cast<double>(hit->Integral()); // from float

    // Create 3D vectors with x,y,z = time,channel,charge
    fit_hits_vec.emplace_back( TVector3(ptime, channel, charge) );

    // Get shower start/end
    double r = Dist2D(beam_vertex.first, beam_vertex.second, channel, ptime);
    if( r < rmin ) rmin = r; // shower start
    if( r > rmax ) rmax = r; // shower end

    double angle = TMath::ACos( beam_vertex.first / r ); // angle wrt x-axis (channel)
    if( angle < anglemin ) anglemin = angle; 
    if( angle > anglemax ) anglemax = angle; 

    // Get the shower total charge
    clusterQ += charge;
  }

  // Fit the (Time, Channel, Charge) points to get a charge-weighted direction for the shower
  TVector3 shower_segment_dir = FitLine( fit_hits_vec );
  double r = Dist2D(beam_vertex.first, beam_vertex.second, shower_segment_dir.X(), shower_segment_dir.Y());
  double cluster_angle = TMath::ACos( shower_segment_dir.X() / r );

  // Shower length
  double len = rmax - rmin;

  // Shower angle span
  double angle_span = anglemax - anglemin;

  // return ClusterProps structure of the shower segement properties
  ClusterProps return_props( rmin, shower_segment_dir.X(), shower_segment_dir.Y(), cluster_angle, angle_span, len, clusterQ );
  return_props.hits = clusterHits;
  return_props.dirXY = shower_segment_dir;
  return return_props;

}

// Function to label the PDG of the Hits
std::vector<int> protoana::Pi0Shower::GetHitPdg( const art::Event &evt, std::vector<const recob::Hit*> &hitvec ) {

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

  // Function that loops over all hits in an event and returns those that an MCParticle
  // contributed to.
  bool use_eve = true;
  std::vector<int> hit_pdg_vec;
 
  // Backtrack all hits to verify whether they belong to the current MCParticle.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // Loop over each hit and attempt to get the MC particle which createed it. If no MC particle give it PDG = -999
  for(const recob::Hit* hit : hitvec) {
    if (use_eve) {
      for(const sim::TrackIDE & ide : bt_serv->HitToEveTrackIDEs(clockData, *hit)) {
        int trackId = ide.trackID;
        const simb::MCParticle* mcpart = 0x0;
        mcpart = pi_serv->TrackIdToParticle_P(trackId);
        if( mcpart != 0x0 ) hit_pdg_vec.emplace_back( mcpart->PdgCode() );
        else hit_pdg_vec.emplace_back( -999 );
      }
    } else {
      for(const int trackId : bt_serv->HitToTrackIds(clockData, *hit)) {
        const simb::MCParticle* mcpart = 0x0;
        mcpart = pi_serv->TrackIdToParticle_P(trackId);
        if( mcpart != 0x0 ) hit_pdg_vec.emplace_back( mcpart->PdgCode() );
        else hit_pdg_vec.emplace_back( -999 );
      }
    }
  }

  return hit_pdg_vec;

}

// Transform a point to a shifted (shower_start origin) and rotated (z-axis parallel to shower_dir)
void protoana::Pi0Shower::TransformPoint( TVector3& point, const TVector3& shower_start, const TVector3& shower_dir ) {

  // Original coordinate system,
  // x-axis: horizontal, transverse to beam direction; y-axis: vertical; z-axis: beam direction

  // Now shift coordinate system
  point -= shower_start;               // shift point to origin which is shower-start
  point.RotateX( shower_dir.Theta() ); // rotate in polar angle 
  point.RotateY( -shower_dir.Phi() );  // rotate in azimuthal angle (note: we're on the negative side of x-axis)

}

// 2D distance helper function
double protoana::Dist2D( double x1, double y1, double x2, double y2 ) {
  return sqrt( pow((x1 - x2), 2) + pow((y1 - y2), 2) );
}

// Fit a line to short shower segments with least square fit.
// Copied from PDSP_Analyzer, thanks Jake!

void protoana::line(double t, double *p, double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

// calculate distance line-point
double protoana::distance2(double x,double y,double z, double *p) {
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   ROOT::Math::XYZVector xp(x,y,z);
   ROOT::Math::XYZVector x0(p[0], p[2], 0.);
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1.);
   ROOT::Math::XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();
   return d2;
}

// function to be minimized
void protoana::SumDistance2(int &, double *, double & sum, double * par, int) {
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
   }
}

// Returns a short vector pointing in the direction of the shower segment
TVector3 protoana::FitLine(const std::vector<TVector3> & input) {

  TGraph2D * gr = new TGraph2D();
  for (size_t i = 0; i < input.size(); ++i) {
    gr->SetPoint(i, input[i].X(), input[i].Y(), input[i].Z());
  }

  TVirtualFitter * min = TVirtualFitter::Fitter(0,4);
  min->SetObjectFit(gr);
  min->SetFCN(SumDistance2);

  double arglist[10];

  // set minimum print level
  arglist[0] = -1;
  min->ExecuteCommand("SET PRINT",arglist,1);


  double pStart[4] = {1,1,1,1};
  min->SetParameter(0,"x0",pStart[0],0.01,0,0);
  min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  min->SetParameter(2,"y0",pStart[2],0.01,0,0);
  min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  min->ExecuteCommand("MIGRAD", arglist, 2);

  // get fit parameters
  double parFit[4];
  for (int i = 0; i < 4; ++i) parFit[i] = min->GetParameter(i);

  double startX1, startY1, startZ1;
  double startX2, startY2, startZ2;
  line(0, parFit, startX1, startY1, startZ1);
  line(1, parFit, startX2, startY2, startZ2);

  TVector3 diff(startX2 - startX1, startY2 - startY1, startZ2 - startZ1);

  delete gr;
  delete min;
  return diff;
}

void protoana::Pi0Shower::beginJob()
{

  gROOT->SetBatch(1);
 
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("points","points tree");
 
  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
 
  /// Beam slice hits
  fTree->Branch("beam_vertex_hit_channel", &beam_vertex_hit_channel);
  fTree->Branch("beam_vertex_hit_time", &beam_vertex_hit_time);
  fTree->Branch("selected_hits_channel", &selected_hits_channel);
  fTree->Branch("selected_hits_time", &selected_hits_time);

  /// Hits ///
  
  // Hits per cluster for each plane
  fTree->Branch("coll_hits_channel_initial", &coll_hits_channel_initial);
  fTree->Branch("coll_hits_time_initial", &coll_hits_time_initial);
  fTree->Branch("coll_hits_charge_initial", &coll_hits_charge_initial);
  fTree->Branch("ind0_hits_channel_initial", &ind0_hits_channel_initial);
  fTree->Branch("ind0_hits_time_initial", &ind0_hits_time_initial);
  fTree->Branch("ind0_hits_charge_initial", &ind0_hits_charge_initial);
  fTree->Branch("ind1_hits_channel_initial", &ind1_hits_channel_initial);
  fTree->Branch("ind1_hits_time_initial", &ind1_hits_time_initial);
  fTree->Branch("ind1_hits_charge_initial", &ind1_hits_charge_initial);

  // DBScan step
  // Hit PDG by plane                                                                                                                   
  fTree->Branch("coll_hit_pdg_dbscan", &coll_hit_pdg_dbscan);
  fTree->Branch("ind0_hit_pdg_dbscan", &ind0_hit_pdg_dbscan);
  fTree->Branch("ind1_hit_pdg_dbscan", &ind1_hit_pdg_dbscan);
  // Hits per cluster for each plane
  fTree->Branch("coll_hits_channel_dbscan", &coll_hits_channel_dbscan);
  fTree->Branch("coll_hits_time_dbscan", &coll_hits_time_dbscan);
  fTree->Branch("coll_hits_charge_dbscan", &coll_hits_charge_dbscan);
  fTree->Branch("ind0_hits_channel_dbscan", &ind0_hits_channel_dbscan);
  fTree->Branch("ind0_hits_time_dbscan", &ind0_hits_time_dbscan);
  fTree->Branch("ind0_hits_charge_dbscan", &ind0_hits_charge_dbscan);
  fTree->Branch("ind1_hits_channel_dbscan", &ind1_hits_channel_dbscan);
  fTree->Branch("ind1_hits_time_dbscan", &ind1_hits_time_dbscan);
  fTree->Branch("ind1_hits_charge_dbscan", &ind1_hits_charge_dbscan);

  // Polar Clustering step
  // Hit PDG by plane
  fTree->Branch("coll_hit_pdg_polar_merge", &coll_hit_pdg_polar_merge);
  fTree->Branch("ind0_hit_pdg_polar_merge", &ind0_hit_pdg_polar_merge);
  fTree->Branch("ind1_hit_pdg_polar_merge", &ind1_hit_pdg_polar_merge);
  // Hits per cluster for each plane
  fTree->Branch("coll_hits_channel_polar_merge", &coll_hits_channel_polar_merge);
  fTree->Branch("coll_hits_time_polar_merge", &coll_hits_time_polar_merge);
  fTree->Branch("coll_hits_charge_polar_merge", &coll_hits_charge_polar_merge);
  fTree->Branch("ind0_hits_channel_polar_merge", &ind0_hits_channel_polar_merge);
  fTree->Branch("ind0_hits_time_polar_merge", &ind0_hits_time_polar_merge);
  fTree->Branch("ind0_hits_charge_polar_merge", &ind0_hits_charge_polar_merge);
  fTree->Branch("ind1_hits_channel_polar_merge", &ind1_hits_channel_polar_merge);
  fTree->Branch("ind1_hits_time_polar_merge", &ind1_hits_time_polar_merge);
  fTree->Branch("ind1_hits_charge_polar_merge", &ind1_hits_charge_polar_merge);

}

void protoana::Pi0Shower::endJob()
{

}

void protoana::Pi0Shower::reset() 
{

  run = -999;
  subrun = -999;
  event = -999;

  /// Beam slice hits
  beam_vertex_hit_channel = -999;
  beam_vertex_hit_time = -999;
  selected_hits_channel.clear();
  selected_hits_time.clear();

  // Pre-cut
  coll_hits_channel_initial.clear(); coll_hits_time_initial.clear(); coll_hits_charge_initial.clear();
  ind0_hits_channel_initial.clear(); ind0_hits_time_initial.clear(); ind0_hits_charge_initial.clear();
  ind1_hits_channel_initial.clear(); ind1_hits_time_initial.clear(); ind1_hits_charge_initial.clear();

  // DBScan step
  coll_hit_pdg_dbscan.clear(); ind0_hit_pdg_dbscan.clear();  ind1_hit_pdg_dbscan.clear();
  coll_hits_channel_dbscan.clear(); coll_hits_time_dbscan.clear(); coll_hits_charge_dbscan.clear();
  ind0_hits_channel_dbscan.clear(); ind0_hits_time_dbscan.clear(); ind0_hits_charge_dbscan.clear();
  ind1_hits_channel_dbscan.clear(); ind1_hits_time_dbscan.clear(); ind1_hits_charge_dbscan.clear();

  // Polar cluster step
  coll_hit_pdg_polar_merge.clear(); ind0_hit_pdg_polar_merge.clear();  ind1_hit_pdg_polar_merge.clear();
  coll_hits_channel_polar_merge.clear(); coll_hits_time_polar_merge.clear(); coll_hits_charge_polar_merge.clear();  
  ind0_hits_channel_polar_merge.clear(); ind0_hits_time_polar_merge.clear(); ind0_hits_charge_polar_merge.clear();
  ind1_hits_channel_polar_merge.clear(); ind1_hits_time_polar_merge.clear(); ind1_hits_charge_polar_merge.clear();


}


DEFINE_ART_MODULE(protoana::Pi0Shower)

