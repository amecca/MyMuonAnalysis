// -*- C++ -*-
//
// Package:    MuonGeneralAnalyzer
// Class:      MuonGeneralAnalyzer
// 
/**\class MuonGeneralAnalyzer 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// System include files
#include <memory>
#include <iomanip>
#include <map>
#include <string>

// Framework and utilities include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Data formats 
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
// #include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/Math/interface/deltaPhi.h" 
#include "DataFormats/Math/interface/deltaR.h" 

// Tracking tools 
// #include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
// #include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"

// Muon geometry (some are probably not needed) 
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// Muon ID 

// Custom
#include "MyMuonAnalysis/DisplacedMuons/interface/CommonHitCounter.h"

// ROOT includes 
#include "TFile.h"
#include "TH1.h"
//#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TMath.h"

//
// class declaration
//

using namespace edm;
using std::vector;
using std::string;

#define MuonGeneralAnalyzer_BINS_pu 20,0.,120
#define MuonGeneralAnalyzer_BINS_allpu 20,200.,600

class MuonGeneralAnalyzer : public edm::one::EDAnalyzer<> {
  // Helper struct defining a "row" of associations
  struct AssociationRow {
    // edm::Ref<reco::MuonCollection>
    reco::MuonRef  muon; // muon with inner track
    reco::TrackRef inner;
    reco::TrackRef outer;
    reco::TrackRef global;
    reco::TrackRef oldSda;
    reco::TrackRef oldSdaUpd;
    reco::TrackRef oldGl;
    reco::TrackRef newSda;
    reco::TrackRef newSdaUpd;
    reco::TrackRef newGl;
    reco::GenParticleRef genMuon;
  };

  struct TrackInfo {
    template<class T=reco::TrackRef>
    TrackInfo(const T& t) :
      pt    (t->pt()   ),
      eta   (t->eta()  ),
      aeta  (fabs(eta) ),
      phi   (t->phi()  ),
      charge(t->charge())
    {}

    TrackInfo() = default;

    float pt, eta, aeta, phi, charge;
  };
  
public:
  explicit MuonGeneralAnalyzer(const edm::ParameterSet&);
  ~MuonGeneralAnalyzer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // Internal utilities
  void createPlots_old();
  void createPlots();
  void fillPlotsNumerator(const char* numerator, const char* denominator, const TrackInfo& den, const TrackInfo& num, int npuOOT, int npuIT);
  void fillPlots(const AssociationRow& row, int npuOOT, int npuIT);
  int fillMatchedPlots(const char* label, const reco::TrackCollection& trackCollection, const reco::GenParticle& gp, int npuOOT, int npuIT);
  void fillFakePlots(const char* label, const reco::TrackCollection& trackCollection, const std::vector<bool>& mask, int npuOOT, int npuIT);
  void endPlots(const char* label);
  template <class T> std::string debugRef(const edm::Ref<T>
					  // , size_t (*hashFn)(const Ref<T>&)
					  ) const;

  // ----------member data ---------------------------
  
  // Output file and histos
  std::string outputFileName;
  edm::Service<TFileService> outfile_;
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TGraphAsymmErrors*> graphs_;
  
  // Collections: name of the process (for debug)
  const std::string processTag_;
  
  // Tokens 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> gpToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::TrackCollection> oldSaTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> oldSaUpdTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> oldGlobalTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> newSaTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> newSaUpdTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> newGlobalTracksToken_;
  
  bool useGenMatchMap_;
  edm::EDGetTokenT<reco::GenParticleMatch> genParticleMatchToken_;
  edm::Handle<reco::GenParticleMatch> genParticleMatch_;

  CommonHitCounter commonHitCounter;

  // Debug options
  bool debugPrint;
  bool outputPrint;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonGeneralAnalyzer::MuonGeneralAnalyzer(const edm::ParameterSet& iConfig) :
  processTag_(""),  // (iConfig.getParameter<edm::InputTag>("saTracksTag").process()),
  puInfoTag_( consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag")) ),
  gpToken_( consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gpCollectionTag")) ),
  muonToken_(        consumes<reco::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonTag"       )) ),
  oldSaTracksToken_(    consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("oldSaTracksTag"   )) ),
  oldSaUpdTracksToken_( consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("oldSaUpdTracksTag")) ),
  oldGlobalTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("oldGlobalTracks"  )) ),
  newSaTracksToken_(    consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("newSaTracksTag"   )) ),
  newSaUpdTracksToken_( consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("newSaUpdTracksTag")) ),
  newGlobalTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("newGlobalTracks"  )) ),
  
  useGenMatchMap_(iConfig.existsAs<edm::InputTag>("genParticleMatch")),
  
  commonHitCounter(iConfig.getUntrackedParameterSet("hitCounterParams", edm::ParameterSet())),
  
  debugPrint     (iConfig.getUntrackedParameter<bool>("verbose"    , false)),
  outputPrint    (iConfig.getUntrackedParameter<bool>("printOutput", false))
{
  if(debugPrint) std::cout << "Inside Constructor" << std::endl;
  
  if(useGenMatchMap_)
    genParticleMatchToken_ = consumes<reco::GenParticleMatch>( iConfig.getParameter<edm::InputTag>("genParticleMatch") );
  
  // edm::ParameterSet psetCommonHitCounter;
  // if(iConfig.existsAs<edm::ParameterSet>("hitCounterParams"))
  //   psetCommonHitCounter = iConfig.getParameter<edm::ParameterSet>("hitCounterParams");
  // commonHitCounter = CommonHitCounter(psetCommonHitCounter);
  
  ////// == Plots ==
  // createPlots_old();
  createPlots();
}


MuonGeneralAnalyzer::~MuonGeneralAnalyzer() { 
  if(debugPrint) std::cout << "Inside Destructor" << std::endl;
}


//
// member functions
//
template <class T>
std::string MuonGeneralAnalyzer::debugRef(const edm::Ref<T> ref
					  // , size_t (*hashFn)(const Ref<T>&)
					  ) const{
  std::string result (Form(
			   "(process: %x, product: %x, key: %x)",
			   // hashFn(ref),
			   ref.id().processIndex(),
			   ref.id().productIndex(),
			   ref.key()
			   ));
  return result;
}


// ------------ method called for each event  ------------
void MuonGeneralAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using std::string;

  if(debugPrint || outputPrint) std::cout << "============================  Inside analyze " << processTag_ << ' '
					  // << iEvent.id().run() << "  "
					  // << iEvent.id().luminosityBlock() << "  "
					  // << iEvent.id().event() << "  "
					  << "============================" << std::endl;

  hists_["evt_muo_counter"]->Fill(1.);

  edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
  iEvent.getByToken(puInfoTag_,puInfoH);
  int npuOOT(0),npuIT(0),npuOOTm1(0);
  if(puInfoH.isValid()) {
    for(std::vector<PileupSummaryInfo>::const_iterator it=puInfoH->begin(); it!=puInfoH->end(); ++it) {
      if(it->getBunchCrossing()==0) {
	npuIT += it->getPU_NumInteractions();
      }
      else {
	npuOOT += it->getPU_NumInteractions();
	//npuOOT = 0; // TMP: figure out how to use this!!!
      }

      if(it->getBunchCrossing()<0)
	npuOOTm1 += it->getPU_NumInteractions();
    }
  }

  // GenParticles 
  edm::Handle<reco::GenParticleCollection> genParts; 
  iEvent.getByToken(gpToken_, genParts); 
  if(!genParts.isValid()) {
    throw cms::Exception("GenParticleCollection not valid!"); 
  } 
  else {
    if(debugPrint) std::cout << "Found GenParticleCollection! (" << genParts->size() << ')' << std::endl;
  }
  std::vector<reco::GenParticle> genMuons;
  genMuons.reserve(8); // Start from a reasonable size to skip some reallocations
  std::copy_if(genParts->begin(), genParts->end(), back_inserter(genMuons),
	       [](const reco::GenParticle gp) { return abs(gp.pdgId()) == 13; }
	       );
  if(debugPrint) std::cout << "of which (" << genMuons.size() << ") are gen muons" << std::endl;

  // Muons (all)
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  if(!muons.isValid()) {
    if(outputPrint || debugPrint) std::cout << "Muon collection not valid!" << std::endl;
    return;
  }
  else {
    if(debugPrint) std::cout << "Found Muon collection! (" << muons->size() << ')' << std::endl;
  }
  
  // GL tracks
  edm::Handle<reco::TrackCollection> oldGlTracks;
  iEvent.getByToken(oldGlobalTracksToken_, oldGlTracks);
  if(!oldGlTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "OLD GLOBAL track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found OLD GLOBAL track collection! (" << oldGlTracks->size() << ')' << std::endl;
  
  edm::Handle<reco::TrackCollection> newGlTracks;
  iEvent.getByToken(newGlobalTracksToken_, newGlTracks);
  if(!newGlTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "NEW GLOBAL track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found NEW GLOBAL track collection! (" << newGlTracks->size() << ')' << std::endl;

  // SA tracks
  edm::Handle<reco::TrackCollection> oldSaTracks;
  iEvent.getByToken(oldSaTracksToken_, oldSaTracks);
  if(!oldSaTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "OLD SA track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found OLD SA track collection! (" << oldSaTracks->size() << ')' << std::endl;
  
  edm::Handle<reco::TrackCollection> newSaTracks;
  iEvent.getByToken(newSaTracksToken_, newSaTracks);
  if(!newSaTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "NEW SA track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found NEW SA track collection! (" << newSaTracks->size() << ')' << std::endl;

  // SA-updated at vertex tracks
  edm::Handle<reco::TrackCollection> oldSaUpdTracks;
  iEvent.getByToken(oldSaUpdTracksToken_, oldSaUpdTracks);
  if(!oldSaUpdTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "SA Updated track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found SA Updated track collection! (" << oldSaUpdTracks->size() << ')' << std::endl;

  edm::Handle<reco::TrackCollection> newSaUpdTracks;
  iEvent.getByToken(newSaUpdTracksToken_, newSaUpdTracks);
  if(!newSaUpdTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "SA Updated track collection not valid!" << std::endl;
    return;
  }
  else if(debugPrint) std::cout << "Found SA Updated track collection! (" << newSaUpdTracks->size() << ')' << std::endl;


  // ################################################################################
  // // test
  // for(auto sa : *saTracks){
  //   for(auto saUpd : *saUpdTracks){
  //     auto tuple_hits = commonHitCounter.countMatchingHits(sa, saUpd, true);
  //     std::cout << "tuple_hits: " << std::get<0>(tuple_hits) << ' ' << std::get<1>(tuple_hits) << ' ' << std::get<2>(tuple_hits) << '\n';
  //   }
  // }
  // return;
  // ################################################################################

  //std::cout << "--- 1" << std::endl;
  reco::TrackCollection oldSa0hitsTracks;
  for(auto trk : *oldSaTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      oldSa0hitsTracks.push_back(trk);
  }

  reco::TrackCollection oldSaUpd0hitsTracks;
  for(auto trk : *oldSaUpdTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      oldSaUpd0hitsTracks.push_back(trk);
  }
  
  reco::TrackCollection gl0hitsTracks;
  for(auto trk : *oldGlTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      gl0hitsTracks.push_back(trk);
  }
  
  std::vector<bool>    oldSaMask(   oldSaTracks->size(), false);
  std::vector<bool> oldSaUpdMask(oldSaUpdTracks->size(), false);
  std::vector<bool>    oldGlMask(   oldGlTracks->size(), false);
  
  //std::cout << "--- 2" << std::endl;

  // ********************************************************************************
  // Collections of trackrefs of the muons: match FROM
  std::vector<reco::TrackRef> innerTrackRefs;
  std::vector<reco::TrackRef> outerTrackRefs;
  std::vector<reco::TrackRef> globalTrackRefs;
  innerTrackRefs .reserve(muons->size());
  outerTrackRefs .reserve(muons->size());
  globalTrackRefs.reserve(muons->size());
  for(auto muon : *muons){
    reco::TrackRef innerTrackRef  = muon.innerTrack ();
    reco::TrackRef outerTrackRef  = muon.outerTrack ();
    reco::TrackRef globalTrackRef = muon.globalTrack();
    if(innerTrackRef.isAvailable()  && innerTrackRef.isNonnull() )
      innerTrackRefs.push_back( std::move(innerTrackRef ));
    if(outerTrackRef.isAvailable()  && outerTrackRef.isNonnull() )
      outerTrackRefs.push_back( std::move(outerTrackRef ));
    if(globalTrackRef.isAvailable() && globalTrackRef.isNonnull())
      globalTrackRefs.push_back(std::move(globalTrackRef));
  }

  CommonHitCounter::map_type outerToOldSdA     = commonHitCounter.matchingTrackCollections(outerTrackRefs , *oldSaTracks   );
  CommonHitCounter::map_type outerToOldSdAUpd  = commonHitCounter.matchingTrackCollections(outerTrackRefs , *oldSaUpdTracks);
  CommonHitCounter::map_type globalToOldGl     = commonHitCounter.matchingTrackCollections(globalTrackRefs, *oldGlTracks   );
  // CommonHitCounter::map_type innerToOldGl      = commonHitCounter.matchingTrackCollections(innerTrackRefs , *oldGlTracks   );

  // if(debugPrint){
  //   std::cout << ">>> Dumping outer>SdA map:\n";
  //   for(auto it : outerToSdA)
  //     std::cout << '\t' 
  // 		<< Form("%lx", outerToSdA.hash_function()(it.first)) << ' ' << debugRef(it.first ) << "  " 
  // 		<< Form("%lx", outerToSdA.hash_function()(it.first)) << ' ' << debugRef(it.second)
  // 		<< '\n';    
  // }

  // Contruct the table 
  if(debugPrint) std::cout << ">>> building AssociationTable:\n";
  std::vector<AssociationRow> associationTable;
  associationTable.reserve(muons->size());
  for(unsigned int i = 0; i < muons->size(); ++i){
    AssociationRow row;
    row.muon = reco::MuonRef(muons, i);
    row.inner  = muons->at(i).innerTrack();
    row.outer  = muons->at(i).outerTrack();
    row.global = muons->at(i).globalTrack();

    if(debugPrint)
      std::cout << "\tMuon "<<i<<":  itself: " << debugRef(row.muon);

    if(row.outer.isNonnull()){
      if(debugPrint)
      	std::cout << "  outer: " << Form("%lx", outerToOldSdA.hash_function()(row.outer)) << debugRef(row.outer);
      
      auto it_oldSdA    = outerToOldSdA   .find(row.outer);
      auto it_oldSdAUpd = outerToOldSdAUpd.find(row.outer);
      
      if(debugPrint){
      	std::cout<<"  SdA = ";
      	  if(it_oldSdA    == outerToOldSdA.end())
      	    std::cout << "END";
      	  else{
      	    std::cout << std::distance(outerToOldSdA.begin(), it_oldSdA) << " (";
      	    if(it_oldSdA->second.isNonnull()) std::cout << it_oldSdA->second.key();
      	    else                              std::cout << "NULL";
      	    std::cout << ')';
      	  }
      }
      
      row.oldSda    = it_oldSdA    != outerToOldSdA.end()    ? it_oldSdA   ->second : reco::TrackRef(oldSaTracks.id()   );
      row.oldSdaUpd = it_oldSdAUpd != outerToOldSdAUpd.end() ? it_oldSdAUpd->second : reco::TrackRef(oldSaUpdTracks.id());
    }
    else if(debugPrint) std::cout << "  outer is null!";
    if(debugPrint) std::cout << '\n';
    
    // Global
    if(row.muon->globalTrack().isNonnull()){
      CommonHitCounter::map_type::const_iterator it_oldGl;
      if     ( (it_oldGl = globalToOldGl.find(row.muon->globalTrack())) != globalToOldGl.end() )
	row.oldGl = it_oldGl->second;  // First try matching global tracks
      // else if( (it_oldGl = innerToOldGl .find(row.muon->globalTrack())) != innerToOldGl .end() )
      // 	row.oldGl = it_oldGl->second;  // Then try matching the inner part of the global
      else
	row.oldGl = reco::TrackRef(oldGlTracks.id());
    }
    
    // Gen muon
    float muEta = row.muon->eta();
    float muPhi = row.muon->phi();
    auto it_closestGen = std::min_element(genMuons.begin(), genMuons.end(), [muEta, muPhi](const reco::GenParticle& a, const reco::GenParticle& b) {
	return deltaR(muEta, muPhi, a.eta(), a.phi()) < deltaR(muEta, muPhi, b.eta(), b.phi());
      } );
    if(it_closestGen != genMuons.end() && deltaR(muEta, muPhi, it_closestGen->eta(), it_closestGen->phi()) < 0.2)
      row.genMuon = reco::GenParticleRef(&genMuons, std::distance(genMuons.begin(), it_closestGen));
    else
      row.genMuon = reco::GenParticleRef();

    associationTable.push_back(std::move(row));
  }
  
  
  std::cout << Form("%6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %s\n", "muon", "genMu", "inner", "outer", "global", "oldSA", "newSA", "oldSaU", "newSaU", "oldGl", "newGl");
  for(const AssociationRow& row : associationTable)
    std::cout << Form("%6d | %6d | %6s | %6s | %6s | %6d | %6d | %6d | %6d | %6d | %6d \n",
  		      row.muon.key(),
  		      row.genMuon.isNonnull() ? row.genMuon.key() : -1,
		      row.inner .isNonnull() ? "OK" : "-",
  		      row.outer .isNonnull() ? "OK" : "-",
		      row.global.isNonnull() ? "OK" : "-",
  		      row.oldSda   .key(), // row.sdaTrack   .isNonnull() ? row.sdaTrack   .key() : -1,
  		      row.newSda   .key(),
		      row.oldSdaUpd.key(),
		      row.newSdaUpd.key(),
		      row.oldGl    .key(),
		      row.newGl    .key()
  		      );

  std::cout<<"--------------------------------------------------------------------------------\n";

  for(const AssociationRow& row : associationTable){
    fillPlots(row, npuOOT, npuIT);
  }

  return;

  // ********************************************************************************
  
  reco::TrackCollection genMuonTracks;
  // Fill plots
  for(auto gp : *genParts) {
    if(abs(gp.pdgId())!=13) continue;
    // Denominator
    float geneta  = gp.eta();
    float genaeta = fabs(geneta);
    if(genaeta>2.4) continue;
    hists_["evt_muo_counter"]->Fill(2.);
    float genpt = gp.pt();
    float genphi = gp.phi();
    
    hists_["eff_den_aeta" ]->Fill(genaeta);
    hists_["eff_den_pt"   ]->Fill(genpt);
    hists_["eff_den_phi"  ]->Fill(genphi);
    hists_["eff_den_pu"   ]->Fill(npuIT);
    hists_["eff_den_allpu"]->Fill(npuIT+npuOOT);

    //std::cout << "--- 2.1" << std::endl;
    
    // SA tracks
    int saIdx = fillMatchedPlots("sa"        , *oldSaTracks       , gp, npuOOT, npuIT);
    if(saIdx >= 0)
      oldSaMask.at(saIdx) = true;
    // std::cout << "sa: "<<saIdx<<std::endl;
    // fillMatchedPlots("sa0hits"   , sa0hitsTracks  , gp, npuOOT, npuIT);
    
    // SA-updated tracks
    int saupdIdx = fillMatchedPlots("saupd"  , *oldSaUpdTracks    , gp, npuOOT, npuIT);
    if(saupdIdx >= 0)
      oldSaUpdMask.at(saupdIdx) = true;
    // std::cout << "saupd: "<<saupdIdx<<std::endl;
    // fillMatchedPlots("saupd0hits", saupd0hitsTracks, gp, npuOOT, npuIT);

    // GLOBAL tracks
    int glIdx = fillMatchedPlots("gl"        , *oldGlTracks       , gp, npuOOT, npuIT);
    if(glIdx >= 0)
      oldGlMask.at(glIdx) = true;
    // std::cout << "gl: "<<glIdx<<std::endl;
    // fillMatchedPlots("gl0hits"   , gl0hitsTracks   , gp, npuOOT, npuIT);
  }
  
  fillFakePlots("sa"   , *oldSaTracks   , oldSaMask   , npuOOT, npuIT);
  fillFakePlots("saupd", *oldSaUpdTracks, oldSaUpdMask, npuOOT, npuIT);
  fillFakePlots("gl"   , *oldGlTracks   , oldGlMask   , npuOOT, npuIT);

  return;
}


// ------------ method called once each job just before starting event loop  ------------
void MuonGeneralAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonGeneralAnalyzer::endJob() {
  endPlots("sa");
  endPlots("sa0hits");
  endPlots("saupd");
  endPlots("saupd0hits");
  endPlots("gl");
  endPlots("gl0hits");
} 

// ------------ method called when starting to processes a run  ------------
void MuonGeneralAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // edm::ESHandle<TrackAssociatorBase> theAssociator;
  // for(size_t w=0; w<associatorNames.size(); ++w) { 
  //   iSetup.get<TrackAssociatorRecord>().get(associatorNames[w], theAssociator); 
  //   if(theAssociator.isValid()) associators.push_back(theAssociator.product()); 
  // } 
  // if(associators.size()==0) { 
  //   TString assonames; 
  //   for(size_t w=0; w<associatorNames.size(); ++w) assonames += (associatorNames[w]+" "); 
  //   throw cms::Exception("No Associators found")
  //     << "Cannot find any of the following associators: " << assonames.Data(); 
  // } 
} // end of beginRun


// ------------ method called when ending the processing of a run  ------------
void MuonGeneralAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}


// ------------ method called when starting to processes a luminosity block  ------------
void MuonGeneralAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuonGeneralAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {
}


void MuonGeneralAnalyzer::createPlots_old(){
  // -- Binnings --
  //double pt_bins[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 150., 200.};
  double pt_bins[] = {0., 1., 2., 3., 4., 5., 7.5, 10., 15., 20., 35., 50., 100., 200.};
  int n_pt_bins = sizeof(pt_bins)/sizeof(*pt_bins);

  // //double eta_bins[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};
  // double eta_bins[] = {-3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
  // int n_eta_bins = sizeof(eta_bins)/sizeof(double);

  double aeta_bins[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6};
  int n_aeta_bins = sizeof(aeta_bins)/sizeof(*aeta_bins);

  static const double pig=3.1416;

  //// ~~ Common ~~
  hists_["eff_den_pt"      ] = outfile_->make<TH1F>("eff_den_pt"    , ";muon p_{T} [GeV];events (denominator)"                      , n_pt_bins-1  , pt_bins  );
  hists_["eff_den_aeta"    ] = outfile_->make<TH1F>("eff_den_aeta"  , ";muon |#eta|;events (denominator)"                           , n_aeta_bins-1, aeta_bins);
  hists_["eff_den_phi"     ] = outfile_->make<TH1F>("eff_den_phi"   , ";muon #phi [rad];events (denominator)"                       ,  10, -pig    , pig      );
  hists_["eff_den_pu"      ] = outfile_->make<TH1F>("eff_den_pu"    , ";number of in-time vertices;events (denominator)"            , MuonGeneralAnalyzer_BINS_pu   );
  hists_["eff_den_allpu"   ] = outfile_->make<TH1F>("eff_den_allpu" , ";number of in- and out-of-time vertices;events (denominator)", MuonGeneralAnalyzer_BINS_allpu);
  
  // Event and segment counter 
  hists_["evt_muo_counter"] = outfile_->make<TH1F>("evt_muo_counter", ";;counts", 2, 0.5, 2.5); 
  hists_["evt_muo_counter"]->GetXaxis()->SetBinLabel(1, "N. events"); 
  hists_["evt_muo_counter"]->GetXaxis()->SetBinLabel(2, "N. muons"); 
  
  for(const char* label : {"sa", "saupd", "gl", "sa0hits", "saupd0hits", "gl0hits"}){
    // -- Resolutions --
    hists_[Form("%s_res_pt" , label)] = outfile_->make<TH1F>(Form("%s_res_pt" , label), ";#Delta(q/p_{T})_{REC-GEN}/(q/p_{T})_{GEN};events", 100, -5., 5.);
    hists_[Form("%s_res_eta", label)] = outfile_->make<TH1F>(Form("%s_res_eta", label), ";#Delta#eta_{REC-GEN};events"                     , 60, -.3, .3);
    hists_[Form("%s_res_phi", label)] = outfile_->make<TH1F>(Form("%s_res_phi", label), ";#Delta#phi_{REC-GEN};events"                     , 60, -.3, .3);
    hists_[Form("%s_res_dr" , label)] = outfile_->make<TH1F>(Form("%s_res_dr" , label), ";#DeltaR_{REC-GEN};events"                       , 60,  0., .3);
      
    // -- Track quality --
    hists_[Form("%s_nhits"       , label)] = outfile_->make<TH1F>(Form("%s_nhits"       , label), ";number of muon hits;events"        , 60, 0., 60.);
    hists_[Form("%s_nhits_bar"   , label)] = outfile_->make<TH1F>(Form("%s_nhits_bar"   , label), ";number of barrel muon hits;events" , 60, 0., 60.);
    hists_[Form("%s_nhits_ovl"   , label)] = outfile_->make<TH1F>(Form("%s_nhits_ovl"   , label), ";number of overlap muon hits;events", 60, 0., 60.);
    hists_[Form("%s_nhits_end"   , label)] = outfile_->make<TH1F>(Form("%s_nhits_end"   , label), ";number of endcap muon hits;events" , 60, 0., 60.);
    hists_[Form("%s_nvalhits"    , label)] = outfile_->make<TH1F>(Form("%s_nvalhits"    , label), ";number of valid muon hits;events"        , 60, 0., 60.);
    hists_[Form("%s_nvalhits_bar", label)] = outfile_->make<TH1F>(Form("%s_nvalhits_bar", label), ";number of valid barrel muon hits;events" , 60, 0., 60.);
    hists_[Form("%s_nvalhits_ovl", label)] = outfile_->make<TH1F>(Form("%s_nvalhits_ovl", label), ";number of valid overlap muon hits;events", 60, 0., 60.);
    hists_[Form("%s_nvalhits_end", label)] = outfile_->make<TH1F>(Form("%s_nvalhits_end", label), ";number of valid endcap muon hits;events" , 60, 0., 60.);
    hists_[Form("%s_nsegs"       , label)] = outfile_->make<TH1F>(Form("%s_nsegs"       , label), ";number of segments;events"        ,  5, 0.,  5.);
    hists_[Form("%s_nsegs_bar"   , label)] = outfile_->make<TH1F>(Form("%s_nsegs_bar"   , label), ";number of barrel segments;events" ,  5, 0.,  5.);
    hists_[Form("%s_nsegs_ovl"   , label)] = outfile_->make<TH1F>(Form("%s_nsegs_ovl"   , label), ";number of overlap segments;events",  5, 0.,  5.);
    hists_[Form("%s_nsegs_end"   , label)] = outfile_->make<TH1F>(Form("%s_nsegs_end"   , label), ";number of endcap segments;events" ,  5, 0.,  5.);
    hists_[Form("%s_nvalsegs"    , label)] = outfile_->make<TH1F>(Form("%s_nvalsegs"    , label), ";number of valid segments;events"        ,  5, 0.,  5.);
    hists_[Form("%s_nvalsegs_bar", label)] = outfile_->make<TH1F>(Form("%s_nvalsegs_bar", label), ";number of valid barrel segments;events" ,  5, 0.,  5.);
    hists_[Form("%s_nvalsegs_ovl", label)] = outfile_->make<TH1F>(Form("%s_nvalsegs_ovl", label), ";number of valid overlap segments;events",  5, 0.,  5.);
    hists_[Form("%s_nvalsegs_end", label)] = outfile_->make<TH1F>(Form("%s_nvalsegs_end", label), ";number of valid endcap segments;events" ,  5, 0.,  5.);

    // -- Efficiencies -- 
    hists_[Form("%s_eff_num_pt"   , label)] = outfile_->make<TH1F>(Form("%s_eff_num_pt"   , label), ";muon p_{T} [GeV];events (numerator)"                      , n_pt_bins-1  , pt_bins  );
    hists_[Form("%s_eff_num_aeta" , label)] = outfile_->make<TH1F>(Form("%s_eff_num_aeta" , label), ";muon |#eta|;events (numerator)"                           , n_aeta_bins-1, aeta_bins);
    hists_[Form("%s_eff_num_phi"  , label)] = outfile_->make<TH1F>(Form("%s_eff_num_phi"  , label), ";muon #phi [rad];events (numerator)"                       ,  10, -pig    , pig      );
    hists_[Form("%s_eff_num_pu"   , label)] = outfile_->make<TH1F>(Form("%s_eff_num_pu "  , label), ";number of in-time vertices;events (numerator)"            ,MuonGeneralAnalyzer_BINS_pu   );
    hists_[Form("%s_eff_num_allpu", label)] = outfile_->make<TH1F>(Form("%s_eff_num_allpu", label), ";number of in- and out-of-time vertices;events (numerator)",MuonGeneralAnalyzer_BINS_allpu);
    
    graphs_[Form("%s_eff_pt_err"   , label)] = outfile_->make<TGraphAsymmErrors>();
    graphs_[Form("%s_eff_aeta_err" , label)] = outfile_->make<TGraphAsymmErrors>();
    graphs_[Form("%s_eff_phi_err"  , label)] = outfile_->make<TGraphAsymmErrors>();
    graphs_[Form("%s_eff_pu_err"   , label)] = outfile_->make<TGraphAsymmErrors>();
    graphs_[Form("%s_eff_allpu_err", label)] = outfile_->make<TGraphAsymmErrors>();
    
    graphs_[Form("%s_eff_pt_err"   , label)]->SetNameTitle(Form("%s_eff_pt_err"   , label), ";muon p_{T} [GeV];efficiency"                      );
    graphs_[Form("%s_eff_aeta_err" , label)]->SetNameTitle(Form("%s_eff_aeta_err" , label), ";muon |#eta|;efficiency"                           );
    graphs_[Form("%s_eff_phi_err"  , label)]->SetNameTitle(Form("%s_eff_phi_err"  , label), ";muon #phi [rad];efficiency"                       );
    graphs_[Form("%s_eff_pu_err"   , label)]->SetNameTitle(Form("%s_eff_pu_err"   , label), ";number of in-time vertices;efficiency"            );
    graphs_[Form("%s_eff_allpu_err", label)]->SetNameTitle(Form("%s_eff_allpu_err", label), ";number of in- and out-of-time vertices;efficiency");
  
    // -- Fakes --
    hists_[Form("%s_fake_pt_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_pt_perevt"   , label), ";muon p_{T} [GeV];n. fakes/event"                      , n_pt_bins-1  , pt_bins  );
    hists_[Form("%s_fake_aeta_perevt" , label)] = outfile_->make<TH1F>(Form("%s_fake_aeta_perevt" , label), ";muon |#eta|;n. fakes/event"                           , n_aeta_bins-1, aeta_bins);
    hists_[Form("%s_fake_phi_perevt"  , label)] = outfile_->make<TH1F>(Form("%s_fake_phi_perevt"  , label), ";muon #phi [rad];n. fakes/event"                       ,  10    , -pig, pig      );
    hists_[Form("%s_fake_pu_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_pu_perevt"   , label), ";number of in-time vertices;n. fakes/event"            , MuonGeneralAnalyzer_BINS_pu   );
    hists_[Form("%s_fake_allpu_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_allpu_perevt", label), ";number of in- and out-of-time vertices;n. fakes/event", MuonGeneralAnalyzer_BINS_allpu);
  
    hists_[Form("%s_fake_nhits_perevt"       , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_perevt"       , label), ";number of muon hits;n. fakes/event"            , 60, 0., 60.);
    hists_[Form("%s_fake_nhits_bar_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_bar_perevt"   , label), ";number of barrel muon hits;n. fakes/event"     , 60, 0., 60.);
    hists_[Form("%s_fake_nhits_ovl_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_ovl_perevt"   , label), ";number of overlap muon hits;n. fakes/event"    , 60, 0., 60.);
    hists_[Form("%s_fake_nhits_end_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_end_perevt"   , label), ";number of endcap muon hits;n. fakes/event"     , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_perevt"    , label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_perevt"    , label), ";number of valid muon hits;n. fakes/event"        , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_bar_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_bar_perevt", label), ";number of valid barrel muon hits;n. fakes/event" , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_ovl_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_ovl_perevt", label), ";number of valid overlap muon hits;n. fakes/event", 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_end_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_end_perevt", label), ";number of valid endcap muon hits;n. fakes/event" , 60, 0., 60.);
    hists_[Form("%s_fake_nsegs_perevt"       , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_perevt"       , label), ";number of segments;n. fakes/event"        ,  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_bar_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_bar_perevt"   , label), ";number of barrel segments;n. fakes/event" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_ovl_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_ovl_perevt"   , label), ";number of overlap segments;n. fakes/event",  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_end_perevt"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_end_perevt"   , label), ";number of endcap segments;n. fakes/event" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_perevt"    , label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_perevt"    , label), ";number of valid segments;n. fakes/event"        ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_bar_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_bar_perevt", label), ";number of valid barrel segments;n. fakes/event" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_ovl_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_ovl_perevt", label), ";number of valid overlap segments;n. fakes/event",  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_end_perevt", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_end_perevt", label), ";number of valid endcap segments;n. fakes/event" ,  5, 0.,  5.);
  
    hists_[Form("%s_fake_pt_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_pt_permuo"   , label), ";muon p_{T} [GeV];n. fakes/muon"                      , n_pt_bins-1  , pt_bins  );
    hists_[Form("%s_fake_aeta_permuo" , label)] = outfile_->make<TH1F>(Form("%s_fake_aeta_permuo" , label), ";muon |#eta|;n. fakes/muon"                           , n_aeta_bins-1, aeta_bins);
    hists_[Form("%s_fake_phi_permuo"  , label)] = outfile_->make<TH1F>(Form("%s_fake_phi_permuo"  , label), ";muon #phi [rad];n. fakes/muon"                       ,  10    , -pig, pig      );
    hists_[Form("%s_fake_pu_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_pu_permuo"   , label), ";number of in-time vertices;n. fakes/muon"            , MuonGeneralAnalyzer_BINS_pu   );
    hists_[Form("%s_fake_allpu_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_allpu_permuo", label), ";number of in- and out-of-time vertices;n. fakes/muon", MuonGeneralAnalyzer_BINS_allpu); 
	   
    hists_[Form("%s_fake_nhits_permuo"       , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_permuo"       , label), ";number of muon hits;n. fakes/muon"        , 60, 0., 60.);
    hists_[Form("%s_fake_nhits_bar_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_bar_permuo"   , label), ";number of barrel muon hits;n. fakes/muon" , 60, 0., 60.);
    hists_[Form("%s_fake_nhits_ovl_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_ovl_permuo"   , label), ";number of overlap muon hits;n. fakes/muon", 60, 0., 60.);
    hists_[Form("%s_fake_nhits_end_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nhits_end_permuo"   , label), ";number of endcap muon hits;n. fakes/muon" , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_permuo"    , label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_permuo"    , label), ";number of valid muon hits;n. fakes/muon"        , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_bar_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_bar_permuo", label), ";number of valid barrel muon hits;n. fakes/muon" , 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_ovl_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_ovl_permuo", label), ";number of valid overlap muon hits;n. fakes/muon", 60, 0., 60.);
    hists_[Form("%s_fake_nvalhits_end_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalhits_end_permuo", label), ";number of valid endcap muon hits;n. fakes/muon" , 60, 0., 60.);
    hists_[Form("%s_fake_nsegs_permuo"       , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_permuo"       , label), ";number of segments;n. fakes/muon"        ,  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_bar_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_bar_permuo"   , label), ";number of barrel segments;n. fakes/muon" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_ovl_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_ovl_permuo"   , label), ";number of overlap segments;n. fakes/muon",  5, 0.,  5.);
    hists_[Form("%s_fake_nsegs_end_permuo"   , label)] = outfile_->make<TH1F>(Form("%s_fake_nsegs_end_permuo"   , label), ";number of endcap segments;n. fakes/muon" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_permuo"    , label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_permuo"    , label), ";number of valid segments;n. fakes/muon"        ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_bar_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_bar_permuo", label), ";number of valid barrel segments;n. fakes/muon" ,  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_ovl_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_ovl_permuo", label), ";number of valid overlap segments;n. fakes/muon",  5, 0.,  5.);
    hists_[Form("%s_fake_nvalsegs_end_permuo", label)] = outfile_->make<TH1F>(Form("%s_fake_nvalsegs_end_permuo", label), ";number of valid endcap segments;n. fakes/muon" ,  5, 0.,  5.);
  }
}


void MuonGeneralAnalyzer::createPlots(){
  std::cout << ">>> CREATING PLOTS\n";
  // -- Binnings --
  //double pt_bins[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 150., 200.};
  double pt_bins[] = {0., 1., 2., 3., 4., 5., 7.5, 10., 15., 20., 35., 50., 100., 200.};
  int n_pt_bins = sizeof(pt_bins)/sizeof(*pt_bins);

  double aeta_bins[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6};
  int n_aeta_bins = sizeof(aeta_bins)/sizeof(*aeta_bins);

  static const double pig = TMath::Pi();

  // Event and segment counter
  hists_["evt_muo_counter"] = outfile_->make<TH1F>("evt_muo_counter", ";;counts", 2, 0.5, 2.5);
  hists_["evt_muo_counter"]->GetXaxis()->SetBinLabel(1, "N. events");
  hists_["evt_muo_counter"]->GetXaxis()->SetBinLabel(2, "N. muons");

  //// ~~ Common ~~
  // plot category   :   all  | gen_matched | no_match
  // reference track :  inner |     gen     |  inner

  // Efficiency
  for(const char* denominator : {"inner", "gen"}){
    std::cout << '\t' << denominator << ' ';
    std::string name_den_pt   ( Form("eff_den_%s_pt"   , denominator) );
    std::string name_den_aeta ( Form("eff_den_%s_aeta" , denominator) );
    std::string name_den_phi  ( Form("eff_den_%s_phi"  , denominator) );
    std::string name_den_pu   ( Form("eff_den_%s_pu"   , denominator) );
    std::string name_den_allpu( Form("eff_den_%s_allpu", denominator) );
    std::cout << "\tcreating " << name_den_pt << std::endl;
    hists_[name_den_pt   .c_str()] = outfile_->make<TH1F>(name_den_pt   .c_str(), Form(";%s track p_{T} [GeV];events (denominator)", denominator), n_pt_bins-1  , pt_bins  );
    hists_[name_den_aeta .c_str()] = outfile_->make<TH1F>(name_den_aeta .c_str(), Form(";%s track |#eta|;events (denominator)"     , denominator), n_aeta_bins-1, aeta_bins);
    hists_[name_den_phi  .c_str()] = outfile_->make<TH1F>(name_den_phi  .c_str(), Form(";%s track #phi [rad];events (denominator)" , denominator),  10, -pig    , pig      );
    hists_[name_den_pu   .c_str()] = outfile_->make<TH1F>(name_den_pu   .c_str(),      ";# of IT vertices;events (denominator)"    , MuonGeneralAnalyzer_BINS_pu   );
    hists_[name_den_allpu.c_str()] = outfile_->make<TH1F>(name_den_allpu.c_str(),      ";# of IT+OOT vertices;events (denominator)", MuonGeneralAnalyzer_BINS_allpu);

    for(const char* numerator: {"outer", "global", "oldSda", "oldSdaUpd", "oldGl", "newSda", "newSdaUpd", "newGl"}){
      std::cout << "\t\t" << numerator << ' ';
      std::string name_pt   ( Form("eff_num_%s_%s_pt"   , denominator, numerator) );
      std::cout << "\t\tcreating " << name_pt << std::endl;
      std::string name_aeta ( Form("eff_num_%s_%s_aeta" , denominator, numerator) );
      std::string name_phi  ( Form("eff_num_%s_%s_phi"  , denominator, numerator) );
      std::string name_pu   ( Form("eff_num_%s_%s_pu"   , denominator, numerator) );
      std::string name_allpu( Form("eff_num_%s_%s_allpu", denominator, numerator) );
      std::cout << name_pt << '\n';
      hists_[name_pt   .c_str()] = outfile_->make<TH1F>(name_pt   .c_str(), Form(";%s track p_{T} [GeV];events (numerator)", denominator), n_pt_bins-1  , pt_bins  );
      hists_[name_aeta .c_str()] = outfile_->make<TH1F>(name_aeta .c_str(), Form(";%s track |#eta|;events (numerator)"     , denominator), n_aeta_bins-1, aeta_bins);
      hists_[name_phi  .c_str()] = outfile_->make<TH1F>(name_phi  .c_str(), Form(";%s track #phi [rad];events (numerator)" , denominator),  10, -pig    , pig      );
      hists_[name_pu   .c_str()] = outfile_->make<TH1F>(name_pu   .c_str(),      ";# of IT vertices;events (numerator)"    , MuonGeneralAnalyzer_BINS_pu   );
      hists_[name_allpu.c_str()] = outfile_->make<TH1F>(name_allpu.c_str(),      ";# of IT+OOT vertices;events (numerator)", MuonGeneralAnalyzer_BINS_allpu);

      std::string res_pt  = Form("res_%s_%s_pt" , numerator, denominator);
      std::string res_eta = Form("res_%s_%s_eta", numerator, denominator);
      std::string res_phi = Form("res_%s_%s_phi", numerator, denominator);
      std::string res_dr  = Form("res_%s_%s_dr" , numerator, denominator);
      hists_[res_pt .c_str()] = outfile_->make<TH1F>(res_pt .c_str(), ";#Delta(q/p_{T})_{REC-GEN}/(q/p_{T})_{GEN};events", 100, -5., 5.);
      hists_[res_eta.c_str()] = outfile_->make<TH1F>(res_eta.c_str(), ";#Delta#eta_{REC-GEN};events"                     , 60, -.3, .3);
      hists_[res_phi.c_str()] = outfile_->make<TH1F>(res_phi.c_str(), ";#Delta#phi_{REC-GEN};events"                     , 60, -.3, .3);
      hists_[res_dr .c_str()] = outfile_->make<TH1F>(res_dr .c_str(), ";#DeltaR_{REC-GEN};events"                        , 60,  0., .3);
    }
  }
  std::cout << "... DONE !\n";
}


void MuonGeneralAnalyzer::fillPlotsNumerator(const char* numerator, const char* denominator, const MuonGeneralAnalyzer::TrackInfo& den, const MuonGeneralAnalyzer::TrackInfo& num, int npuOOT, int npuIT){
  hists_[ Form("eff_num_%s_%s_pt"   , denominator, numerator) ]->Fill(den.pt  );
  hists_[ Form("eff_num_%s_%s_aeta" , denominator, numerator) ]->Fill(den.aeta);
  hists_[ Form("eff_num_%s_%s_phi"  , denominator, numerator) ]->Fill(den.phi );
  hists_[ Form("eff_num_%s_%s_pu"   , denominator, numerator) ]->Fill(npuIT   );
  hists_[ Form("eff_num_%s_%s_allpu", denominator, numerator) ]->Fill(npuIT+npuOOT);

  // if     (num.pt == 0)
  //   std::cout << "ERROR: " << numerator   << " pt is null!" << std::endl;
  // else if(den.pt == 0)
  //   std::cout << "ERROR: " << denominator << " pt is null!" << std::endl;
  // else
  hists_[Form("res_%s_%s_pt" , numerator, denominator)]->Fill((num.charge/den.charge)*(den.pt/num.pt) - 1.);
  hists_[Form("res_%s_%s_eta", numerator, denominator)]->Fill(num.eta - den.eta);
  hists_[Form("res_%s_%s_phi", numerator, denominator)]->Fill(deltaPhi(num.phi, den.phi));
  hists_[Form("res_%s_%s_dr" , numerator, denominator)]->Fill(deltaR(den.eta, den.phi, num.eta, num.phi));
}


void MuonGeneralAnalyzer::fillPlots(const MuonGeneralAnalyzer::AssociationRow& row, int npuOOT, int npuIT){
  char denominator[8];
  TrackInfo den;
  if(row.genMuon.isNonnull()){
    den = TrackInfo(row.genMuon);
    snprintf(denominator, 8, "gen");
  }
  else if(row.inner.isNonnull()){
    den = TrackInfo(row.inner);
    snprintf(denominator, 8, "inner");
  }
  else{
    std::cout << "WARN: null inner track!";
    return;
  }

  hists_[ Form("eff_den_%s_pt"   , denominator) ]->Fill(den.pt  );
  hists_[ Form("eff_den_%s_aeta" , denominator) ]->Fill(den.aeta);
  hists_[ Form("eff_den_%s_phi"  , denominator) ]->Fill(den.phi );
  hists_[ Form("eff_den_%s_pu"   , denominator) ]->Fill(npuIT   );
  hists_[ Form("eff_den_%s_allpu", denominator) ]->Fill(npuIT+npuOOT);

  if(row.outer.isNonnull()    ) fillPlotsNumerator("outer"    , denominator, den, TrackInfo(row.outer    ), npuOOT, npuIT);
  if(row.global.isNonnull()   ) fillPlotsNumerator("global"   , denominator, den, TrackInfo(row.global   ), npuOOT, npuIT);
  if(row.oldSda.isNonnull()   ) fillPlotsNumerator("oldSda"   , denominator, den, TrackInfo(row.oldSda   ), npuOOT, npuIT);
  if(row.oldSdaUpd.isNonnull()) fillPlotsNumerator("oldSdaUpd", denominator, den, TrackInfo(row.oldSdaUpd), npuOOT, npuIT);
  if(row.oldGl.isNonnull()    ) fillPlotsNumerator("oldGl"    , denominator, den, TrackInfo(row.oldGl    ), npuOOT, npuIT);
  if(row.newSda.isNonnull()   ) fillPlotsNumerator("newSda"   , denominator, den, TrackInfo(row.newSda   ), npuOOT, npuIT);
  if(row.newSdaUpd.isNonnull()) fillPlotsNumerator("newSdaUpd", denominator, den, TrackInfo(row.newSdaUpd), npuOOT, npuIT);
  if(row.newGl.isNonnull()    ) fillPlotsNumerator("newGl"    , denominator, den, TrackInfo(row.newGl    ), npuOOT, npuIT);
}


int MuonGeneralAnalyzer::fillMatchedPlots(const char* label, const reco::TrackCollection& trackCollection, const reco::GenParticle& gp, int npuOOT, int npuIT){
  float geneta    = gp.eta();
  float genaeta   = fabs(geneta);
  // if(genaeta>2.4) continue;
  float genpt     = gp.pt();
  float genphi    = gp.phi();
  float gencharge = gp.charge();
  
  
  int closestIndex = -1;
  // std::cout << label << ": ";
  auto closestRecTrk = std::min_element(trackCollection.begin(), trackCollection.end(), [geneta, genphi](const reco::Track& a, const reco::Track& b){
      return deltaR(geneta, genphi, a.eta(), a.phi()) < deltaR(geneta, genphi, b.eta(), b.phi());
    });
  if(closestRecTrk == trackCollection.end()){
    // std::cout << "no match" << std::endl;
    return -1;
  }
  
  // commonHitCounter.matchingTracks()
  
  auto track = *closestRecTrk;
  
  float receta = track.eta();
  float recphi = track.phi();
  float recpt = track.pt();
  auto reccharge = track.charge();
  
  float dR = deltaR(geneta, genphi, receta, recphi);
  if(dR <= 0.3){
    closestIndex = std::distance(closestRecTrk, trackCollection.begin());
    // std::cout << closestIndex << std::endl;

    char etaRegion[8];
    if     (genaeta<0.9) strncpy(etaRegion, "bar", sizeof(etaRegion));
    else if(genaeta<1.2) strncpy(etaRegion, "ovl", sizeof(etaRegion));
    else                 strncpy(etaRegion, "end", sizeof(etaRegion));
    
    unsigned short nhits    = track.hitPattern().numberOfMuonHits();
    unsigned short nvalhits = track.hitPattern().numberOfValidMuonHits();
    unsigned short nsegs(0), nvalsegs(0);
    for(auto rh : track.recHits()) {
      DetId rhDetId = rh->geographicalId();
      if(rhDetId.det()==DetId::Detector::Muon && (rhDetId.subdetId()==MuonSubdetId::DT || rhDetId.subdetId()==MuonSubdetId::CSC)) {
	nsegs++;
	if(rh->isValid()) nvalsegs++;
      }
    }
    
    // Numerator
    hists_[Form("%s_eff_num_aeta" , label)]->Fill(genaeta);
    hists_[Form("%s_eff_num_pt"   , label)]->Fill(genpt);
    hists_[Form("%s_eff_num_phi"  , label)]->Fill(genphi);
    hists_[Form("%s_eff_num_pu"   , label)]->Fill(npuIT);
    hists_[Form("%s_eff_num_allpu", label)]->Fill(npuIT+npuOOT);
      
    // Resolution
    hists_[Form("%s_res_pt" , label)]->Fill((reccharge/gencharge)*(genpt/recpt) - 1.);
    hists_[Form("%s_res_eta", label)]->Fill(receta - geneta);
    hists_[Form("%s_res_phi", label)]->Fill(deltaPhi(recphi, genphi));
    hists_[Form("%s_res_dr" , label)]->Fill(dR);
      
    // Track quality
    hists_[Form("%s_nhits"   , label)]->Fill(nhits);
    hists_[Form("%s_nvalhits", label)]->Fill(nvalhits);
    hists_[Form("%s_nsegs"   , label)]->Fill(nsegs);
    hists_[Form("%s_nvalsegs", label)]->Fill(nvalsegs);
    
    hists_[Form("%s_nhits_%s"   , label, etaRegion)]->Fill(nhits);
    hists_[Form("%s_nvalhits_%s", label, etaRegion)]->Fill(nvalhits);
    hists_[Form("%s_nsegs_%s"   , label, etaRegion)]->Fill(nsegs);
    hists_[Form("%s_nvalsegs_%s", label, etaRegion)]->Fill(nvalsegs);
  }
  
  return closestIndex;
}


void MuonGeneralAnalyzer::fillFakePlots(const char* label, const reco::TrackCollection& trackCollection, const std::vector<bool>& mask, int npuOOT, int npuIT){
  for(size_t i = 0; i < mask.size(); ++i){
    // std::cout << "Fake " << label << ": " << i << ' ';
    if(mask.at(i)){
      // std::cout << "continue" << std::endl;
      continue;
    }
    // std::cout<<"fill plots" << std::endl;
    
    const reco::Track& track = trackCollection.at(i);
    float receta = track.eta();
    float recaeta = fabs(receta);
    float recphi = track.phi();
    float recpt = track.pt();
    // auto reccharge = track.charge();
    
    char etaRegion[8];
    if     (recaeta<0.9) strncpy(etaRegion, "bar", sizeof(etaRegion));
    else if(recaeta<1.2) strncpy(etaRegion, "ovl", sizeof(etaRegion));
    else                 strncpy(etaRegion, "end", sizeof(etaRegion));
    
    unsigned short nhits    = track.hitPattern().numberOfMuonHits();
    unsigned short nvalhits = track.hitPattern().numberOfValidMuonHits();
    unsigned short nsegs(0), nvalsegs(0);
    for(auto rh : track.recHits()) {
      DetId rhDetId = rh->geographicalId();
      if(rhDetId.det()==DetId::Detector::Muon && (rhDetId.subdetId()==MuonSubdetId::DT || rhDetId.subdetId()==MuonSubdetId::CSC)) {
	nsegs++;
	if(rh->isValid()) nvalsegs++;
      }
    }
    
    for(const char* suffix : {"perevt", "permuo"}){
    	hists_[Form("%s_fake_pt_%s"   , label, suffix)]->Fill(recpt);
    	hists_[Form("%s_fake_aeta_%s" , label, suffix)]->Fill(fabs(receta));
    	hists_[Form("%s_fake_phi_%s"  , label, suffix)]->Fill(recphi);
    	hists_[Form("%s_fake_pu_%s"   , label, suffix)]->Fill(npuIT);
    	hists_[Form("%s_fake_allpu_%s", label, suffix)]->Fill(npuIT+npuOOT);
	
    	hists_[Form("%s_fake_nhits_%s"   , label, suffix)]->Fill(nhits);
    	hists_[Form("%s_fake_nvalhits_%s", label, suffix)]->Fill(nvalhits);
    	hists_[Form("%s_fake_nsegs_%s"   , label, suffix)]->Fill(nsegs);
    	hists_[Form("%s_fake_nvalsegs_%s", label, suffix)]->Fill(nvalsegs);
	
    	hists_[Form("%s_fake_nhits_%s_%s"   , label, etaRegion, suffix)]->Fill(nhits);
    	hists_[Form("%s_fake_nvalhits_%s_%s", label, etaRegion, suffix)]->Fill(nvalhits);
    	hists_[Form("%s_fake_nsegs_%s_%s"   , label, etaRegion, suffix)]->Fill(nsegs);
    	hists_[Form("%s_fake_nvalsegs_%s_%s", label, etaRegion, suffix)]->Fill(nvalsegs);
    }
  }
}


void MuonGeneralAnalyzer::endPlots(const char* label){
  const char divopt[] = "cl=0.683 b(1,1) mode";
  // try{ graphs_[Form("%s_eff_pt_err"   , label)]->Divide(hists_[Form("%s_eff_num_pt"   , label)], hists_["eff_den_pt"   ], divopt); } catch(cms::Exception& ex) {}
  // try{ graphs_[Form("%s_eff_aeta_err" , label)]->Divide(hists_[Form("%s_eff_num_aeta" , label)], hists_["eff_den_aeta" ], divopt); } catch(cms::Exception& ex) {}
  // try{ graphs_[Form("%s_eff_phi_err"  , label)]->Divide(hists_[Form("%s_eff_num_phi"  , label)], hists_["eff_den_phi"  ], divopt); } catch(cms::Exception& ex) {}
  // try{ graphs_[Form("%s_eff_pu_err"   , label)]->Divide(hists_[Form("%s_eff_num_allpu", label)], hists_["eff_den_allpu"], divopt); } catch(cms::Exception& ex) {}
  // try{ graphs_[Form("%s_eff_allpu_err", label)]->Divide(hists_[Form("%s_eff_num_allpu", label)], hists_["eff_den_allpu"], divopt); } catch(cms::Exception& ex) {}
  
  float n_events = hists_["evt_muo_counter"]->GetBinContent(1);
  float n_muons  = hists_["evt_muo_counter"]->GetBinContent(2);
  
  // hists_[Form("%s_fake_pt_perevt"          , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_aeta_perevt"        , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_phi_perevt"         , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_pu_perevt"          , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_allpu_perevt"       , label)]->Scale( 1./n_events );

  // hists_[Form("%s_fake_nhits_perevt"       , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nhits_bar_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nhits_ovl_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nhits_end_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalhits_perevt"    , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalhits_bar_perevt", label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalhits_ovl_perevt", label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalhits_end_perevt", label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nsegs_perevt"       , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nsegs_bar_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nsegs_ovl_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nsegs_end_perevt"   , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalsegs_perevt"    , label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalsegs_bar_perevt", label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalsegs_ovl_perevt", label)]->Scale( 1./n_events );
  // hists_[Form("%s_fake_nvalsegs_end_perevt", label)]->Scale( 1./n_events );

  // hists_[Form("%s_fake_pt_permuo"          , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_aeta_permuo"        , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_phi_permuo"         , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_pu_permuo"          , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_allpu_permuo"       , label)]->Scale( 1./n_muons );

  // hists_[Form("%s_fake_nhits_permuo"       , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nhits_bar_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nhits_ovl_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nhits_end_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalhits_permuo"    , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalhits_bar_permuo", label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalhits_ovl_permuo", label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalhits_end_permuo", label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nsegs_permuo"       , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nsegs_bar_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nsegs_ovl_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nsegs_end_permuo"   , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalsegs_permuo"    , label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalsegs_bar_permuo", label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalsegs_ovl_permuo", label)]->Scale( 1./n_muons );
  // hists_[Form("%s_fake_nvalsegs_end_permuo", label)]->Scale( 1./n_muons );
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuonGeneralAnalyzer);
