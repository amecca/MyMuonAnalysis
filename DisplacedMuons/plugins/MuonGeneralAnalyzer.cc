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
//#include "TCanvas.h"
#include "TString.h"

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
    reco::GenParticleRef genMuon;
    reco::TrackRef outerTrack;
    reco::TrackRef sdaTrack;
    reco::TrackRef sdaUpdTrack;
    
    // reco::TrackRef oldSdaTrack;
    // reco::TrackRef oldSdaUpdTrack;
    // reco::TrackRef newSdaTrack;
    // reco::TrackRef newSdaUpdTrack;
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
  void createPlots();
  int fillMatchedPlots(const char* label, const reco::TrackCollection& trackCollection, const reco::GenParticle& gp, int npuOOT, int npuIT);
  void fillFakePlots(const char* label, const reco::TrackCollection& trackCollection, const std::vector<bool>& mask, int npuOOT, int npuIT);
  void endPlots(const char* label);

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
  edm::EDGetTokenT<reco::TrackCollection> saTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> saUpdTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> globalTracksToken_;
  
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
  processTag_(iConfig.getParameter<edm::InputTag>("saTracksTag").process()),
  puInfoTag_( consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag")) ),
  gpToken_( consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gpCollectionTag")) ),
  muonToken_(        consumes<reco::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonTag"       )) ),
  saTracksToken_(    consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("saTracksTag"   )) ),
  saUpdTracksToken_( consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("saUpdTracksTag")) ),
  globalTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("globalTracks"  )) ),
  
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
  createPlots();
}


MuonGeneralAnalyzer::~MuonGeneralAnalyzer() { 
  if(debugPrint) std::cout << "Inside Destructor" << std::endl;
}


//
// member functions
//

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
    if(debugPrint) std::cout << "Found Muon collection!" << muons->size() << ')' << std::endl;
  }
  
  // GL tracks
  edm::Handle<reco::TrackCollection> glTracks;
  iEvent.getByToken(globalTracksToken_, glTracks);
  if(!glTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "GLOBAL track collection not valid!" << std::endl;
    return;
  }
  else {
    if(debugPrint) std::cout << "Found GLOBAL track collection! (" << glTracks->size() << ')' << std::endl;
  }

  // SA tracks
  edm::Handle<reco::TrackCollection> saTracks;
  iEvent.getByToken(saTracksToken_, saTracks);
  if(!saTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "SA track collection not valid!" << std::endl;
    return;
  }
  else {
    if(debugPrint) std::cout << "Found SA track collection! (" << saTracks->size() << ')' << std::endl;
  }

  // SA-updated at vertex tracks
  edm::Handle<reco::TrackCollection> saUpdTracks;
  iEvent.getByToken(saUpdTracksToken_, saUpdTracks);
  if(!saUpdTracks.isValid()) {
    if(outputPrint || debugPrint) std::cout << "SA Updated track collection not valid!" << std::endl;
    return;
  }
  else {
    if(debugPrint) std::cout << "Found SA Updated track collection! (" << saUpdTracks->size() << ')' << std::endl;
  }
  
  // ################################################################################
  // test
  for(auto sa : *saTracks){
    for(auto saUpd : *saUpdTracks){
      auto tuple_hits = commonHitCounter.countMatchingHits(sa, saUpd, true);
      std::cout << "tuple_hits: " << std::get<0>(tuple_hits) << ' ' << std::get<1>(tuple_hits) << ' ' << std::get<2>(tuple_hits) << '\n';
    }
  }
  return;
  // ################################################################################

  //std::cout << "--- 1" << std::endl;
  reco::TrackCollection sa0hitsTracks;
  for(auto trk : *saTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      sa0hitsTracks.push_back(trk);
  }

  reco::TrackCollection saupd0hitsTracks;
  for(auto trk : *saUpdTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      saupd0hitsTracks.push_back(trk);
  }
  
  reco::TrackCollection gl0hitsTracks;
  for(auto trk : *glTracks){
    if(trk.hitPattern().numberOfValidMuonHits() == 0)
      gl0hitsTracks.push_back(trk);
  }
  
  std::vector<bool>    saMask(   saTracks->size(), false);
  std::vector<bool> saupdMask(saUpdTracks->size(), false);
  std::vector<bool>    glMask(   glTracks->size(), false);
  
  //std::cout << "--- 2" << std::endl;
  
  // ********************************************************************************
  // Make individual associations
  reco::TrackCollection outerTracks;
  outerTracks.reserve(muons->size());
  for(auto muon : *muons){
    reco::TrackRef outerTrack = muon.outerTrack();
    if(outerTrack.isNonnull())  // isAvailable()
      outerTracks.push_back(*outerTrack);
  }

  CommonHitCounter::map_type outerToSdA    = commonHitCounter.matchingTrackCollections(outerTracks, *saTracks   );
  CommonHitCounter::map_type outerToSdAUpd = commonHitCounter.matchingTrackCollections(outerTracks, *saUpdTracks);
  
  // Contruct the table 
  std::vector<AssociationRow> associationTable;
  associationTable.reserve(muons->size());
  for(unsigned int i = 0; i < muons->size(); ++i){
    AssociationRow row;
    row.muon = reco::MuonRef(muons, i);
    row.outerTrack = row.muon->outerTrack();

    if(row.outerTrack.isNonnull()){
      auto it_SdA    = outerToSdA   .find(row.outerTrack);
      auto it_SdAUpd = outerToSdAUpd.find(row.outerTrack);
      row.sdaTrack    = it_SdA    != outerToSdA.end()    ? it_SdA->val    : reco::TrackRef(saTracks.id()   );
      row.sdaUpdTrack = it_SdAUpd != outerToSdAUpd.end() ? it_SdAUpd->val : reco::TrackRef(saUpdTracks.id());
    }

    float muEta = row.muon->eta();
    float muPhi = row.muon->phi();
    auto it_closestGen = std::min_element(genMuons.begin(), genMuons.end(), [muEta, muPhi](const reco::GenParticle& a, const reco::GenParticle& b) {
	return deltaR(muEta, muPhi, a.eta(), a.phi()) < deltaR(muEta, muPhi, b.eta(), b.phi());
      } );
    row.genMuon = it_closestGen != genMuons.end() && deltaR(muEta, muPhi, it_closestGen->eta(), it_closestGen->phi()) < 0.2 ? 
      reco::GenParticleRef(&genMuons, std::distance(genMuons.begin(), it_closestGen)) : 
      reco::GenParticleRef();

    associationTable.push_back(std::move(row));
  }
  
  
  std::cout << Form("%6s | %6s | %6s | %6s | %6s\n", "muon", "genMu", "outer", "SdA", "SdAUpd");
  for(const AssociationRow& row : associationTable)
    std::cout << Form("%6d | %6d | %6s | %6d | %6d\n",
		      row.muon.key(),
		      row.genMuon.isNonnull() ? row.genMuon.key() : -1,
		      row.outerTrack.isNonnull() ? "OK" : "-",
		      row.sdaTrack   .isNonnull() ? row.sdaTrack   .key() : -1,
		      row.sdaUpdTrack.isNonnull() ? row.sdaUpdTrack.key() : -1
		      );
  
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
    int saIdx = fillMatchedPlots("sa"        , *saTracks       , gp, npuOOT, npuIT);
    if(saIdx >= 0)
      saMask.at(saIdx) = true;
    // std::cout << "sa: "<<saIdx<<std::endl;
    // fillMatchedPlots("sa0hits"   , sa0hitsTracks  , gp, npuOOT, npuIT);
    
    // SA-updated tracks
    int saupdIdx = fillMatchedPlots("saupd"     , *saUpdTracks    , gp, npuOOT, npuIT);
    if(saupdIdx >= 0)
      saupdMask.at(saupdIdx) = true;
    // std::cout << "saupd: "<<saupdIdx<<std::endl;
    // fillMatchedPlots("saupd0hits", saupd0hitsTracks, gp, npuOOT, npuIT);

    // GLOBAL tracks
    int glIdx = fillMatchedPlots("gl"        , *glTracks       , gp, npuOOT, npuIT);
    if(glIdx >= 0)
      glMask.at(glIdx) = true;
    // std::cout << "gl: "<<glIdx<<std::endl;
    // fillMatchedPlots("gl0hits"   , gl0hitsTracks   , gp, npuOOT, npuIT);
  }
  
  fillFakePlots("sa"   , *saTracks   , saMask   , npuOOT, npuIT);
  fillFakePlots("saupd", *saUpdTracks, saupdMask, npuOOT, npuIT);
  fillFakePlots("gl"   , *glTracks   , glMask   , npuOOT, npuIT);

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


void MuonGeneralAnalyzer::createPlots(){
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
  try{ graphs_[Form("%s_eff_pt_err"   , label)]->Divide(hists_[Form("%s_eff_num_pt"   , label)], hists_["eff_den_pt"   ], divopt); } catch(cms::Exception& ex) {}
  try{ graphs_[Form("%s_eff_aeta_err" , label)]->Divide(hists_[Form("%s_eff_num_aeta" , label)], hists_["eff_den_aeta" ], divopt); } catch(cms::Exception& ex) {}
  try{ graphs_[Form("%s_eff_phi_err"  , label)]->Divide(hists_[Form("%s_eff_num_phi"  , label)], hists_["eff_den_phi"  ], divopt); } catch(cms::Exception& ex) {}
  try{ graphs_[Form("%s_eff_pu_err"   , label)]->Divide(hists_[Form("%s_eff_num_allpu", label)], hists_["eff_den_allpu"], divopt); } catch(cms::Exception& ex) {}
  try{ graphs_[Form("%s_eff_allpu_err", label)]->Divide(hists_[Form("%s_eff_num_allpu", label)], hists_["eff_den_allpu"], divopt); } catch(cms::Exception& ex) {}
  
  float n_events = hists_["evt_muo_counter"]->GetBinContent(1);
  float n_muons  = hists_["evt_muo_counter"]->GetBinContent(2);
  
  hists_[Form("%s_fake_pt_perevt"          , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_aeta_perevt"        , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_phi_perevt"         , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_pu_perevt"          , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_allpu_perevt"       , label)]->Scale( 1./n_events );
  
  hists_[Form("%s_fake_nhits_perevt"       , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nhits_bar_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nhits_ovl_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nhits_end_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalhits_perevt"    , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalhits_bar_perevt", label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalhits_ovl_perevt", label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalhits_end_perevt", label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nsegs_perevt"       , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nsegs_bar_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nsegs_ovl_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nsegs_end_perevt"   , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalsegs_perevt"    , label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalsegs_bar_perevt", label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalsegs_ovl_perevt", label)]->Scale( 1./n_events );
  hists_[Form("%s_fake_nvalsegs_end_perevt", label)]->Scale( 1./n_events );
  
  hists_[Form("%s_fake_pt_permuo"          , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_aeta_permuo"        , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_phi_permuo"         , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_pu_permuo"          , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_allpu_permuo"       , label)]->Scale( 1./n_muons );
  
  hists_[Form("%s_fake_nhits_permuo"       , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nhits_bar_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nhits_ovl_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nhits_end_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalhits_permuo"    , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalhits_bar_permuo", label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalhits_ovl_permuo", label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalhits_end_permuo", label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nsegs_permuo"       , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nsegs_bar_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nsegs_ovl_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nsegs_end_permuo"   , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalsegs_permuo"    , label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalsegs_bar_permuo", label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalsegs_ovl_permuo", label)]->Scale( 1./n_muons );
  hists_[Form("%s_fake_nvalsegs_end_permuo", label)]->Scale( 1./n_muons );
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuonGeneralAnalyzer);
