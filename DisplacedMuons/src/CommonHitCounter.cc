#include "MyMuonAnalysis/DisplacedMuons/interface/CommonHitCounter.h"

CommonHitCounter::CommonHitCounter(const edm::ParameterSet& iConfig) :
  debug_( iConfig.exists("debug") && iConfig.getUntrackedParameter<bool>("debug") )
  , matchingFractionCut_(iConfig.getParameter<double>("matchingFractionCut"))
{}


reco::TrackRef CommonHitCounter::matchHits(const std::vector<TrackingRecHit*>& hitsFrom, const reco::TrackCollection& tracksTo, bool flatten, const std::unordered_set<DetId::Detector>& detectors) const{
  std::vector<TrackRefProb> vProbTo;

  for(size_t i = 0; i < tracksTo.size(); ++i){
    const reco::TrackRef trackRefTo(&tracksTo, i);
    std::vector<TrackingRecHit*> hitsTo = hitsFromTrack(*trackRefTo, flatten, detectors);

    int nMatch, nFrom, nTo;
    std::tie (nMatch, nFrom, nTo) = countMatchingHits(hitsFrom, hitsTo);

    if(nMatch){
      // ########################################
      float prob = 2.*nMatch/(nTo + nFrom);
      // ########################################
      vProbTo.push_back(make_pair(trackRefTo, prob));
    }
  }

  if(vProbTo.size() > 0){
    auto best = std::max_element(vProbTo.begin(), vProbTo.end(),
				 [](const TrackRefProb& a, const TrackRefProb& b){ return a.second > b.second; }
				 );
    if(best->second > matchingFractionCut_)
      return best->first;
  }

  return reco::TrackRef();
}


reco::TrackRef CommonHitCounter::matchTrack(const reco::Track& trackFrom, const reco::TrackCollection& tracksTo, bool flatten, const std::unordered_set<DetId::Detector>& detectors) const{
  std::vector<TrackingRecHit*> hitsFrom = hitsFromTrack(trackFrom, flatten, detectors);
  return matchHits(hitsFrom, tracksTo, flatten, detectors);
}


std::vector<CommonHitCounter::TrackrefHitsPair> CommonHitCounter::prepareTrackRefToHits(const std::vector<reco::TrackRef>& trackRefs, bool flatten, const std::unordered_set<DetId::Detector>& detectors) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> out;
  out.reserve(trackRefs.size());

  for(const reco::TrackRef& ref : trackRefs)
    out.push_back( std::make_pair(ref, hitsFromTrack(*ref, flatten, detectors)) );

  return out;
}


std::vector<CommonHitCounter::TrackrefHitsPair> CommonHitCounter::prepareTrackRefToHits(const reco::TrackCollection& tracks         , bool flatten, const std::unordered_set<DetId::Detector>& detectors) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> out;
  out.reserve(tracks.size());

  for(size_t i = 0; i < tracks.size(); ++i)
    out.push_back( std::make_pair(reco::TrackRef(&tracks, i), hitsFromTrack(tracks.at(i), flatten, detectors)) );

  return out;
}


CommonHitCounter::map_type CommonHitCounter::doMatchTrackCollections(std::vector<CommonHitCounter::TrackrefHitsPair>& trackToHitsFrom, std::vector<CommonHitCounter::TrackrefHitsPair>& trackToHitsTo) const{
  typedef std::pair<std::vector<TrackrefHitsPair>::iterator, float> IterProbPair;
  
  map_type associationMap;

  for(auto& [trackFrom, hitsFrom] : trackToHitsFrom){
    std::vector<IterProbPair> vIterProbsTo;  // <<Track,Hits>, prob>, prob is the probability that a given Track (to) shares the hits from this Track (from)
    
    for(auto iterIdxHitsTo = trackToHitsTo.begin(); iterIdxHitsTo != trackToHitsTo.end(); ++iterIdxHitsTo){ //auto& [idxTo, hitsTo]:trackToHitsTo
      const reco::TrackRef& trackTo = iterIdxHitsTo->first;
      const std::vector<TrackingRecHit*>& hitsTo = iterIdxHitsTo->second;
      
      int nMatch, nFrom, nTo;
      std::tie (nMatch, nFrom, nTo) = countMatchingHits(hitsFrom, hitsTo);
      
      if(nMatch){
	float prob = nMatch/(std::min(nTo, nFrom));
	if(debug_)
	  std::cout << Form("\tTracks: %u, %u \tmatch: %d, 1: %d, 2: %d --> %.1f%%\n",
			    trackFrom.key(), trackTo.key(), nMatch, nFrom, nTo, 100*prob
			    );
	vIterProbsTo.push_back(make_pair(iterIdxHitsTo, prob));
      }
    }
    
    if(vIterProbsTo.size() > 0){
      std::sort(vIterProbsTo.begin(), vIterProbsTo.end(),
		[](const IterProbPair& a, const IterProbPair& b){ return a.second > b.second; }
		);
      auto winner = vIterProbsTo.front();
      if(debug_)
	std::cout << "Choosen: "<<vIterProbsTo.front().second<<" last: "<<vIterProbsTo.back().second<<'\n';
      if(winner.second > matchingFractionCut_){
	const reco::TrackRef& trackTo = winner.first->first;
	associationMap.insert(std::make_pair(trackFrom, trackTo));
	trackToHitsTo.erase(winner.first);
      }
    }
  }

  return associationMap;
}


CommonHitCounter::map_type CommonHitCounter::matchingTrackCollections(const reco::TrackCollection&      tracksFrom, const reco::TrackCollection& tracksTo      , bool flatten, const std::unordered_set<DetId::Detector>& d) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsFrom = prepareTrackRefToHits(tracksFrom, flatten, d);
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsTo   = prepareTrackRefToHits(tracksTo  , flatten, d);

  return doMatchTrackCollections(trackToHitsFrom, trackToHitsTo);
}

CommonHitCounter::map_type CommonHitCounter::matchingTrackCollections(const reco::TrackCollection&      tracksFrom, const std::vector<reco::TrackRef>& tracksTo, bool flatten, const std::unordered_set<DetId::Detector>& d) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsFrom = prepareTrackRefToHits(tracksFrom, flatten, d);
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsTo   = prepareTrackRefToHits(tracksTo  , flatten, d);

  return doMatchTrackCollections(trackToHitsFrom, trackToHitsTo);
}

CommonHitCounter::map_type CommonHitCounter::matchingTrackCollections(const std::vector<reco::TrackRef>& tracksFrom, const reco::TrackCollection& tracksTo     , bool flatten, const std::unordered_set<DetId::Detector>& d) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsFrom = prepareTrackRefToHits(tracksFrom, flatten, d);
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsTo   = prepareTrackRefToHits(tracksTo  , flatten, d);

  return doMatchTrackCollections(trackToHitsFrom, trackToHitsTo);
}

CommonHitCounter::map_type CommonHitCounter::matchingTrackCollections(const std::vector<reco::TrackRef>& tracksFrom, const std::vector<reco::TrackRef>& tracksTo, bool flatten, const std::unordered_set<DetId::Detector>& d) const{
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsFrom = prepareTrackRefToHits(tracksFrom, flatten, d);
  std::vector<CommonHitCounter::TrackrefHitsPair> trackToHitsTo   = prepareTrackRefToHits(tracksTo  , flatten, d);

  return doMatchTrackCollections(trackToHitsFrom, trackToHitsTo);
}


// ------------ Utility functions used internally ------------
std::vector<TrackingRecHit*> CommonHitCounter::flattenTrackingRecHits(const std::vector<TrackingRecHit*>& hits) const{
  std::vector<TrackingRecHit*> out;
  out.reserve(hits.size()*2);
  for(auto it = hits.begin(); it != hits.end(); ++it){
    if(dynamic_cast<InvalidTrackingRecHit*>(*it) != nullptr)  // check correspondence with isValid()
      continue;  // Skip broken hits
    
    std::vector<TrackingRecHit*> nested  = (*it)->recHits();
    if(nested.size() == 0){
      out.push_back(*it);
    }
    else{
      std::vector<TrackingRecHit*> nested_flat = flattenTrackingRecHits(nested);
      out.insert(out.end(), nested_flat.begin(), nested_flat.end());
    }
  }
  return out;
}


std::vector<TrackingRecHit*> CommonHitCounter::hitsFromTrack(const reco::Track& t, bool flatten, const std::unordered_set<DetId::Detector>& detectors) const{
  // This copy is done also to not modify the track's oiginal hits
  auto begin = t.recHitsBegin();
  auto end   = t.recHitsEnd();
  std::vector<TrackingRecHit*> hits;
  hits.reserve(t.recHitsSize()); //std::distance(begin, end)

  std::copy_if(begin, end, std::back_inserter(hits), [detectors](auto it_hit){
      return dynamic_cast<InvalidTrackingRecHit*>(it_hit) == nullptr                  // skip invalid hits
	&& ( detectors.empty() || detectors.count(it_hit->geographicalId().det()) );  // If given a set of detectors, use oly hits in those detectors
    });

  if(flatten){
    std::vector<TrackingRecHit*> out = flattenTrackingRecHits(hits);
    return out;
  }
  else
    return hits;
}


std::tuple<int, int, int> CommonHitCounter::countMatchingHits(const std::vector<TrackingRecHit*>& hits_1, const std::vector<TrackingRecHit*>& hits_2) const{
  int matches = 0;
  for(const TrackingRecHit* h1 : hits_1){
    // if(! recHitCondition(h1)) continue;
    
    for(const TrackingRecHit* h2 : hits_2){
      // if(! recHitCondition(h2)) continue;
      
      if(h1->rawId() == h2->rawId() && h1->localPosition() == h2->localPosition() )
  	++matches;
    }
  }
  return std::make_tuple(matches, hits_1.size(), hits_2.size());
}


std::tuple<int, int, int> CommonHitCounter::countMatchingHits(const reco::Track& t1, const reco::Track& t2, bool flatten, const std::unordered_set<DetId::Detector>& d) const{
  std::vector<TrackingRecHit*> hits_1 = hitsFromTrack(t1, flatten, d);
  std::vector<TrackingRecHit*> hits_2 = hitsFromTrack(t2, flatten, d);
  return countMatchingHits(hits_1, hits_2);
}


void CommonHitCounter::dumpHit(const TrackingRecHit& hit, std::ostream& ostr) const{
  DetId id = hit.geographicalId();
  char detName[16] = "";
  switch (id.det()){
  case DetId::Detector::Tracker:
    sprintf(detName, "Tracker"); break;
  case DetId::Detector::Muon:
    sprintf(detName, "Muon"); break;
  case DetId::Detector::Ecal:
    sprintf(detName, "Ecal"); break;
  case DetId::Detector::Hcal:
    sprintf(detName, "Hcal"); break;
  case DetId::Detector::Calo:
    sprintf(detName, "Calo"); break;
  case DetId::Detector::Forward:
    sprintf(detName, "Forward"); break;
  case DetId::Detector::VeryForward:
    sprintf(detName, "VeryForward"); break;
  case DetId::Detector::HGCalEE:
    sprintf(detName, "HGCalEE"); break;
  case DetId::Detector::HGCalHSi:
    sprintf(detName, "HGCalHSi"); break;
  case DetId::Detector::HGCalHSc:
    sprintf(detName, "HGCalHSc"); break;
  case DetId::Detector::HGCalTrigger:
    sprintf(detName, "HGCalTrigger"); break;
  default:
    snprintf(detName, 16, "unknown_%d", id.det());
  }
  
  char type[16] = "";
  switch (hit.type()){
  case TrackingRecHit::Type::valid:
    sprintf(type, "valid"); break;
  case TrackingRecHit::Type::missing:
    sprintf(type, "missing"); break;
  case TrackingRecHit::Type::inactive:
    sprintf(type, "inactive"); break;
  case TrackingRecHit::Type::bad:
    sprintf(type, "bad"); break;
  case TrackingRecHit::Type::missing_inner:
    sprintf(type, "missing_inner"); break;
  case TrackingRecHit::Type::missing_outer:
    sprintf(type, "missing_outer"); break;
  case TrackingRecHit::Type::inactive_inner:
    sprintf(type, "inactive_inner"); break;
  case TrackingRecHit::Type::inactive_outer:
    sprintf(type, "inactive_outer"); break;
  default:
    snprintf(type, 16, "unknown_%d", hit.type());
  }

  // TrackingRecHit::Type type = ;
  
  ostr<<"DetId="<<detName<<':'<<id.subdetId()<<"  type="<<type<<"  #recHits="<<hit.recHits().size()<<"  isValid="<<hit.isValid()<<"  rawId="<<hit.rawId()<<"  localPosition="<<hit.localPosition()<<'\n';
}
