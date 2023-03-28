#include "MyMuonAnalysis/DisplacedMuons/interface/CommonHitCounter.h"

CommonHitCounter::CommonHitCounter(const edm::ParameterSet& iConfig) :
  debug_( iConfig.exists("debug") && iConfig.getUntrackedParameter<bool>("debug") )
{}


std::vector<CommonHitCounter::TrackRefProb> CommonHitCounter::matchingTracks(const reco::Track& trackFrom, reco::TrackCollection& tracksTo, bool flatten) const{
  std::vector<TrackRefProb> result;
  
  std::vector<TrackingRecHit*> hitsFrom = hitsFromTrack(trackFrom, flatten);
  for(size_t i = 0; i < tracksTo.size(); ++i){
    const reco::TrackRef trackRefTo(&tracksTo, i);
    std::vector<TrackingRecHit*> hitsTo = hitsFromTrack(*trackRefTo, flatten);
    
    int nMatch, nFrom, nTo;
    std::tie (nMatch, nFrom, nTo) = countMatchingHits(hitsFrom, hitsTo);
    
    if(nMatch){
      // ########################################
      float prob = 2.*nMatch/(nTo + nFrom);
      // ########################################
      result.push_back(make_pair(trackRefTo, prob));
    }
  }
  
  std::sort(result.begin(), result.end(),
	    [](const TrackRefProb& a, const TrackRefProb& b){ return a.second > b.second; }
	    );
  return result;
}


// std::vector<CommonHitCounter::TrackRefProb> CommonHitCounter::matchingTracks(const reco::TrackRef& trackRefFrom, reco::TrackCollection& tracksTo, bool flatten) const{
//   return matchingTracks(*trackRefFrom, tracksTo, flatten);
// }


// CommonHitCounter::map_type CommonHitCounter::associateTracks(const reco::TrackCollection& tracksFrom, reco::TrackCollection& tracksTo, bool flatten) const{
//   typedef std::pair<reco::TrackRef, std::vector<TrackingRecHit*>> TrackHitsPair;
//   typedef std::pair<std::vector<TrackHitsPair>::iterator, float> IterProbPair;
  
//   // Tracks and their hits
//   std::vector<TrackHitsPair> trackToHitsFrom;
//   std::vector<TrackHitsPair> trackToHitsTo  ;
  
//   for(size_t i = 0; i < tracksFrom.size(); ++i)
//     trackToHitsFrom.push_back(std::make_pair(
// 					     reco::TrackRef(&tracksFrom, i),
// 					     hitsFromTrack(tracksFrom.at(i), flatten))
// 			      );
//   for(size_t i = 0; i < tracksTo.size(); ++i)
//     trackToHitsTo.push_back  (std::make_pair(
// 					     reco::TrackRef(&tracksTo  , i),
// 					     hitsFromTrack(tracksTo.at(i)  , flatten))
// 			      );
  
//   map_type associationMap;
//   // std::vector<std::tuple<size_t, size_t, float>> associationMap;  // TEMP
  
//   for(auto& [trackFrom, hitsFrom] : trackToHitsFrom){
//     std::vector<IterProbPair> vIterProbsTo;  // <<Track,Hits>, prob>, prob is the probability that a given Track (to) shares the hits from this Track (from)
    
//     for(auto iterIdxHitsTo = trackToHitsTo.begin(); iterIdxHitsTo != trackToHitsTo.end(); ++iterIdxHitsTo){ //auto& [idxTo, hitsTo]:trackToHitsTo
//       const reco::TrackRef& trackTo = iterIdxHitsTo->first;
//       std::vector<TrackingRecHit*>& hitsTo = iterIdxHitsTo->second;
      
//       int nMatch, nFrom, nTo;
//       std::tie (nMatch, nFrom, nTo) = countMatchingHits(hitsFrom, hitsTo);
      
//       if(nMatch){
// 	float prob = nMatch/(std::min(nTo, nFrom));
// 	if(debug_)
// 	  std::cout << Form("\tTracks: %u, %u \tmatch: %d, 1: %d, 2: %d --> %.1f%%\n",
// 			    trackFrom.key(), trackTo.key(), nMatch, nFrom, nTo, 100*prob
// 			    );
// 	vIterProbsTo.push_back(make_pair(iterIdxHitsTo, prob));
//       }
//     }
    
//     if(vIterProbsTo.size() > 0){
//       std::sort(vIterProbsTo.begin(), vIterProbsTo.end(),
// 		[](const IterProbPair& a, const IterProbPair& b){ return a.second > b.second; }
// 		);
//       auto winner = vIterProbsTo.front();
//       if(debug_)
// 	std::cout << "Choosen: "<<vIterProbsTo.front().second<<" last: "<<vIterProbsTo.back().second<<'\n';
//       if(winner.second > 0.5){
// 	const reco::TrackRef& trackTo = winner.first->first;
// 	// associationMap.push_back( std::make_tuple(idxFrom, idxTo, winner.second) );
// 	associationMap.insert(trackFrom, trackTo);
// 	trackToHitsTo.erase(winner.first);
//       }
//     }
//   }
//   return associationMap;
// }


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


std::vector<TrackingRecHit*> CommonHitCounter::hitsFromTrack(const reco::Track& t, bool flatten) const{
  // This copy is done also to not modify the track's oiginal hits
  auto begin = t.recHitsBegin();
  auto end   = t.recHitsEnd();
  size_t n = std::distance(begin, end);
  std::vector<TrackingRecHit*> hits(n);
  std::copy(begin, end, hits.begin());
  
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
      
      if(h1->rawId() == h2->rawId() /*&& TODO compare local position */ )
  	++matches;
    }
  }
  return std::make_tuple(matches, hits_1.size(), hits_2.size());
}

std::tuple<int, int, int> CommonHitCounter::countMatchingHits(const reco::Track& t1, const reco::Track& t2, bool flatten) const{
  // Negative numbers correspond to errors --> move to matchProbability()
  // if(! (t1->extra().isAvailable() && t2->extra().isAvailable()) )
  //   return std::make_tuple(-1, 0, 0);
  
  std::vector<TrackingRecHit*> hits_1 = hitsFromTrack(t1, flatten);
  std::vector<TrackingRecHit*> hits_2 = hitsFromTrack(t2, flatten);
  
  std::cout << ">>>>>>>>>> Dumping hits_1\n";
  for(auto hit : hits_1){
    dumpHit(*hit, std::cout);
  }
  std::cout << ">>>>>>>>>> Dumping hits_2\n";
  for(auto hit : hits_2){
    dumpHit(*hit, std::cout);
  }
  
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
  
  ostr<<"DetId="<<detName<<':'<<id.subdetId()<<"  type="<<type<<"  #recHits="<<hit.recHits().size()<<"  isValid="<<hit.isValid()<<"  rawId="<<hit.rawId()<<'\n';
}
