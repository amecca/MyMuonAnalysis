#ifndef COMMONHITCOUNTER_H
#define COMMONHITCOUNTER_H

// -*- C++ -*-
//
// Package:    RecoMuon/TracksAssociator
// Class:      TracksAssociator
//
/**\class TracksAssociator TracksAssociator.cc RecoMuon/TracksAssociator/src/CommonHitCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Mecca
//         Created:  Fri, 20 Jan 2023 11:06:57 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/* #include "DataFormats/Common/interface/OneToManyWithQuality.h" */
#include "DataFormats/Common/interface/OneToOne.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

//
// class declaration
//

class CommonHitCounter {
 public:
  template <class T>
  struct HashRef {
    // Hash for edm::Ref
    // bytes  :   7 | 6 | 5 | 4 | 3 | 2 | 1 | 0
    // content:    proc |  prod |      key

    static constexpr unsigned int bitshift_half    = (sizeof(size_t)/2) * 8;
    static constexpr unsigned int bitshift_quarter = (sizeof(size_t)/4) * 8;
    size_t operator()(const edm::Ref<T>& ref) const{
      edm::ProductID id = ref.id();
      unsigned short process = id.processIndex();  // bytes 7-6
      unsigned short product = id.productIndex();  // bytes 5-4
      unsigned int upper_half = (process << bitshift_quarter) + product; // bytes 7-4
      // ref.key() is added without shift, and allowed to overflow into the upper half
      return (size_t(upper_half) << bitshift_half) + ref.key();
    }
  };

  typedef std::pair<reco::TrackRef, std::vector<TrackingRecHit*>> TrackrefHitsPair;
  typedef std::pair<reco::TrackRef, float> TrackRefProb;
  /* typedef edm::AssociationMap<edm::OneToOne<reco::TrackCollection, reco::TrackCollection>> map_type; */
  typedef std::unordered_map<reco::TrackRef, reco::TrackRef, HashRef<reco::TrackCollection>> map_type;

  explicit CommonHitCounter(const edm::ParameterSet&);

  reco::TrackRef matchHits (const std::vector<TrackingRecHit*>& hitsFrom, const reco::TrackCollection& tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;
  reco::TrackRef matchTrack(const reco::Track& trackFrom                , const reco::TrackCollection& tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;

  map_type matchingTrackCollections(const reco::TrackCollection&       tracksFrom, const reco::TrackCollection&       tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;
  map_type matchingTrackCollections(const reco::TrackCollection&       tracksFrom, const std::vector<reco::TrackRef>& tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;
  map_type matchingTrackCollections(const std::vector<reco::TrackRef>& tracksFrom, const reco::TrackCollection&       tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;
  map_type matchingTrackCollections(const std::vector<reco::TrackRef>& tracksFrom, const std::vector<reco::TrackRef>& tracksTo, bool flatten=true, const std::unordered_set<DetId::Detector>& d={}) const;

  std::tuple<int, int, int> countMatchingHits(const std::vector<TrackingRecHit*>&, const std::vector<TrackingRecHit*>&) const;
  std::tuple<int, int, int> countMatchingHits(const reco::Track&    t1, const reco::Track& t2   , bool flatten, const std::unordered_set<DetId::Detector>& d) const;
  std::tuple<int, int, int> countMatchingHits(const reco::TrackRef& t1, const reco::TrackRef& t2, bool flatten, const std::unordered_set<DetId::Detector>& d) const{
    return countMatchingHits(t1, t2, flatten, d);
  }

  std::vector<TrackingRecHit*> flattenTrackingRecHits(const std::vector<TrackingRecHit*>&) const;

 protected:
  map_type doMatchTrackCollections(std::vector<TrackrefHitsPair>& trackToHitsFrom, std::vector<TrackrefHitsPair>& trackToHitsTo) const;

  std::vector<TrackrefHitsPair> prepareTrackRefToHits(const std::vector<reco::TrackRef>& trackRefs, bool flatten, const std::unordered_set<DetId::Detector>&) const;
  std::vector<TrackrefHitsPair> prepareTrackRefToHits(const reco::TrackCollection& tracks         , bool flatten, const std::unordered_set<DetId::Detector>&) const;

  std::vector<TrackingRecHit*> hitsFromTrack(const reco::Track& t, bool flatten, const std::unordered_set<DetId::Detector>& d) const;
  
  void dumpHit(const TrackingRecHit& hit, std::ostream& ostr) const;
  
 private:
  bool debug_;
  float matchingFractionCut_;  // Fraction of common hits required to declare that two tracks are equal
};

#endif
