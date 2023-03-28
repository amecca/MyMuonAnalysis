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
  typedef std::pair<reco::TrackRef, float> TrackRefProb;
  /* typedef edm::AssociationMap<edm::OneToOne<reco::TrackCollection, reco::TrackCollection>> map_type; */

  explicit CommonHitCounter(const edm::ParameterSet&);

  /* std::vector<TrackRefProb> matchingTracks(const reco::TrackRef& trackRefFrom, reco::TrackCollection& tracksTo, bool flatten) const; */
  std::vector<TrackRefProb> matchingTracks(const reco::Track&       trackFrom, reco::TrackCollection& tracksTo, bool flatten) const;
  /* map_type associateTracks(const reco::TrackCollection&, reco::TrackCollection&, bool flatten=false) const; */
  
  std::tuple<int, int, int> countMatchingHits(const std::vector<TrackingRecHit*>&, const std::vector<TrackingRecHit*>&) const;
  std::tuple<int, int, int> countMatchingHits(const reco::Track& t1, const reco::Track& t2, bool flatten=false) const;

 protected:
  std::vector<TrackingRecHit*> hitsFromTrack(const reco::Track& t, bool flatten=false) const;
  std::vector<TrackingRecHit*> flattenTrackingRecHits(const std::vector<TrackingRecHit*>&) const;
  
  void dumpHit(const TrackingRecHit& hit, std::ostream& ostr) const;
  
 private:
  bool debug_;
};

#endif
