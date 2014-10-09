#ifndef FLASHgg_DiPhotonTagBase_h
#define FLASHgg_DiPhotonTagBase_h

#include "flashgg/MicroAODFormats/interface/DiPhotonCandidate.h"
#include "flashgg/TagFormats/interface/DiPhotonMVAResult.h"

namespace flashgg {

  class DiPhotonTagBase : public reco::CompositeCandidate {
  public:
    DiPhotonTagBase();
    ~DiPhotonTagBase();
    DiPhotonTagBase(edm::Ptr<DiPhotonCandidate>,edm::Ptr<DiPhotonMVAResult>);
    const DiPhotonCandidate* diPhoton() const;
    const DiPhotonMVAResult diPhotonMVA() const;
  private:
    DiPhotonMVAResult mva_result_;
  };

}

#endif
