//
//  AliCustomPIDResponse.cxx
//
//  Created by Maximiliano Puccio on 14/11/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AliCustomPIDResponse.h"

#include <limits>
#define NEPS std::numeric_limits<float>::epsilon()

#include "AliExternalTrackParam.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

typedef AliPID::EParticleType Particle;
#define MD 1.875612859f
#define M2D MD*MD

ClassImp(AliCustomPIDResponse)

Float_t AliCustomPIDResponse::NumberOfSigmasTOF(const AliVParticle *track, Particle type) const {
  if (type != AliPID::kDeuteron) {
    return AliPIDResponse::NumberOfSigmasTOF(track, type);
  } else {
    const AliVTrack *vT = static_cast<const AliVTrack*>(track);
    const Double_t p = vT->GetTPCmomentum();
    Float_t time0 = fTOFResponse.GetStartTime(vT->P());
    Float_t time = vT->GetTOFsignal() - time0;
    Float_t beta = vT->GetIntegratedLength() / (2.99792457999999984e-02 * time);
    if (beta < (1.f - NEPS)) {
      Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
      const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
      return dm / 0.15f;
    } else return 999.f;
  }
}

Float_t AliCustomPIDResponse::NumberOfSigmasTPC(const AliVParticle *track, Particle type) const {
  if (type != AliPID::kDeuteron) {
    return AliPIDResponse::NumberOfSigmasTOF(track, type);
  } else {
    const AliVTrack *vT = static_cast<const AliVTrack*>(track);
    const Double_t p = vT->GetTPCmomentum();
    Double_t est = AliExternalTrackParam::BetheBlochAleph(p/1.875612, 4.69637, 7.51827,
                                                          0.0183746,2.60,2.7);
    return (vT->GetTPCsignal() - est) / 0.1f;
  }
}
