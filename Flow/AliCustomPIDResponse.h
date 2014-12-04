//
//  AliCustomPIDResponse.h
//
//  Created by Maximiliano Puccio on 14/11/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef AliCustomPIDResponse_h
#define AliCustomPIDResponse_h

#include "AliESDpid.h"

class AliVParticle;

class AliCustomPIDResponse : public AliESDpid {
public:
  AliCustomPIDResponse(const AliESDpid& pid) : AliESDpid(pid) {}
  virtual Float_t NumberOfSigmasTOF  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type) const;
  ClassDef(AliCustomPIDResponse, 1)
};

#endif
