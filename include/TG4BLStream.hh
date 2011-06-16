//-*-C++-*-
// $Id: TG4BLStream.hh 149 2008-07-02 06:33:39Z comet $

#ifndef TG4BLStream_HH
#define TG4BLStream_HH 1

static const char* TG4BLStream_hh =
"@(#) $Id: TG4BLStream.hh 149 2008-07-02 06:33:39Z comet $";

#include <iostream>      // C++
#include <fstream>       // C++

#include "G4ThreeVector.hh"
#include "globals.hh"


class TG4BLStream
{

private:

  G4int fFileStatus;
  G4String fFileName;

  std::ifstream fFileIn;

private:

  G4bool fIsValid;

  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
  G4double fMomentum;
  G4double fProperTime;
  G4double fWeight;
  G4int fPDGEncoding;
  G4int fEventID;
  G4int fTrackID;
  G4int fParentID;

public:
  TG4BLStream();
  TG4BLStream(G4String val);
  ~TG4BLStream();

public:
  void Open(G4String val);

  void Close();

  void Next();

  G4bool IsValid() {return fIsValid;};

  G4bool IsOpen() {return fFileIn.is_open();};

  const G4ThreeVector& GetPosition() {return fPosition;};
  const G4ThreeVector& GetDirection() {return fDirection;};
  const G4double GetMomentum() {return fMomentum;};
  const G4double GetProperTime() {return fProperTime;};
  const G4double GetWeight() {return fWeight;};
  const G4int GetPDGEncoding() {return fPDGEncoding;};
  const G4int GetEventID() {return fEventID;};
  const G4int GetTrackID() {return fTrackID;};
  const G4int GetParentID() {return fParentID;};

};

#endif
