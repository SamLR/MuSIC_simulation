//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef scintPMTHit_h
#define scintPMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

class G4VTouchable;

class scintPMTHit : public G4VHit
{
public:
  
  scintPMTHit();
  ~scintPMTHit();
  scintPMTHit(const scintPMTHit &right);

  const scintPMTHit& operator=(const scintPMTHit &right);
  G4int operator==(const scintPMTHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

  inline void SetDrawit(G4bool b){drawit=b;}
  inline G4bool GetDrawit(){return drawit;}

  inline void IncPhotonCount(){photons++;}
  inline G4int GetPhotonCount(){return photons;}

  inline void SetPMTNumber(G4int n) { pmtNumber = n; }
  inline G4int GetPMTNumber() { return pmtNumber; }

  inline void SetPMTPhysVol(G4VPhysicalVolume* physVol){this->physVol=physVol;}
  inline G4VPhysicalVolume* GetPMTPhysVol(){return physVol;}

  inline void SetPMTPos(G4double x,G4double y,G4double z){
    pos=G4ThreeVector(x,y,z);
  }
  
private:

  G4double time;

public:

  inline G4ThreeVector GetPMTPos(){return pos;}
  //inline G4double GetTime(){return time;}
  
  inline void SetTime(G4double t)             {time=t;};
  inline G4double GetTime()                   const {return time;};
  
private:
  G4int pmtNumber;
  G4int photons;
  G4ThreeVector pos;
  G4VPhysicalVolume* physVol;
  G4bool drawit;

};

typedef G4THitsCollection<scintPMTHit> scintPMTHitsCollection;

extern G4Allocator<scintPMTHit> scintPMTHitAllocator;

inline void* scintPMTHit::operator new(size_t){
  void *aHit;
  aHit = (void *) scintPMTHitAllocator.MallocSingle();
  return aHit;
}

inline void scintPMTHit::operator delete(void *aHit){
  scintPMTHitAllocator.FreeSingle((scintPMTHit*) aHit);
}

#endif


