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

#ifndef scintScintHit_h
#define scintScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

class scintScintHit : public G4VHit
{
public:
  
  scintScintHit();
  scintScintHit(G4VPhysicalVolume* pVol);
  ~scintScintHit();
  scintScintHit(const scintScintHit &right);
  const scintScintHit& operator=(const scintScintHit &right);
  G4int operator==(const scintScintHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

  inline void SetEdep(G4double de) { edep = de; }  
  inline void AddEdep(G4double de) { edep += de; }
  inline G4double GetEdep() { return edep; }
  
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
  inline G4ThreeVector GetPos() { return pos; }

  inline const G4VPhysicalVolume * GetPhysV() { return physVol; }

  //inline G4double GetTime(){return time;}
  
  void SetTime           (G4double t2)       { time = t2; };
  G4double GetTime()                         { return time; };

  void SetParticle (G4String name) {particleName = name; };
  G4String GetParticle() {return particleName; };

  void SetParticleEnergy (G4double e1)       { particleEnergy = e1; };
  G4double GetParticleEnergy()               { return particleEnergy;};      
  
private:
  G4double edep;
  G4ThreeVector pos;
  const G4VPhysicalVolume* physVol;
  G4double time;
  G4String particleName;
  G4double particleEnergy;

};

typedef G4THitsCollection<scintScintHit> scintScintHitsCollection;

extern G4Allocator<scintScintHit> scintScintHitAllocator;

inline void* scintScintHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) scintScintHitAllocator.MallocSingle();
  return aHit;
}

inline void scintScintHit::operator delete(void *aHit)
{
  scintScintHitAllocator.FreeSingle((scintScintHit*) aHit);
}

#endif


