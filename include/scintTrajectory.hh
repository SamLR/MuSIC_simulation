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
#ifndef scintTrajectory_h
#define scintTrajectory_h 1

#include "G4Trajectory.hh"
#include "G4Allocator.hh"
#include "G4ios.hh" 
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;                   // Forward declaration.

class scintTrajectory : public G4Trajectory
{
public:
  scintTrajectory();
  scintTrajectory(const G4Track* aTrack);
  scintTrajectory(scintTrajectory &);
  virtual ~scintTrajectory();
  
  virtual void DrawTrajectory(G4int i_mode=0) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void SetDrawTrajectory(G4bool b){drawit=b;}
  void WLS(){wls=true;}
  void SetForceDrawTrajectory(G4bool b){forceDraw=b;}
  void SetForceNoDrawTrajectory(G4bool b){forceNoDraw=b;}
private:
  G4bool wls;
  G4bool drawit;
  G4bool forceNoDraw;
  G4bool forceDraw;
  G4ParticleDefinition* particleDefinition;
};

extern G4Allocator<scintTrajectory> scintTrajectoryAllocator;

inline void* scintTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)scintTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void scintTrajectory::operator delete(void* aTrajectory)
{
  scintTrajectoryAllocator.FreeSingle((scintTrajectory*)aTrajectory);
}

#endif
