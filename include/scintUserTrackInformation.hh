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
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

#ifndef scintUserTrackInformation_h
#define scintUserTrackInformation_h 1

enum scintTrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,
                      hitSphere=16, inactive=14};

/*scintTrackStatus:
  active: still being tracked
  hitPMT: stopped by being detected in a PMT
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  hitSphere: track hit the sphere at some point
  inactive: track is stopped for some reason
   -This is the sum of all stopped flags so can be used to remove stopped flags
  
 */

class scintUserTrackInformation : public G4VUserTrackInformation
{
public:
  scintUserTrackInformation();
  ~scintUserTrackInformation();
  
  //Sets the track status to s (does not check validity of flags)
  void SetTrackStatusFlags(int s){status=s;}
  //Does a smart add of track status flags (disabling old flags that conflict)
  //If s conflicts with itself it will not be detected
  void AddTrackStatusFlag(int s);
  
  int GetTrackStatus()const {return status;}
  
  void IncReflections(){reflections++;}
  G4int GetReflectionCount()const {return reflections;}

  void SetForceDrawTrajectory(G4bool b){forcedraw=b;}
  G4bool GetForceDrawTrajectory(){return forcedraw;}

  inline void Print()const{};
private:
  int status;
  G4int reflections;
  G4bool forcedraw;
};

#endif
