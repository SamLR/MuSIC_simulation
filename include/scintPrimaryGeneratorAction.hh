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

#ifndef scintPrimaryGeneratorAction_h
#define scintPrimaryGeneratorAction_h 1

static const char* scintPrimaryGeneratorAction_hh =
"@(#) $Id: scintPrimaryGeneratorAction.hh 228 2008-12-24 08:20:56Z comet $";

#include <iostream>      // C++
#include <fstream>       // C++

#include "G4Navigator.hh"
#include "G4TouchableHistory.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "TG4BLStream.hh"

class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;
class scintPrimaryGeneratorMessenger;

class scintPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  scintPrimaryGeneratorAction();
  ~scintPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  void FileOpen(G4String val);
  void FileClose();

  void SetPositionMode(G4String val) { fPositionMode = val;};
  void SetEnergyMode(G4String val);

  void SetEnergySpread(G4double val) {fEnergySpread = val;};
  G4double GetEnergySpread() {return fEnergySpread;};

  void SetPositionSpread(const G4ThreeVector& val) {fPositionSpread = val;};
  const G4ThreeVector& GetPositionSpread() {return fPositionSpread;};

  void SetDirectionSpread(const G4double val) {fDirectionSpread = val;};
  const G4double GetDirectionSpread() {return fDirectionSpread;};

  void SetBeamDispersion(const G4double val) {
    fBeamDispersion = val;
  };
  const G4double GetBeamDispersion() {
    return fBeamDispersion;
  };
  
  void SetCutOffEnergy(const G4double val) {fCutOffEnergy = val;};
  const G4double GetCutOffEnergy() {return fCutOffEnergy;};
  
private:
  G4ParticleGun* fParticleGun;
  scintPrimaryGeneratorMessenger* fGunMessenger;

  G4int fFileStatus;
  TG4BLStream fG4BLStream;

  G4String fPositionMode;
  G4String fEnergyMode;

  G4double fEnergySpread;
  G4ThreeVector fPositionSpread;
  G4double fDirectionSpread;

  G4double fBeamDispersion;
  G4RandGauss* fGauss;
  G4double fCutOffEnergy;


};

#endif


