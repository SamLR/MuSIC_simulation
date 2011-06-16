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
#include "scintPhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4IonPhysics.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "scintGeneralPhysics.hh"
#include "scintEMPhysics.hh"
#include "scintMuonPhysics.hh"
#include "scintHadronPhysics.hh"
#include "scintOpticalPhysics.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"

//#include "G4EmStandardPhysics.hh"

scintPhysicsList::scintPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  // SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new scintGeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new scintEMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new scintMuonPhysics("muon"));

   // Optical Physics
  RegisterPhysics(  new scintOpticalPhysics("optical"));

  //Pion Physics
  //RegisterPhysics(  new HadronPhysicsQGSP_BERT_HP("hadron"));
  RegisterPhysics(  new scintHadronPhysics("pion"));
  //RegisterPhysics(  new G4EmStandardPhysics("hadron"));

  //Hadron Physics
  RegisterPhysics(  new HadronPhysicsQGSP_BERT_HP("hadron"));

  //Ion Physics
  RegisterPhysics(  new G4IonPhysics("ion"));
}

scintPhysicsList::~scintPhysicsList()
{
}

void scintPhysicsList::SetCuts(){
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}




