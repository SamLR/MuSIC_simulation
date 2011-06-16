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
#include "scintHadronPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


scintHadronPhysics::scintHadronPhysics(const G4String& name)
                   :  G4VPhysicsConstructor(name)
{
}

scintHadronPhysics::~scintHadronPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

void scintHadronPhysics::ConstructParticle()
{
  // Pi
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
}

#include "G4ProcessManager.hh"

void scintHadronPhysics::ConstructProcess()
{
  fPiPlusIonisation = new G4hIonisation();
  fPiPlusMultipleScattering = new G4MultipleScattering();
  fPionPlusInelasticProcess = new G4PionPlusInelasticProcess();
  fPiPlusLEPModel = new G4LEPionPlusInelastic();
  fPiPlusHEPModel = new G4HEPionPlusInelastic();

  fPiMinusIonisation = new G4hIonisation();
  fPiMinusMultipleScattering = new G4MultipleScattering;
  fPionMinusInelasticProcess = new G4PionMinusInelasticProcess();
  fPiMinusLEPModel = new G4LEPionMinusInelastic();
  fPiMinusHEPModel = new G4HEPionMinusInelastic();

  G4ProcessManager * pManager = 0;

  // Pion Plus Physics
  pManager = G4PionPlus::PionPlus()->GetProcessManager();

  fPionPlusInelasticProcess->RegisterMe(fPiPlusLEPModel);
  fPionPlusInelasticProcess->RegisterMe(fPiPlusHEPModel);
  pManager->AddDiscreteProcess(fPionPlusInelasticProcess);
   
  pManager->AddProcess(fPiPlusMultipleScattering,-1,  1, 1);
  pManager->AddProcess(fPiPlusIonisation,        -1,  2, 2);

  // Pion Minus Physics
  pManager = G4PionMinus::PionMinus()->GetProcessManager();

  fPionMinusInelasticProcess->RegisterMe(fPiMinusLEPModel);
  fPionMinusInelasticProcess->RegisterMe(fPiMinusHEPModel);
  pManager->AddDiscreteProcess(fPionMinusInelasticProcess);
   
  pManager->AddProcess(fPiMinusMultipleScattering,-1,  1, 1);
  pManager->AddProcess(fPiMinusIonisation,        -1,  2, 2);

  // Proton Physics
  pManager = G4Proton::Proton()->GetProcessManager();

  fProtonLEPModel = new G4LEProtonInelastic();
  fProtonHEPModel = new G4HEProtonInelastic();
  fProtonInelasticProcess = new G4ProtonInelasticProcess();
  fProtonIonisation = new G4hIonisation();
  fProtonMultipleScattering = new G4MultipleScattering;

  fProtonInelasticProcess->RegisterMe(fProtonLEPModel);
  fProtonInelasticProcess->RegisterMe(fProtonHEPModel);
  pManager->AddDiscreteProcess(fProtonInelasticProcess);
  
  pManager->AddProcess(fProtonMultipleScattering,-1,  1, 1);
  pManager->AddProcess(fProtonIonisation,        -1,  2, 2);
  
}



