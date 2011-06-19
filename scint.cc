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
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4String.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactory.hh"

#include "scintDetectorConstruction.hh"
#include "scintPhysicsList.hh"
#include "scintPrimaryGeneratorAction.hh"
#include "scintEventAction.hh"
#include "scintStackingAction.hh"
#include "scintSteppingAction.hh"
#include "scintTrackingAction.hh"
#include "scintRunAction.hh"
#include "scintSteppingVerbose.hh"

#include "RecorderBase.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include <fstream>
std::ofstream hitsfile;
//std::ifstream timing;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
int main(int argc, char** argv)
{
    // TODO check: is this a stepping manager or a function to be carried out for stepping?

  G4VSteppingVerbose::SetInstance(new scintSteppingVerbose);

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new scintDetectorConstruction);
  runManager->SetUserInitialization(new scintPhysicsList);

  //G4PhysListFactory factory;
  //G4VModularPhysicsList* phys = factory.ReferencePhysList();
  //runManager->SetUserInitialization(phys);

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // TODO check doc: WTF is a recorder?
  RecorderBase* recorder = NULL;//No recording is done in this example

  runManager->SetUserAction(new scintPrimaryGeneratorAction);
  runManager->SetUserAction(new scintStackingAction);
  
  runManager->SetUserAction(new scintRunAction(recorder));
  runManager->SetUserAction(new scintEventAction(recorder));
  runManager->SetUserAction(new scintTrackingAction(recorder));
  runManager->SetUserAction(new scintSteppingAction(recorder));

  runManager->Initialize();
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if(argc==1){
    G4UIsession* session=0;
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);
#else
    session = new G4UIterminal();
#endif

    //execute vis.mac
    UI->ApplyCommand("/control/execute vis.mac");

    session->SessionStart();
    delete session;

  }
  else{
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UI->ApplyCommand(command+filename);
  }

  if(recorder)delete recorder;

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;
  return 0;
}


