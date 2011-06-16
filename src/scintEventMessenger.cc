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
#include "scintEventMessenger.hh"
#include "scintEventAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintEventMessenger::scintEventMessenger(scintEventAction* event)
:scintEvent(event)
{
  saveThresholdCmd = new G4UIcmdWithAnInteger("/scint/saveThreshold",this);
  saveThresholdCmd->SetGuidance("Set the photon count threshold for saving the random number seed for an event.");
  saveThresholdCmd->SetParameterName("photons",true);
  saveThresholdCmd->SetDefaultValue(4500);
  saveThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/scint/eventVerbose",this);
  verboseCmd->SetGuidance("Set the verbosity of event data.");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);

  pmtThresholdCmd = new G4UIcmdWithAnInteger("/scint/pmtThreshold",this);
  pmtThresholdCmd->SetGuidance("Set the pmtThreshold (in # of photons)");

  forceDrawPhotonsCmd=new G4UIcmdWithABool("/scint/forceDrawPhotons",this);
  forceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
  forceDrawPhotonsCmd
    ->SetGuidance("(Higher priority than /scint/forceDrawNoPhotons)");

  forceDrawNoPhotonsCmd=new G4UIcmdWithABool("/scint/forceDrawNoPhotons",this);
  forceDrawNoPhotonsCmd->SetGuidance("Force no drawing of photons.");
  forceDrawNoPhotonsCmd
    ->SetGuidance("(Lower priority than /scint/forceDrawPhotons)");

  SaveHitsCmd = new G4UIcmdWithABool("/scint/saveHits",this);
  SaveHitsCmd->SetGuidance("Set flag to save hits in each run");
  SaveHitsCmd->SetGuidance("into file 'hits.out'");
  SaveHitsCmd->SetGuidance("Default = true");
  SaveHitsCmd->SetParameterName("saveHitsFlag", false);

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintEventMessenger::~scintEventMessenger(){
  delete saveThresholdCmd;
  delete verboseCmd;
  delete pmtThresholdCmd;
  delete forceDrawPhotonsCmd;
  delete forceDrawNoPhotonsCmd;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventMessenger::SetNewValue(G4UIcommand* command, G4String newValue){ 
  if( command == saveThresholdCmd ){ 
    scintEvent->SetSaveThreshold(saveThresholdCmd->GetNewIntValue(newValue));
  }
  else if( command == verboseCmd ){
    scintEvent->SetEventVerbose(verboseCmd->GetNewIntValue(newValue));
  }
  else if( command == pmtThresholdCmd ){
    scintEvent->SetPMTThreshold(pmtThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == forceDrawPhotonsCmd){
    scintEvent->SetForceDrawPhotons(forceDrawPhotonsCmd
				  ->GetNewBoolValue(newValue));
  }
  else if(command == forceDrawNoPhotonsCmd){
    scintEvent->SetForceDrawNoPhotons(forceDrawNoPhotonsCmd
				  ->GetNewBoolValue(newValue));
    G4cout<<"TEST"<<G4endl;
  }
}





