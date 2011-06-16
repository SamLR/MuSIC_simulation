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

#include "scintRunActionMessenger.hh"
#include "scintRunAction.hh"

#include <sstream>

//#include "scintRunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

scintRunActionMessenger::scintRunActionMessenger(scintRunAction* run)
:scintRun(run)
{ 
  SaveHitsCmd = new G4UIcmdWithAString("/scint/hitsfile",this);
  SaveHitsCmd->SetGuidance("output file for hits collection (txt)");
  SaveHitsCmd->SetGuidance("Default = hits.out");
  SaveHitsCmd->SetParameterName("savehitsFile", false);
  SaveHitsCmd->SetDefaultValue("hits.out");
  /*
  SavePmtCmd = new G4UIcmdWithAString("/dmx/pmtfile",this);
  SavePmtCmd->SetGuidance("output file for pmt hits (txt)");
  SavePmtCmd->SetGuidance("Default = pmt.out");
  SavePmtCmd->SetParameterName("savepmtFile", false);
  SavePmtCmd->SetDefaultValue("pmt.out");

  SaveHistFileCmd = new G4UIcmdWithAString("/dmx/histogramfile",this);
  SaveHistFileCmd->SetGuidance("output file for histograms");
  SaveHistFileCmd->SetGuidance("Default = dmx.his");
  //  SaveHistFileCmd->SetParameterName("savehistFile", false);
  SaveHistFileCmd->SetParameterName("histFile", false);
  SaveHistFileCmd->SetDefaultValue("dmx.his");

  PlotEventCmd = new G4UIcmdWithABool("/dmx/PlotDuringRun",this);
  PlotEventCmd->SetGuidance("Flag for Interactive plotting during run");
  PlotEventCmd->SetGuidance("Default = false");
  PlotEventCmd->SetParameterName("plotevent", false);
  PlotEventCmd->SetDefaultValue(false);

  InteractPlotCmd = new G4UIcmdWithABool("/dmx/PlotAtEnd",this);
  InteractPlotCmd->SetGuidance("Flag for Interactive plotting at end of run");
  InteractPlotCmd->SetGuidance("Default = false");
  InteractPlotCmd->SetParameterName("interactplot", false);
  InteractPlotCmd->SetDefaultValue(false);
  */
  //  FileCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

scintRunActionMessenger::~scintRunActionMessenger()
{
  delete SaveHitsCmd;  
  //delete SavePmtCmd;  
  //delete SaveHistFileCmd;  
  //delete PlotEventCmd;
  //delete InteractPlotCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void scintRunActionMessenger::SetNewValue(G4UIcommand * command, G4String val)
{ 
  if(command == SaveHitsCmd)
    scintRun->SetsavehitsFile(val);

  /*
  if(command == SavePmtCmd)
    scintRun->SetsavepmtFile(newValue);

  if(command == SaveHistFileCmd)
    scintRun->SetsavehistFile(newValue);

  if(command == PlotEventCmd) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    scintRun->Setplotevent(vl!=0);
  }

  if(command == InteractPlotCmd) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    scintRun->Setinteractplot(vl!=0);
  }
  */
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





