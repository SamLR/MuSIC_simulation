//-*-C++-*-
// $Id: scintPrimaryGeneratorMessenger.cc 228 2008-12-24 08:20:56Z comet $
//

static const char* scintPrimaryGeneratorMessenger_cc = "@(#) $Id: scintPrimaryGeneratorMessenger.cc 228 2008-12-24 08:20:56Z comet $";

#include "scintPrimaryGeneratorMessenger.hh"

#include "scintPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"



scintPrimaryGeneratorMessenger::
scintPrimaryGeneratorMessenger(scintPrimaryGeneratorAction* myGun)
  :fMyAction(myGun)
{ 

  G4cout << "Calling scintPrimaryGeneratorMessenger" << G4endl;
  //
  // Command sub Tree directory of the Gun Mode settings.
  //      
  fModeDir = new G4UIdirectory("/gun/mode/");
  fModeDir->SetGuidance("Gun mode control.");

  fFileOpenCmd = new G4UIcmdWithAString("/gun/mode/fileOpen",this);
  fFileOpenCmd->SetGuidance("Open an external file of primary guns.");
  fFileOpenCmd->SetParameterName("fileName",true);
  fFileOpenCmd->SetDefaultValue("default.gun");
  fFileOpenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  
  fFileCloseCmd = new G4UIcmdWithoutParameter("/gun/mode/fileClose",this);
  fFileCloseCmd->SetGuidance("Close the external file.");
  fFileCloseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMPosCmd = new G4UIcmdWithAString("/gun/mode/position",this);
  fMPosCmd->SetGuidance("Set the position mode.");
  fMPosCmd->SetGuidance("  Choice : point, random, target");
  fMPosCmd->SetGuidance("    point  -> point source,");
  fMPosCmd->SetGuidance("    random -> randomly spread around the point source,");
  fMPosCmd->SetGuidance("    target -> uniformly distributed in target disks.");
  fMPosCmd->SetParameterName("choice",true);
  fMPosCmd->SetDefaultValue("point");
  fMPosCmd->SetCandidates("point random target");
  fMPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMEngCmd = new G4UIcmdWithAString("/gun/mode/energy",this);
  fMEngCmd->SetGuidance("Set the energy mode.");
  fMEngCmd->SetGuidance("  Choice : mono, random, dio, rpc, pmc, dispersivebeam, cosmic");
  fMEngCmd->SetGuidance("    mono   -> monochromatic,");
  fMEngCmd->SetGuidance("    random -> randomly spread around the specified energy,");
  fMEngCmd->SetGuidance("    dio    -> Watanabe-Shanker's DIO spectrum.");
  fMEngCmd->SetGuidance("              Production threshold energy should be defined by /gun/energy beforehand.");
  fMEngCmd->SetGuidance("    rpc    -> Radiative Pion Capture spectrum (photon).");
  fMEngCmd->SetGuidance("    pmc    -> Muon-Captured Proton spectrum.");
  fMEngCmd->SetGuidance("    dispersivebeam -> Dispersion along Y.");
  fMEngCmd->SetGuidance("    cosmic -> Cosmic Muon Spectrum.");
  fMEngCmd->SetGuidance("              Both energy and angle will be generated.");
  fMEngCmd->SetGuidance("              Maximum \"momentum\" should be defined by /gun/energy beforehand.");
  fMEngCmd->SetParameterName("choice",true);
  fMEngCmd->SetDefaultValue("mono");
  fMEngCmd->SetCandidates("mono random dio rpc pmc dispersivebeam cosmic");
  fMEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //
  // Command sub Tree directory of the Spreads of guns.
  //

  fSpreadDir = new G4UIdirectory("/gun/spread/");
  fSpreadDir->SetGuidance("Spreads of gun.");

  fDEngCmd = new G4UIcmdWithADoubleAndUnit("/gun/spread/energy",this);
  fDEngCmd->SetGuidance("Set the kinetic energy spread.");
  fDEngCmd->SetParameterName("engSpread",false);
  fDEngCmd->SetRange("engSpread>=0.");
  fDEngCmd->SetUnitCategory("Energy");
  fDEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDPosCmd = new G4UIcmdWith3VectorAndUnit("/gun/spread/position",this);
  fDPosCmd->SetGuidance("Set position spread.");
  fDPosCmd->SetParameterName("X0","Y0","Z0",true,true);
  fDPosCmd->SetRange("X0>=0. && Y0>=0 && Z0>=0");
  fDPosCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fDPosCmd->SetUnitCategory("Length");
  fDPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDDirectionCmd = new G4UIcmdWithADoubleAndUnit("/gun/spread/direction",this);
  fDDirectionCmd->SetGuidance("Set direction spread.");
  fDDirectionCmd->SetParameterName("DPolar",false);
  fDDirectionCmd->SetDefaultValue(0.0);
  fDDirectionCmd->SetUnitCategory("Angle");
  fDDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDDispersionCmd = new G4UIcmdWithADouble("/gun/spread/dispersion",this);
  fDDispersionCmd->SetGuidance("Set the beam dispersion along y (MeV/cm).");
  fDDispersionCmd->SetParameterName("dispSpread",false);
  fDDispersionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDCutOffEnergyCmd =
    new G4UIcmdWithADoubleAndUnit("/gun/spread/cutoffenergy",this);
  fDCutOffEnergyCmd->SetGuidance("Set the cut-off kinetic energy.");
  fDCutOffEnergyCmd->SetParameterName("cutoffeng",false);
  fDCutOffEnergyCmd->SetRange("cutoffeng>=0.");
  fDCutOffEnergyCmd->SetUnitCategory("Energy");
  fDCutOffEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

scintPrimaryGeneratorMessenger::~scintPrimaryGeneratorMessenger()
{
  delete fDCutOffEnergyCmd;
  delete fDDispersionCmd;
  delete fDDirectionCmd;
  delete fDPosCmd;
  delete fDEngCmd;
  delete fSpreadDir;

  delete fMEngCmd;
  delete fMPosCmd;
  delete fFileCloseCmd;
  delete fFileOpenCmd;
  delete fModeDir;

}

void scintPrimaryGeneratorMessenger::
SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fFileOpenCmd )
    { fMyAction->FileOpen(newValue); }
  
  if( command == fFileCloseCmd )
    { fMyAction->FileClose(); }
  
  if( command == fMPosCmd )
    { fMyAction->SetPositionMode(newValue);}

  if( command == fMEngCmd )
    { fMyAction->SetEnergyMode(newValue);}

  if( command == fDEngCmd ) {
    fMyAction->SetEnergySpread(fDEngCmd->GetNewDoubleValue(newValue));
  }

  if( command == fDPosCmd ) {
    fMyAction->SetPositionSpread(fDPosCmd->GetNew3VectorValue(newValue));
  }

  if( command == fDDirectionCmd ) {
    fMyAction->
      SetDirectionSpread(fDDirectionCmd->GetNewDoubleValue(newValue));
  }

  if( command == fDDispersionCmd ) {
    fMyAction->
      SetBeamDispersion(fDDispersionCmd->GetNewDoubleValue(newValue));
  }

  if( command == fDCutOffEnergyCmd ) {
    fMyAction->
      SetCutOffEnergy(fDCutOffEnergyCmd->GetNewDoubleValue(newValue));
  }
  
}
