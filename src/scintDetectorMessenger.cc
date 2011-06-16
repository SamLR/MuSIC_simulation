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
#include "scintDetectorMessenger.hh"
#include "scintDetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintDetectorMessenger::scintDetectorMessenger(scintDetectorConstruction* scintDetect)
:scintDetector(scintDetect)
{
  //Setup a command directory for detector controls with guidance
  detectorDir = new G4UIdirectory("/scint/detector/");
  detectorDir->SetGuidance("Detector geometry control");

  volumesDir = new G4UIdirectory("/scint/detector/volumes/");
  volumesDir->SetGuidance("Enable/disable volumes");
   
  //Various commands for modifying detector geometry
  dimensionsCmd = 
    new G4UIcmdWith3VectorAndUnit("/scint/detector/dimensions",this);
  dimensionsCmd->SetGuidance("Set the dimensions of the detector volume.");
  dimensionsCmd->SetParameterName("scint_x","scint_y","scint_z",false);
  dimensionsCmd->SetDefaultUnit("cm");

  housingThicknessCmd = new G4UIcmdWithADoubleAndUnit
    ("/scint/detector/housingThickness",this);
  housingThicknessCmd->SetGuidance("Set the thickness of the housing.");
  housingThicknessCmd->SetParameterName("d_mtl",false);
  housingThicknessCmd->SetDefaultUnit("cm");

  pmtRadiusCmd = new G4UIcmdWithADoubleAndUnit
    ("/scint/detector/pmtRadius",this);
  pmtRadiusCmd->SetGuidance("Set the radius of the PMTs.");
  pmtRadiusCmd->SetParameterName("radius",false);
  pmtRadiusCmd->SetDefaultUnit("cm");

  nxCmd = new G4UIcmdWithAnInteger("/scint/detector/nx",this);
  nxCmd->SetGuidance("Set the number of PMTs along the x-dimension.");
  nxCmd->SetParameterName("nx",false);

  nyCmd = new G4UIcmdWithAnInteger("/scint/detector/ny",this);
  nyCmd->SetGuidance("Set the number of PMTs along the y-dimension.");
  nyCmd->SetParameterName("ny",false);
  
  nzCmd = new G4UIcmdWithAnInteger("/scint/detector/nz",this);
  nzCmd->SetGuidance("Set the number of PMTs along the z-dimension.");
  nzCmd->SetParameterName("nz",false);

  sphereCmd = new G4UIcmdWithABool("/scint/detector/volumes/sphere",this);
  sphereCmd->SetGuidance("Enable/Disable the sphere."); 

  reflectivityCmd = new G4UIcmdWithADouble("/scint/detector/reflectivity",this);
  reflectivityCmd->SetGuidance("Set the reflectivity of the housing.");

  wlsCmd = new G4UIcmdWithABool("/scint/detector/volumes/wls",this);
  wlsCmd->SetGuidance("Enable/Disable the WLS slab");

  lxeCmd = new G4UIcmdWithABool("/scint/detector/volumes/lxe",this);
  wlsCmd->SetGuidance("Enable/Disable the main detector volume.");

  nFibersCmd = new G4UIcmdWithAnInteger("/scint/detector/nfibers",this);
  nFibersCmd->SetGuidance("Set the number of WLS fibers in the WLS slab.");

  updateCmd = new G4UIcommand("/scint/detector/update",this);
  updateCmd->SetGuidance("Update the detector geometry with changed values.");
  updateCmd->SetGuidance
    ("Must be run before beamOn if detector has been changed.");
  
  defaultsCmd = new G4UIcommand("/scint/detector/defaults",this);
  defaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  defaultsCmd->SetGuidance("(Update still required)");

  MainScintYield=new G4UIcmdWithADouble("/scint/detector/MainScintYield",this);
  MainScintYield->SetGuidance("Set scinitillation yield of main volume.");
  MainScintYield->SetGuidance("Specified in photons/MeV");

  WLSScintYield = new G4UIcmdWithADouble("/scint/detector/WLSScintYield",this);
  WLSScintYield->SetGuidance("Set scintillation yield of WLS Slab");
  WLSScintYield->SetGuidance("Specified in photons/MeV");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintDetectorMessenger::~scintDetectorMessenger()
{
  delete dimensionsCmd;
  delete housingThicknessCmd;
  delete pmtRadiusCmd;
  delete nxCmd;
  delete nyCmd;
  delete nzCmd;
  delete updateCmd;
  delete detectorDir;
  delete volumesDir;
  delete defaultsCmd;
  delete sphereCmd;
  delete wlsCmd;
  delete lxeCmd;
  delete nFibersCmd;
  delete reflectivityCmd;
  delete MainScintYield;
  delete WLSScintYield;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == dimensionsCmd ){ 
    scintDetector->SetDimensions(dimensionsCmd->GetNew3VectorValue(newValue));
  }
  else if (command == housingThicknessCmd){
    scintDetector->SetHousingThickness(housingThicknessCmd
				     ->GetNewDoubleValue(newValue));
  }
  else if (command == pmtRadiusCmd){
    scintDetector->SetPMTRadius(pmtRadiusCmd->GetNewDoubleValue(newValue));
  }
  else if (command == nxCmd){
    scintDetector->SetNX(nxCmd->GetNewIntValue(newValue));
  }
  else if (command == nyCmd){
    scintDetector->SetNY(nyCmd->GetNewIntValue(newValue));
  }
  else if (command == nzCmd){
    scintDetector->SetNZ(nzCmd->GetNewIntValue(newValue));
  }
  else if (command == updateCmd){
    scintDetector->UpdateGeometry();
  }
  else if (command == defaultsCmd){
    scintDetector->SetDefaults();
  }
  else if (command == sphereCmd){
    scintDetector->SetSphereOn(sphereCmd->GetNewBoolValue(newValue));
  }
  else if (command == reflectivityCmd){
    scintDetector
      ->SetHousingReflectivity(reflectivityCmd->GetNewDoubleValue(newValue));
  }
  else if (command == wlsCmd){
    scintDetector->SetWLSSlabOn(wlsCmd->GetNewBoolValue(newValue));
  }
  else if (command == lxeCmd){
    scintDetector->SetMainVolumeOn(lxeCmd->GetNewBoolValue(newValue));
  }
  else if (command == nFibersCmd){
    scintDetector->SetNFibers(nFibersCmd->GetNewIntValue(newValue));
  }
  else if (command == MainScintYield){
   scintDetector->SetMainScintYield(MainScintYield->GetNewDoubleValue(newValue));
  }
  else if (command == WLSScintYield){
    scintDetector->SetWLSScintYield(WLSScintYield->GetNewDoubleValue(newValue));
  }
}



