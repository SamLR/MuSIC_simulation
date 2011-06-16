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

#ifndef scintDetectorMessenger_h
#define scintDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class scintDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;

class scintDetectorMessenger: public G4UImessenger
{
public:
  scintDetectorMessenger(scintDetectorConstruction*);
  ~scintDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  scintDetectorConstruction*   scintDetector;
  G4UIdirectory*               detectorDir;
  G4UIdirectory*               volumesDir;
  G4UIcmdWith3VectorAndUnit*   dimensionsCmd;
  G4UIcmdWithADoubleAndUnit*   housingThicknessCmd;
  G4UIcmdWithADoubleAndUnit*   pmtRadiusCmd;
  G4UIcmdWithAnInteger*        nxCmd;
  G4UIcmdWithAnInteger*        nyCmd;
  G4UIcmdWithAnInteger*        nzCmd;
  G4UIcmdWithABool*            sphereCmd;
  G4UIcmdWithADouble*          reflectivityCmd;
  G4UIcmdWithABool*            wlsCmd;
  G4UIcmdWithABool*            lxeCmd;
  G4UIcmdWithAnInteger*        nFibersCmd;
  G4UIcommand*                 updateCmd;
  G4UIcommand*                 defaultsCmd;
  G4UIcmdWithADouble*          MainScintYield;
  G4UIcmdWithADouble*          WLSScintYield;
};

#endif

