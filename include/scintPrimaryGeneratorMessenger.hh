//-*-C++-*-
// $Id: scintPrimaryGeneratorMessenger.hh 189 2008-10-31 09:24:18Z comet $

#ifndef scintPrimaryGeneratorMessenger_HH
#define scintPrimaryGeneratorMessenger_HH 1

static const char* scintPrimaryGeneratorMessenger_hh =
"@(#) $Id: scintPrimaryGeneratorMessenger.hh 189 2008-10-31 09:24:18Z comet $";

#include "G4UImessenger.hh"
#include "globals.hh"

class scintPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;



class scintPrimaryGeneratorMessenger: public G4UImessenger
{
public:

  scintPrimaryGeneratorMessenger(scintPrimaryGeneratorAction*);
  ~scintPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  scintPrimaryGeneratorAction* fMyAction; 

  G4UIdirectory*      fModeDir;

  G4UIcmdWithAString* fFileOpenCmd;
  G4UIcmdWithoutParameter* fFileCloseCmd;
  G4UIcmdWithAString* fMPosCmd;
  G4UIcmdWithAString* fMEngCmd;

  G4UIdirectory*      fSpreadDir;

  G4UIcmdWithADoubleAndUnit* fDEngCmd;
  G4UIcmdWith3VectorAndUnit* fDPosCmd;
  G4UIcmdWithADoubleAndUnit* fDDirectionCmd;
  G4UIcmdWithADouble* fDDispersionCmd;
  G4UIcmdWithADoubleAndUnit* fDCutOffEnergyCmd;


};

#endif
