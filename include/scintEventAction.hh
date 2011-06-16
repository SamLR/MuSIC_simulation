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

#ifndef scintEventAction_h
#define scintEventAction_h 1

#include "scintEventMessenger.hh"
#include "scintScintHit.hh"
#include "scintRunActionMessenger.hh"
#include "scintRunAction.hh"

#include "G4ios.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Event;
class RecorderBase;
class scintRunAction;
class scintPrimaryGeneratorAction;

class scintEventAction : public G4UserEventAction
{
public:
  scintEventAction(RecorderBase*); //hello git
  scintEventAction(RecorderBase*, scintRunAction*);
  ~scintEventAction();
  
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  void SetSaveThreshold(G4int save);

  void SetEventVerbose(G4int v){verbose=v;}

  void SetPMTThreshold(G4int t){pmtThreshold=t;}

  void SetForceDrawPhotons(G4bool b){forcedrawphotons=b;}
  void SetForceDrawNoPhotons(G4bool b){forcenophotons=b;}

  void SetSaveHitsFlag (G4int val)  {saveHitsFlag = val;};

  void writeScintHitsToFile(void);

  void writePmtHitsToFile(void);

  void writeEnergyToFile(void);

private:
  G4int event_id;
  G4int n_hit;
  G4int pmts;
  G4int totalHits;

  RecorderBase* recorder;
  scintEventMessenger* eventMessenger;

  G4int              saveThreshold;

  G4int              scintCollID;
  G4int              pmtCollID;

  G4int              verbose;
  
  G4int              pmtThreshold;

  G4int saveHitsFlag;
  
  G4bool forcedrawphotons;
  G4bool forcenophotons;

  G4String particleName;
  G4String firstParticleName;

  scintRunAction*    runAct;  //pointer to run action

  G4double firstScintHitTime;
  G4double timePMTHit;
  G4double aveTimePMTHit;
  G4double totEnergy;

  G4double particleEnergy;
  G4double leftAve;
  G4double rightAve;

  G4double one, two, three, four;

  G4double max_pmt0, max_pmt1, max_pmt2, max_pmt3;
  G4double min_pmt0, min_pmt1, min_pmt2, min_pmt3;

};

#endif

    
