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

#include "scintRunAction.hh"
#include "scintRunActionMessenger.hh"
#include "RecorderBase.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

extern   std::ofstream hitsfile;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintRunAction::scintRunAction(RecorderBase* r)
//scintRunAction::scintRunAction()
//  :recorder(r)
{
  recorder = r;
  runMessenger = new scintRunActionMessenger(this);
  savehitsFile = "hits.out";
  //savepmtFile  = "pmt.out";
  //savehistFile = "dmx.his";
  //plotevent    = false;
  //interactplot = false;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintRunAction::~scintRunAction()
{
  runMessenger = 0;
  delete runMessenger;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintRunAction::BeginOfRunAction(const G4Run* aRun){
  if(recorder)recorder->RecordBeginOfRun(aRun);

  hitsfile.open("timing.out", std::ios::out);
  if (! hitsfile.good()) {
    G4String errorMessage="*** fail to open a file (timing.out).";
    G4Exception(errorMessage);  
  }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintRunAction::EndOfRunAction(const G4Run* aRun){
  if(recorder)recorder->RecordEndOfRun(aRun);
  hitsfile.close();
}
