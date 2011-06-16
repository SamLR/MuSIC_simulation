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
#include "scintSteppingAction.hh"
#include "scintEventAction.hh"
#include "scintTrackingAction.hh"
#include "scintTrajectory.hh"
#include "scintPMTSD.hh"
#include "scintUserTrackInformation.hh"
#include "scintUserEventInformation.hh"
#include "scintSteppingMessenger.hh"
#include "RecorderBase.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintSteppingAction::scintSteppingAction(RecorderBase* r)
  :recorder(r),oneStepPrimaries(false)
{
  steppingMessenger = new scintSteppingMessenger(this);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintSteppingAction::~scintSteppingAction()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintSteppingAction::UserSteppingAction(const G4Step * theStep){

  //G4cout << "UserSteppingAction" << G4endl;

  G4Track* theTrack = theStep->GetTrack();
  
  scintUserTrackInformation* trackInformation
    =(scintUserTrackInformation*)theTrack->GetUserInformation();
  scintUserEventInformation* eventInformation
    =(scintUserEventInformation*)G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetUserInformation();

  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundary=NULL;
  
  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm 
      = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }

  if(theTrack->GetParentID()==0){
    //This is a primary track
    
    G4TrackVector* fSecondary=fpSteppingManager->GetfSecondary();   
    G4int tN2ndariesTot = fpSteppingManager->GetfN2ndariesAtRestDoIt()
      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
      + fpSteppingManager->GetfN2ndariesPostStepDoIt();
    
    //If we havent already found the conversion position and there were 
    //secondaries generated, then search for it
    if(!eventInformation->IsConvPosSet() && tN2ndariesTot>0 ){    
      for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
	  lp1<(*fSecondary).size(); lp1++){
	const G4VProcess* creator=(*fSecondary)[lp1]->GetCreatorProcess();
	if(creator){
	  G4String creatorName=creator->GetProcessName();
	  if(creatorName=="phot"||creatorName=="compt"||creatorName=="conv"){ 
	    //since this is happening before the secondary is being tracked
	    //the Vertex position has not been set yet(set in initial step)
	    eventInformation->SetConvPos((*fSecondary)[lp1]->GetPosition());
	  } 
	}
      }
    }

    if(oneStepPrimaries&&thePrePV->GetName()=="scintillator")
      theTrack->SetTrackStatus(fStopAndKill);
  }

  if(!thePostPV){//out of world
    return;
  }

  //G4cout << "UserSteppingAction::TestBefore" << G4endl;

  G4ParticleDefinition* particleType = theTrack->GetDefinition();

  //G4cout << "Particle Type : " << particleType->GetParticleName() << G4endl;

  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    //G4cout << "UserSteppingAction::OpticalPhotonOnly" << G4endl;

    if(thePrePV->GetName()=="Slab")
      //force drawing of photons in WLS slab
      trackInformation->SetForceDrawTrajectory(true);
    else if(thePostPV->GetName()=="expHall")
      //Kill photons entering expHall from something other than Slab
      theTrack->SetTrackStatus(fStopAndKill);
        
    //Was the photon absorbed by the absorption process
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName()
       =="OpAbsorption"){
      eventInformation->IncAbsorption();
      trackInformation->AddTrackStatusFlag(absorbed);
    }
   
    boundaryStatus=boundary->GetStatus();
 
    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      switch(boundaryStatus){
      case Absorption:
	trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
	eventInformation->IncBoundaryAbsorption();
	break;
      case Detection: //Note, this assumes that the volume causing detection
                      //is the photocathode because it is the only one with
	              //non-zero efficiency
	{
	  //Trigger sensitive detector manually since photon is
	  //absorbed but status was Detection
	  G4SDManager* SDman = G4SDManager::GetSDMpointer();
	  G4String sdName="/scintDet/pmtSD";
	  scintPMTSD* pmtSD = (scintPMTSD*)SDman
	    ->FindSensitiveDetector(sdName);
	  if(pmtSD)
	    pmtSD->ProcessHits_constStep(theStep,NULL);
	  trackInformation->AddTrackStatusFlag(hitPMT);
	  break;
	}
      case FresnelReflection:
      case TotalInternalReflection:
      case SpikeReflection:
	trackInformation->IncReflections();
	break;
      default:
	break;
      }
      if(thePostPV->GetName()=="sphere")
	trackInformation->AddTrackStatusFlag(hitSphere);
    }
  }
  
  //G4cout << "UserSteppingAction::TestAfter" << G4endl;

  if(recorder)recorder->RecordStep(theStep);
}


















