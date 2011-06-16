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
#include "scintPMTSD.hh"
#include "scintPMTHit.hh"
#include "scintDetectorConstruction.hh"
#include "scintUserTrackInformation.hh"
#include "scintEventAction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

extern  std::ofstream hitsfile;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintPMTSD::scintPMTSD(G4String name)
  :G4VSensitiveDetector(name),pmtHitCollection(0),pmtPositionsX(0)
  ,pmtPositionsY(0),pmtPositionsZ(0)
{
  collectionName.insert("pmtHitCollection");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintPMTSD::~scintPMTSD()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintPMTSD::Initialize(G4HCofThisEvent* HCE){
  pmtHitCollection = new scintPMTHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  //Store collection with event and keep ID
  static G4int HCID = -1;
  if(HCID<0){ 
    HCID = GetCollectionID(0); 
  }
  HCE->AddHitsCollection( HCID, pmtHitCollection );
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool scintPMTSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist){
  return ProcessHits_constStep(aStep,ROhist);
}

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool scintPMTSD::ProcessHits_constStep(const G4Step* aStep,
					 G4TouchableHistory*){

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() 
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
 
  //User replica number 1 since photocathode is a daughter volume
  //to the pmt which was replicated
  G4int pmtNumber=
    aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  G4VPhysicalVolume* physVol=
    aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);
  G4int event_id=
    G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetEventID();
    
  //G4cout << "Copy number : " << currentCopyNo << G4endl;
  //G4cout << "PMT number : " << pmtNumber << G4endl;
  //G4cout << "Event ID : " << event_id << G4endl;

  TSlice = aStep->GetPostStepPoint()->GetGlobalTime();

  //G4String* processName=
  //aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(1);

  //Find the correct hit collection
  G4int n=pmtHitCollection->entries();
  //G4cout << "HitCollection Entries : " << n << G4endl;
  scintPMTHit* hit=NULL;
  for(G4int i=0;i<n;i++){
    if((*pmtHitCollection)[i]->GetPMTNumber()==pmtNumber){
      hit=(*pmtHitCollection)[i];
      break;
    }
  }
  
  //G4cout << pmtNumber << G4endl;
  
  if(hit==NULL){//this pmt wasnt previously hit in this event
    hit = new scintPMTHit(); //so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetPMTPhysVol(physVol);
    pmtHitCollection->insert(hit);
    hit->SetPMTPos((*pmtPositionsX)[pmtNumber],(*pmtPositionsY)[pmtNumber],
		   (*pmtPositionsZ)[pmtNumber]);
  }

  hit->IncPhotonCount(); //increment hit for the selected pmt
    
  if(!scintDetectorConstruction::GetSphereOn()){
    hit->SetDrawit(true);
    //If the sphere is disabled then this hit is automaticaly drawn
  }
  else{//sphere enabled
    scintUserTrackInformation* trackInfo=
      (scintUserTrackInformation*)aStep->GetTrack()->GetUserInformation();
    if(trackInfo->GetTrackStatus()&hitSphere) 
      //only draw this hit if the photon has hit the sphere first
      hit->SetDrawit(true);
  }

  //G4cout << n << "    " << TSlice
  
  if(hitsfile.is_open()) {
  
    hitsfile << std::setiosflags(std::ios::fixed)
	     << std::setprecision(4)
	     << std::setiosflags(std::ios::left)
	     << std::setw(6)
	     << event_id << "\t    "
	     << pmtNumber << "\t    "
	     << TSlice/nanosecond << "\r"
	     << G4endl;
  }
   
  return true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintPMTSD::EndOfEvent(G4HCofThisEvent* ){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintPMTSD::clear(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintPMTSD::DrawAll(){
} 

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintPMTSD::PrintAll(){
} 

