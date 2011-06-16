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
#include "scintEventAction.hh"
#include "scintRunAction.hh"
#include "scintScintHit.hh"
#include "scintPrimaryGeneratorAction.hh"
#include "scintPMTHit.hh"
#include "scintUserEventInformation.hh"
#include "scintTrajectory.hh"
#include "RecorderBase.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

#include <fstream>
#include <iomanip>
#include <vector>

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintEventAction::scintEventAction(RecorderBase* r)//, scintRunAction* scintRun)
  :recorder(r),saveThreshold(0),scintCollID(-1),pmtCollID(-1),verbose(0),
   pmtThreshold(1),forcedrawphotons(false),forcenophotons(false)//,runAct(scintRun)
{
  //create a messenger
  eventMessenger=new scintEventMessenger(this);

  //defaults for messenger
  //savePmtFlag  = 0;
  saveHitsFlag = 1;

}
 
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintEventAction::~scintEventAction(){}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventAction::BeginOfEventAction(const G4Event* anEvent){
  
  //New event, add the user information object
  G4EventManager::
    GetEventManager()->SetUserInformation(new scintUserEventInformation);
  
  //G4int event_id = anEvent->GetEventID();

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(scintCollID<0)
    scintCollID=SDman->GetCollectionID("scintCollection");
  if(pmtCollID<0)
    pmtCollID=SDman->GetCollectionID("pmtHitCollection");

  if(recorder)recorder->RecordBeginOfEvent(anEvent);
}
 
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventAction::EndOfEventAction(const G4Event* anEvent){
  
  event_id = anEvent->GetEventID();

  firstScintHitTime = 0.;
  //G4double ScintHitTime      = 0.;
  timePMTHit        = 0.;
  aveTimePMTHit     = 0.;
  totEnergy         = 0.;
  //firstParticleE    = 0.;
  particleEnergy    = 0.;
  leftAve = 0;
  rightAve = 0;

  firstParticleName = "";
  particleName = "";

  scintUserEventInformation* eventInformation
    =(scintUserEventInformation*)anEvent->GetUserInformation();
  
  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
  
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
 
  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance()){
    for (G4int i=0; i<n_trajectories; i++){ 
      scintTrajectory* trj = (scintTrajectory*)
	((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
	trj->SetForceDrawTrajectory(forcedrawphotons);
	trj->SetForceNoDrawTrajectory(forcenophotons);
      }
      trj->DrawTrajectory(50);
    }
  }

  scintScintHitsCollection* SHC = 0;
  scintPMTHitsCollection* PHC = 0;
  G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
  
  //Get the hit collections
  if(HCE){
    if(scintCollID>=0) SHC = (scintScintHitsCollection*)(HCE->GetHC(scintCollID));
    if(pmtCollID>=0) PHC = (scintPMTHitsCollection*)(HCE->GetHC(pmtCollID));
  }

  //G4cout << "-------------------------------------------" << G4endl;
  
  //Hits in scintillator
  if(SHC){
    n_hit = SHC->entries();
    G4ThreeVector eWeightPos(0.);     
    G4double edep;
    G4double edepMax=0;

    //G4cout << "\tSHC Entries : " << n_hit << G4endl;
    G4cout << "-------------------------------------------" << G4endl;

    for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
      if(i==0) {
	firstScintHitTime = (*SHC)[0]->GetTime(); //get time of first scintillator hit
	//firstParticleE = (*SHC)[0]->GetParticleEnergy();
	firstParticleName = (*SHC)[0]->GetParticle();
      }

      //ScintHitTime =  (*SHC)[i]->GetTime()/nanosecond;
      //G4cout << "Scintillation Time : " << ScintHitTime << "ns" <<  G4endl;

      //if(PHC) {
      //G4int pmtNumber=PHC->entries();
      //	for(G4int j=0;j<pmtNumber;j++){
      //  G4double timePMTHit_scint = ( (*PHC)[i]->GetTime() - ScintHitTime );
	  //G4cout << "PMT Hit Time : " << timePMTHit_scint << G4endl;
      //}
      //}
      
      edep=(*SHC)[i]->GetEdep();
      totEnergy += edep;
      eventInformation->IncEDep(edep); //sum up the edep
      eWeightPos += (*SHC)[i]->GetPos()*edep;//calculate energy weighted pos
      particleName = (*SHC)[i]->GetParticle();
      particleEnergy = (*SHC)[i]->GetParticleEnergy();
      
      G4cout << "Particle Name : " << particleName << G4endl;
      G4cout << "Particle Energy : " << particleEnergy/MeV << "MeV" <<  G4endl;
      G4cout << "Deposited Energy : " << edep/MeV << "MeV" <<  G4endl;
      G4cout << "Total Deposited Energy : " << totEnergy/MeV << "MeV" <<  G4endl;
      G4cout << "-------------------------------------------" << G4endl;
      
      if(edep>edepMax){
	edepMax=edep;//store max energy deposit
	G4ThreeVector posMax=(*SHC)[i]->GetPos();
	eventInformation->SetPosMax(posMax,edep);
      }
    }
    if(eventInformation->GetEDep()==0.){
      if(verbose>0)G4cout<<"No hits in the scintillator this event."<<G4endl;
    }    
    else{
      //Finish calculation of energy weighted position
      eWeightPos/=eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
      if(verbose>0){
	G4cout << "\tEnergy weighted position of hits in scint : "
	       << eWeightPos/mm << G4endl;
      }
    }
    if(verbose>0){
    G4cout << "\tTotal energy deposition in scintillator : "
	   << eventInformation->GetEDep() / keV << " (keV)" << G4endl;
    }
  }

  if(PHC){
    G4ThreeVector reconPos=0.;
    pmts=PHC->entries();
    G4double left = 0;
    G4double right = 0;

    one = 0;
    two = 0;
    three = 0;
    four = 0;

    //G4cout << "\tPHC Entries : " << pmts << G4endl;

    //G4double hitCount = eventInformation->GetHitCount();
    //Gather info from all PMTs

    for(G4int i=0;i<pmts;i++){
      timePMTHit = ( (*PHC)[i]->GetTime() - firstScintHitTime );
      aveTimePMTHit += timePMTHit/(G4double)pmts;

      eventInformation->IncHitCount((*PHC)[i]->GetPhotonCount());
      reconPos+=(*PHC)[i]->GetPMTPos()*(*PHC)[i]->GetPhotonCount();

      //G4cout << "Time of the hit : " << timePMTHit << G4endl;

      if(((*PHC)[i]->GetPMTPos()).x()/cm < 0) { 
	//G4cout << "\tTime of hit at MPPC("<<i<<") " 
	//     << "(-ve,  "<<((*PHC)[i]->GetPMTPos()).y()/cm<<") : "
	//     << timePMTHit/nanosecond << "ns" << G4endl;
	//	G4cout << "\tHits at MPPC("<<i<<")"
	//     << "(-ve,  "<<((*PHC)[i]->GetPMTPos()).y()/cm<<") : "
	//     << (*PHC)[i]->GetPhotonCount() << G4endl;
	rightAve += timePMTHit/(G4double)pmts;
      }

      if(((*PHC)[i]->GetPMTPos()).x()/cm > 0) { 
	//G4cout << "\tTime of hit at MPPC("<<i<<") " 
	//     << "(+ve, "<<((*PHC)[i]->GetPMTPos()).y()/cm<<") : "
	//     << timePMTHit/nanosecond << "ns" << G4endl;

	//	G4cout << "\tHits at MPPC("<<i<<")"
	//     << "(+ve, "<<((*PHC)[i]->GetPMTPos()).y()/cm<<") : "
	//     << (*PHC)[i]->GetPhotonCount() << G4endl;
	leftAve += timePMTHit/(G4double)pmts;
      }
      
      if((*PHC)[i]->GetPhotonCount()>=pmtThreshold){
	eventInformation->IncPMTSAboveThreshold();
      }
      else{//wasnt above the threshold, turn it back off
	(*PHC)[i]->SetDrawit(false);
      }
      if(((*PHC)[i]->GetPMTPos()).x()/cm < 0) {   //if loop to get the PMT hits on the right with -ve position
	right += (*PHC)[i]->GetPhotonCount();
	//G4cout << "right : " << right << G4endl;	
      }
      if(((*PHC)[i]->GetPMTPos()).x()/cm > 0) {   //if loop to get the PMT hits on the left with +ve position
	left += (*PHC)[i]->GetPhotonCount();
	// G4cout << "left : " << left << G4endl;
      }

      if(((*PHC)[i]->GetPMTPos()).x()/cm < 0 && ((*PHC)[i]->GetPMTPos()).y()/cm < 0) {   //if loop to get the PMT hits on the right with -ve position
	one += (*PHC)[i]->GetPhotonCount();
	//G4cout << "right : " << right << G4endl;	
      }
      if(((*PHC)[i]->GetPMTPos()).x()/cm < 0 && ((*PHC)[i]->GetPMTPos()).y()/cm > 0) {   //if loop to get the PMT hits on the right with -ve position
	two += (*PHC)[i]->GetPhotonCount();
	//G4cout << "right : " << right << G4endl;	
      }
      if(((*PHC)[i]->GetPMTPos()).x()/cm > 0 && ((*PHC)[i]->GetPMTPos()).y()/cm < 0) {   //if loop to get the PMT hits on the right with -ve position
	three += (*PHC)[i]->GetPhotonCount();
	//G4cout << "right : " << right << G4endl;	
      }
      if(((*PHC)[i]->GetPMTPos()).x()/cm > 0 && ((*PHC)[i]->GetPMTPos()).y()/cm > 0) {   //if loop to get the PMT hits on the right with -ve position
	four += (*PHC)[i]->GetPhotonCount();
	//G4cout << "right : " << right << G4endl;	
      }
    }

    //G4cout << "Average Time Hit: " << aveTimePMTHit << G4endl;

    /*
    G4double ratio1 = left/(left+right);
    G4double pos1 = 20.0 - ratio1*40.0;

    G4double ratio2 = leftAve/(leftAve+rightAve);
    G4double pos2 = 20.0 - ratio2*40.0;

    G4cout << "\tPosition of incidence from position method : " << pos1 << "cm" << G4endl;
    G4cout << "\tPosition of incidence from time method     : " << pos2 << "cm" << G4endl;
    */
    
    if(eventInformation->GetHitCount()>0){//dont bother unless there were hits
      reconPos/=eventInformation->GetHitCount();
      if(verbose>0){
	G4cout << "\tReconstructed position of hits in scint : "
	       << reconPos/mm << G4endl;
      }
      eventInformation->SetReconPos(reconPos);
    }
    PHC->DrawAllHits();
  }
  
  totalHits = eventInformation->GetHitCount();

  //G4cout << "Total number of photons hitting the MPPCs : " << totalHits << G4endl;

  if(verbose>0 | verbose==0){
    //End of event output. later to be controlled by a verbose level
    /*
    G4cout << "\tNumber of photons that hit MPPCs in this event : "
	   << eventInformation->GetHitCount() << G4endl;
    G4cout << "\tNumber of MPPCs above threshold("<<pmtThreshold<<") : "
	   << eventInformation->GetPMTSAboveThreshold() << G4endl;
    G4cout << "\tNumber of photons produced by scintillation in this event : "
	   << eventInformation->GetPhotonCount_Scint() << G4endl;
    G4cout << "\tNumber of photons produced by cerenkov in this event : "
	   << eventInformation->GetPhotonCount_Ceren() << G4endl;
    G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
	   << eventInformation->GetAbsorptionCount() << G4endl;
    G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
	   << "this event : " << eventInformation->GetBoundaryAbsorptionCount() 
	   << G4endl;
    */
    G4cout << "\tCompleted event ID : " << event_id << G4endl;
    //G4cout << "Unacounted for photons in this event : " 
    //<< (eventInformation->GetPhotonCount_Scint() + 
    //	       eventInformation->GetPhotonCount_Ceren() - 
    //	       eventInformation->GetAbsorptionCount() -
    //	       eventInformation->GetHitCount() - 
    //	       eventInformation->GetBoundaryAbsorptionCount())
    //	   << G4endl;
  }
  //If we have set the flag to save 'special' events, save here
  if(saveThreshold&&eventInformation->GetPhotonCount() <= saveThreshold)
    G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(recorder)recorder->RecordEndOfEvent(anEvent);

  // write out event summary
  if(saveHitsFlag) 
    writeScintHitsToFile();

  //------------------------------------------------------------------------------------------

  // calculate the first hit at each MPPC
  std::ifstream timing("../timing.out");

  std::vector<G4float> event, pmtNumber, time, pmt0, pmt1, pmt2, pmt3;
  G4float a, b, c;
  
  max_pmt0 = max_pmt1 = max_pmt2 = max_pmt3 = 0;
  min_pmt0 = min_pmt1 = min_pmt2 = min_pmt3 = 1e15;
  
  //  timing.open("timing.out", std::ios::out);
  timing.open("timing.out");
  
  if(timing.is_open()) {
    
    //G4cout << "\tOpening timing.out" << G4endl;

    while (!timing.eof()) {
      if((timing >> a) && (timing >> b) && (timing >> c)){
	event.push_back(a);
	pmtNumber.push_back(b);
	time.push_back(c);
      }
      else{
	//G4cout << "\tFailed" << G4endl;
	break;
      }
    }

    //G4cout << "\tParsed in values to file" << G4endl;
    
    G4int size = pmtNumber.size();
    
    for (G4int i=0; i<size; i++) {
      if (event.at(i)==event_id) {
	if (pmtNumber.at(i)==0) {
	  pmt0.push_back(time[i]);
	}
	if (pmtNumber.at(i)==1) {
	  pmt1.push_back(time[i]);
	}
	if (pmtNumber.at(i)==2) {
	  pmt2.push_back(time[i]);
	}
	if (pmtNumber.at(i)==3) {
	  pmt3.push_back(time[i]);
	}
      }
    }
      
    //G4cout << "\tPut values in vector form" << G4endl;

    G4int pmt0Size = pmt0.size();
    G4int pmt1Size = pmt1.size();
    G4int pmt2Size = pmt2.size();
    G4int pmt3Size = pmt3.size();
    
    G4double aveTimePMT0, aveTimePMT1, aveTimePMT2, aveTimePMT3;
    aveTimePMT0 = aveTimePMT1 = aveTimePMT2 = aveTimePMT3 = 0;
    
    for (G4int i=0; i<pmt0Size; i++){
      aveTimePMT0 += pmt0.at(i);
	if (max_pmt0 < pmt0.at(i)) {
	  max_pmt0 = pmt0.at(i);
	}
	if (min_pmt0 > pmt0.at(i)) {
	  min_pmt0 = pmt0.at(i);
	}
      }
    aveTimePMT0 = aveTimePMT0/pmt0Size;
    
    for (G4int i=0; i<pmt1Size; i++){
      aveTimePMT1 += pmt1.at(i);
	if (max_pmt1 < pmt1.at(i)) {
	  max_pmt1 = pmt1.at(i);
	}
	if (min_pmt1 > pmt1.at(i)) {
	  min_pmt1 = pmt1.at(i);
	}
      }
    aveTimePMT1 = aveTimePMT1/pmt1Size;
    
    for (G4int i=0; i<pmt2Size; i++){
      aveTimePMT2 += pmt2.at(i);
	if (max_pmt2 < pmt2.at(i)) {
	  max_pmt2 = pmt2.at(i);
	}
	if (min_pmt2 > pmt2.at(i)) {
	  min_pmt2 = pmt2.at(i);
	}
      }
    aveTimePMT2 = aveTimePMT2/pmt2Size;
    
    for (G4int i=0; i<pmt3Size; i++){
      aveTimePMT3 += pmt3.at(i);
	if (max_pmt3 < pmt3.at(i)) {
	  max_pmt3 = pmt3.at(i);
	}
	if (min_pmt3 > pmt3.at(i)) {
	  min_pmt3 = pmt3.at(i);
	}
      }
    aveTimePMT3 = aveTimePMT3/pmt3Size;
    
    pmtNumber.clear();
    time.clear();
    pmt0.clear();
    pmt1.clear();
    pmt2.clear();
    pmt3.clear();

    //G4cout << "\tClosing timing.out one" << G4endl;

    timing.close();

    //G4cout << "\tClosing timing.out two" << G4endl;
  }

  //------------------------------------------------------------------------------------------

  /*
  G4cout << "MPPC[0] : " << min_pmt0 << "ns" << G4endl;
  G4cout << "MPPC[1] : " << min_pmt1 << "ns" << G4endl;
  G4cout << "MPPC[2] : " << min_pmt2 << "ns" << G4endl;
  G4cout << "MPPC[3] : " << min_pmt3 << "ns" << G4endl;
  */

  //write out the time of the first hit for each MPPC to the file pmtHit.out
  if(saveHitsFlag) 
    writePmtHitsToFile();

  if(saveHitsFlag)
    writeEnergyToFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void scintEventAction::writeScintHitsToFile(void) {

  G4String filename="avgs.out";
  //G4String filename=runAct->GetsavehitsFile();
  std::ofstream hitsfile(filename, std::ios::app);
  if(!event_id) {
    std::ofstream hitsfile(filename);
    hitsfile <<"Evt\tScintHit   scintillations   pmtHits  totalhits  initial    DeltaE" 
	     << G4endl;
    hitsfile <<"#\tns\t   #                #        photons    particle   MeV"
	     << G4endl
	     << G4endl;
  }
  
  if(hitsfile.is_open()) {
    
    hitsfile << std::setiosflags(std::ios::fixed)
	     << std::setprecision(4)
	     << std::setiosflags(std::ios::left)
	     << std::setw(6)
	     << event_id << "\t"
	     << firstScintHitTime/nanosecond << "\t   "
	     << n_hit << "\t            "
	     << pmts << "\t     "
	     << totalHits << " \t"
	     << firstParticleName << "\t   "
	     << totEnergy/MeV << " "
	     << std::setw(6)
	     << std::setiosflags(std::ios::scientific) 
	     << std::setprecision(2)
	     << std::setiosflags(std::ios::fixed) 
	     << std::setprecision(4)
	     << G4endl;
    
    hitsfile.close();
  }

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventAction::writePmtHitsToFile(void) {

  G4String filename="pmtHit.out";
  //G4String filename=runAct->GetsavepmtFile();
  std::ofstream pmtHit(filename, std::ios::app);

  if(pmtHit.is_open()) {   
    pmtHit << event_id << "    "
	   << one      << "    "
	   << two      << "    "
	   << three    << "    "
	   << four     << "    "
	   << G4endl;
    
    pmtHit.close();
  }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventAction::writeEnergyToFile(void) {

  G4String filename="energy.out";
  //G4String filename=runAct->GetsavepmtFile();
  std::ofstream energy(filename, std::ios::app);

  if(energy.is_open()) {   
    energy << event_id      << "    "
	   << totalHits     << "    "
	   << totEnergy/MeV << "    "
	   << G4endl;
    
    energy.close();
  }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintEventAction::SetSaveThreshold(G4int save){
/*Sets the save threshold for the random number seed. If the number of photons
generated in an event is lower than this, then save the seed for this event
in a file called run###evt###.rndm
*/
  saveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}




