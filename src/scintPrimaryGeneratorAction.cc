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
#include <iostream>
#include <iomanip>

#include "scintPrimaryGeneratorAction.hh"
#include "scintPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Tubs.hh"
#include "Randomize.hh"

#include "globals.hh"
#include "TG4BLStream.hh"

#include "CLHEP/GenericFunctions/Square.hh"
#include "CLHEP/GenericFunctions/Sqrt.hh"
#include "CLHEP/GenericFunctions/Cos.hh"
#include "CLHEP/GenericFunctions/Sin.hh"
#include "CLHEP/GenericFunctions/ACos.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintPrimaryGeneratorAction::scintPrimaryGeneratorAction()
{

  G4cout << "scintPrimaryGeneratorAction::scintPrimaryGeneratorAction" << G4endl;

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //create a messenger for this class
  fGunMessenger = new scintPrimaryGeneratorMessenger(this);

  //create gauss random
  fGauss = new G4RandGauss(*CLHEP::HepRandom::getTheEngine());
   
  //scintPrimaryGeneratorMessenger(this);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  //Changed the particle to an electron
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(1.*MeV);
  //particleGun->SetParticleMomentum(40.0*MeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, -10.0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintPrimaryGeneratorAction::~scintPrimaryGeneratorAction(){
    delete fParticleGun;
    delete fGunMessenger;
    delete fGauss;
}

void scintPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //G4cout << "scintPrimaryGeneratorAction::GeneratePrimaries" << G4endl;

  Genfun::Square square;
  Genfun::Sqrt sqrt;
  Genfun::Cos cos;
  Genfun::Sin sin;
  Genfun::ACos acos;

  //#define BUFSIZE 500
  //const G4int bufsize = BUFSIZE;
  //char buffer[BUFSIZE], buff2[BUFSIZE];

  G4double theWeight;

  //this function is called at the begining of event
  // 

  // Save the initial gun settings.
  G4ParticleDefinition* aDefinition = fParticleGun->GetParticleDefinition();
  G4ThreeVector aPosition = fParticleGun->GetParticlePosition();
  G4ThreeVector aDirection = fParticleGun->GetParticleMomentumDirection();
  G4double aKineticEnergy = fParticleGun->GetParticleEnergy();

  // Primary from the external file.
  if (fG4BLStream.IsOpen()) {
    fG4BLStream.Next();
    if (fG4BLStream.IsValid()) {
      G4double dx0, dy0, dz0;

      // Direction
      G4ThreeVector direction =
	G4ThreeVector(fG4BLStream.GetDirection().x(),
		      fG4BLStream.GetDirection().y(),
		      fG4BLStream.GetDirection().z());
      if (fPositionMode == "random") {
	G4ThreeVector anAxis = G4ThreeVector(0.,0.,1.).cross(direction);
	G4double aNorm = sqrt(anAxis*anAxis);
	if (aNorm == 0) {
	  anAxis = direction;
	} else {
	  anAxis = anAxis / aNorm;
	}
	G4double aCos = direction.z();
	G4double aDelta = acos(aCos);

	G4double cospolar = G4UniformRand()*(1. - cos(fDirectionSpread)) +
	  cos(fDirectionSpread);
	G4double azimuth = G4UniformRand() * 2. * pi;
	direction = G4ThreeVector(sin(acos(cospolar))*cos(azimuth),
				  sin(acos(cospolar))*sin(azimuth),
				  cospolar);
	direction.rotate(aDelta,anAxis);
      } 


      //
      // Position
      if (fPositionMode == "random") {
	do {
	  dx0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.x();
	  dy0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.y();
	} while (pow(dx0/fPositionSpread.x(),2) +
		 pow(dy0/fPositionSpread.y(),2) > 1.);
	dz0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.z();
      } else {
	dx0 = dy0 = dz0 = 0;
      }
      // shift pos[0-1] by aPositionX and aPositionY,
      // and replace pos[2] with aPositionZ.
      G4ThreeVector position =
	G4ThreeVector(fG4BLStream.GetPosition().x() + aPosition.x() + dx0,
		      fG4BLStream.GetPosition().y() + aPosition.y() + dy0,
		      aPosition.z() + dz0);


      // Set a primary.
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle
	= particleTable->FindParticle(fG4BLStream.GetPDGEncoding());
      fParticleGun->SetParticleDefinition(particle);
      fParticleGun->SetParticlePosition(position);
      fParticleGun->SetParticleMomentumDirection(direction);

      G4double p0;
      if (fEnergyMode == "random") {
	// In this part, take the fEnergySpread as fMomentumSpread.
	p0 = 2.*(G4UniformRand()-0.5) * fEnergySpread;
      } else {
	p0 = 0.;
      }
      G4double t0 = sqrt(pow(fG4BLStream.GetMomentum()+p0,2)
			 + pow(particle->GetPDGMass(),2))
	- particle->GetPDGMass();
      fParticleGun->SetParticleEnergy(t0);
      theWeight = fG4BLStream.GetWeight();

    } else {
      G4cerr << "MyPrimaryGeneratorAction: No valid G4BLStream." << G4endl;
      exit(1);
    }

  }

  /*

  else {

    G4bool done = true;
    do {
      G4double x0 = aPosition.x();
      G4double y0 = aPosition.y();
      G4double z0 = aPosition.z();
      G4double u0 = aDirection.x();
      G4double v0 = aDirection.y();
      G4double w0 = aDirection.z();
      G4double t0 = aKineticEnergy;

      //G4double dx0, dy0;

      if (fPositionMode == "target") {
	G4int idisk;
	G4ThreeVector apos;
	G4VPhysicalVolume* aVolume;
	do {
	  //       G4cout << "  ---" << fTDiskZ.size() << G4endl;
	  idisk = G4UniformRand()*fTDiskZ.size();
	  idisk = idisk<fTDiskZ.size()?idisk:(fTDiskZ.size()-1);
	  //       G4cout << idisk << G4endl;
	  apos = fTDiskConstruction->GetRandomPosition();
	  x0 = apos.x();
	  y0 = apos.y();
	  z0 = fTDiskZ[idisk] + apos.z();
	  //
	  // Get the tree list down to the current volume (fG4TouchableHistory)
	  // by using G4Navigator.
	  fG4Navigator->
	    LocateGlobalPointAndUpdateTouchable(G4ThreeVector(x0, y0, z0),
						G4ThreeVector(0.,0.,0.),
						fG4TouchableHistory,
						true);
	} while (!(fG4TouchableHistory->GetVolume(0)
		 ->GetName().contains("TStrip")));
    
      } else if (fPositionMode == "random") {
	if (fPositionSpread.x() == 0. && fPositionSpread.y() == 0.) {
	  dx0 = 0;
	  dy0 = 0;
	} else if (fPositionSpread.x() == 0.) {
	  dx0 = 0;
	  do {
	    dy0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.y();
	  } while (pow(dy0/fPositionSpread.y(),2) > 1.);
	} else if (fPositionSpread.y() == 0.) {
	  dy0 = 0;
	  do {
	    dx0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.x();
	  } while (pow(dx0/fPositionSpread.x(),2) > 1.);
	} else {
	  do {
	    dx0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.x();
	    dy0 = 2.*(G4UniformRand()-0.5) * fPositionSpread.y();
	  } while (pow(dx0/fPositionSpread.x(),2) +
		   pow(dy0/fPositionSpread.y(),2) > 1.);
	}
	x0 += dx0;
	y0 += dy0;
	z0 += 2.*(G4UniformRand()-0.5) * fPositionSpread.z();
      }


      if (fPositionMode == "target" ||
	  fPositionMode == "random" ) {
	G4ThreeVector anAxis = G4ThreeVector(0.,0.,1.).cross(aDirection);
	G4double aNorm = sqrt(anAxis*anAxis);
	if (aNorm == 0) {
	  anAxis = aDirection;
	} else {
	  anAxis = anAxis / aNorm;
	}
	G4double aCos = aDirection.z();
	G4double aDelta = acos(aCos);

	G4double cospolar = G4UniformRand()*(1. - cos(fDirectionSpread)) +
	  cos(fDirectionSpread);
	G4double azimuth = G4UniformRand() * 2. * pi;
	G4ThreeVector direction = G4ThreeVector(sin(acos(cospolar))*cos(azimuth),
						sin(acos(cospolar))*sin(azimuth),
						cospolar);
	direction.rotate(aDelta,anAxis);
	u0 = direction.x();
	v0 = direction.y();
	w0 = direction.z();
      } 


      if (fEnergyMode == "random") {
	t0 += 2.*(G4UniformRand()-0.5) * fEnergySpread;

      } else if (fEnergyMode == "dio") {
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle
	  = particleTable->FindParticle(particleName="e-");
	t0 = fWatanabeShanker->GetRandomEnergy2(G4UniformRand())/MeV -
	  particle->GetPDGMass();
      } else if (fEnergyMode == "rpc") {
	t0 = fRadiativePionCapture->Random2Energy(G4UniformRand())/MeV;
	//      cout << t0*MeV << endl;
      } else if (fEnergyMode == "pmc") {
	t0 = fProtonFromMuonCapture->Random2Energy(G4UniformRand())/MeV;
      } else if (fEnergyMode == "dispersivebeam") {
	t0 += fGauss->shoot(fBeamDispersion*(y0/cm)*MeV,fEnergySpread);

	if (t0 < fCutOffEnergy) {
	  done = false;
	} else {
	  done = true;
	}
      } else if (fEnergyMode == "cosmic") {
	G4double mom = fCosmicMuon->Random2Momentum(G4UniformRand());
	t0 = sqrt(pow(mom*GeV,2)+pow(105.195*MeV,2)) - 105.195*MeV;
	G4double theta =
	  fCosmicMuon->GetRandomAngle(mom,
				      CLHEP::HepRandom::getTheEngine());
	G4double phi = 2.*pi*G4UniformRand();
						
	u0 = sin(theta) * sin(phi);
	v0 = -cos(theta);
	w0 = sin(theta) * cos(phi);
      }



      fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(u0,v0,w0));
      fParticleGun->SetParticleEnergy(t0);
      theWeight = 1.0;

    } while (!done);

  }

  */

  // Generate the primary
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Set the weight parameter
  //for (G4int i=0; i < anEvent->GetNumberOfPrimaryVertex(); i++) {
  //G4PrimaryVertex* aVertex = anEvent->GetPrimaryVertex(i);
  //for (G4int j=0; j < aVertex->GetNumberOfParticle(); j++) {
  //  aVertex->GetPrimary(j)->SetWeight(theWeight);
  //}
  //}
//   G4cout << " weight = " << theWeight << G4endl;

  // Restore the initial gun steeings.
  fParticleGun->SetParticleDefinition(aDefinition);
  fParticleGun->SetParticlePosition(aPosition);
  fParticleGun->SetParticleMomentumDirection(aDirection);
  fParticleGun->SetParticleEnergy(aKineticEnergy);

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

void scintPrimaryGeneratorAction::FileOpen(G4String val)
{
  fG4BLStream.Open(val);
}

void scintPrimaryGeneratorAction::FileClose()
{
  fG4BLStream.Close();
}

void scintPrimaryGeneratorAction::SetEnergyMode(G4String val)
{
  fEnergyMode = val;
}

