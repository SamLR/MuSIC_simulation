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
#ifndef scintHadronPhysics_h
#define scintHadronPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
//#include "G4ProtonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
//#include "G4HEProtonInelastic.hh"

#include "G4BinaryCascade.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

//#include "G4LElastic.hh"

class scintHadronPhysics : public G4VPhysicsConstructor
{
  public: 
    scintHadronPhysics(const G4String& name="pion");
    virtual ~scintHadronPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:

  // Hadron physics

  // Pion physics
  G4hIonisation*                  fPiPlusIonisation;
  G4MultipleScattering*           fPiPlusMultipleScattering;
  G4PionPlusInelasticProcess*     fPionPlusInelasticProcess;
  G4LEPionPlusInelastic*          fPiPlusLEPModel;
  G4HEPionPlusInelastic*          fPiPlusHEPModel;

  G4hIonisation*                  fPiMinusIonisation;
  G4MultipleScattering*           fPiMinusMultipleScattering;
  G4PionMinusInelasticProcess*    fPionMinusInelasticProcess;
  G4LEPionMinusInelastic*         fPiMinusLEPModel;
  G4HEPionMinusInelastic*         fPiMinusHEPModel;

  // Proton physics
  G4hIonisation*                  fProtonIonisation;
  G4MultipleScattering*           fProtonMultipleScattering;
  G4LEProtonInelastic*            fProtonLEPModel;
  G4HEProtonInelastic*            fProtonHEPModel;
  G4ProtonInelasticProcess*       fProtonInelasticProcess;

};


#endif

