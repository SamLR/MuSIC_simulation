//
//  scintBarConstruction.h
//  MuSIC_simulation
//
//  Created by Sam Cook on 19/06/2011.
//  Copyright 2011 UCL. All rights reserved.
//

#ifndef scintBarConstruction_H
#def scintBarConstruction_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "scintDetectorConstruction.hh"

class scintBarConstruction : public G4PVPlacement
{
public:
    scintBarConstruction(G4RotationMatrix *pRot,
                         const G4ThreeVector &tlate,
                         G4LogicalVolume *pMotherLogical,
                         G4bool pMany,
                         G4int pCopyNo,
                         scintDetectorConstruction* c);
    const void init();

private:
    
    
};