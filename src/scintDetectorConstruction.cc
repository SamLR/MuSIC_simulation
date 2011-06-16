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

/*
    TODO organise const correctness; make getters/setters inline move to .hh?
*/
#include "scintDetectorConstruction.hh"
#include "scintPMTSD.hh"
#include "scintScintSD.hh"
#include "scintDetectorMessenger.hh"
#include "scintMainVolume.hh"
#include "scintWLSSlab.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"

G4bool scintDetectorConstruction::sphereOn = true;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintDetectorConstruction::scintDetectorConstruction()
  : scint_mt(NULL), MPTPStyrene(NULL), Pethylene_mt(NULL)
{
  SetDefaults();
  detectorMessenger = new scintDetectorMessenger(this);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
scintDetectorConstruction::~scintDetectorConstruction(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::DefineMaterials(){
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyeth = 1;
  G4int nC_eth = 2*polyeth;
  G4int nH_eth = 4*polyeth;

  //***Elements
  H = new G4Element("H", "H", z=1., a=1.01*g/mole);
  C = new G4Element("C", "C", z=6., a=12.01*g/mole);
  N = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);
  
  //Aluminium
  Al = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
			  density=universe_mean_density,kStateGas,0.1*kelvin,
			  1.e-19*pascal);
  //Air
  Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  //Glass
  Glass = new G4Material("Glass", density=1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);
  //Polyethylene scintillating bar
  Pethylene = new G4Material("Pethylene", density=1.032*g/cm3, 2);
  Pethylene->AddElement(H,nH_eth);
  Pethylene->AddElement(C,nC_eth);
  //Polystyrene scintillating disk
  Pstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  Pstyrene->AddElement(C, 8);
  Pstyrene->AddElement(H, 8);

  //Polyethylene for the scintillating bar
  const G4int Pethylene_NUMENTRIES = 3;
  G4double Pethylene_Energy[Pethylene_NUMENTRIES] = { 2.52*eV , 2.92*eV, 3.22*eV};

  G4double Pethylene_SCINT[Pethylene_NUMENTRIES]  = { 0.1, 1.0, 0.1};
  G4double Pethylene_RIND[Pethylene_NUMENTRIES]   = { 1.58, 1.58, 1.58};
  G4double Pethylene_ABSL[Pethylene_NUMENTRIES]   = { 210.0*cm , 210.0*cm, 210.0*cm};
  
  Pethylene_mt = new G4MaterialPropertiesTable();
  Pethylene_mt->AddProperty("FASTCOMPONENT", Pethylene_Energy, Pethylene_SCINT, Pethylene_NUMENTRIES);
  Pethylene_mt->AddProperty("SLOWCOMPONENT", Pethylene_Energy, Pethylene_SCINT, Pethylene_NUMENTRIES);
  Pethylene_mt->AddProperty("RINDEX",        Pethylene_Energy, Pethylene_RIND,  Pethylene_NUMENTRIES);
  Pethylene_mt->AddProperty("ABSLENGTH",     Pethylene_Energy, Pethylene_ABSL,  Pethylene_NUMENTRIES);
  Pethylene_mt->AddConstProperty("SCINTILLATIONYIELD",10000.0/MeV); 
  Pethylene_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  Pethylene_mt->AddConstProperty("FASTTIMECONSTANT",3.80*ns);
  Pethylene_mt->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  Pethylene_mt->AddConstProperty("YIELDRATIO",1.0);
  Pethylene->SetMaterialPropertiesTable(Pethylene_mt);

  //Polystyrene for the scintillating disk
  const G4int WLS_NUMENTRIES = 4;
  G4double WLS_Energy[WLS_NUMENTRIES] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
    
  G4double RIndexPstyrene[WLS_NUMENTRIES]={ 1.5, 1.5, 1.5, 1.5};
  G4double Absorption1[WLS_NUMENTRIES]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  G4double ScintilFast[WLS_NUMENTRIES]={0.00, 0.00, 1.00, 1.00};
  MPTPStyrene = new G4MaterialPropertiesTable();
  MPTPStyrene->AddProperty("RINDEX",WLS_Energy,RIndexPstyrene,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("ABSLENGTH",WLS_Energy,Absorption1,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("FASTCOMPONENT",WLS_Energy, ScintilFast,
			   WLS_NUMENTRIES);
  MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  MPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  Pstyrene->SetMaterialPropertiesTable(MPTPStyrene);
 
  const G4int NUMENTRIES = 3;

  G4double Glass_RIND[NUMENTRIES]={1.49,1.49,1.49};
  G4double Glass_AbsLength[NUMENTRIES]={420.*cm,420.*cm,420.*cm};
  G4MaterialPropertiesTable *Glass_mt = new G4MaterialPropertiesTable();
  Glass_mt->AddProperty("ABSLENGTH", Pethylene_Energy, Glass_AbsLength, Pethylene_NUMENTRIES);
  Glass_mt->AddProperty("RINDEX", Pethylene_Energy, Glass_RIND, Pethylene_NUMENTRIES);
  Glass->SetMaterialPropertiesTable(Glass_mt);

  G4double Vacuum_RIND[NUMENTRIES]={1.,1.,1.};  
  G4MaterialPropertiesTable *Vacuum_mt = new G4MaterialPropertiesTable();
  Vacuum_mt->AddProperty("RINDEX", Pethylene_Energy, Vacuum_RIND, Pethylene_NUMENTRIES);
  Vacuum->SetMaterialPropertiesTable(Vacuum_mt);
  Air->SetMaterialPropertiesTable(Vacuum_mt);//Give air the same rindex

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* scintDetectorConstruction::Construct(){
  DefineMaterials();
  return ConstructDetector();
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* scintDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = scint_x+d_mtl+1.*m;
  G4double expHall_y = scint_y+d_mtl+1.*m;
  G4double expHall_z = scint_z+d_mtl+1.*m;

  //Create experimental hall
  experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Vacuum,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
			      experimentalHall_log,"expHall",0,false,0);
			      // check this signature may swap name & associated logical volume

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  //Place the main volume
  if(mainVolume){
    new scintMainVolume(0,G4ThreeVector(),experimentalHall_log,false,0,this);
  }
  return experimentalHall_phys;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetDimensions(G4ThreeVector dims){
  this->scint_x=dims[0];
  this->scint_y=dims[1];
  this->scint_z=dims[2];
  updated=true;
}
 
 //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetHousingThickness(G4double d_mtl){
  this->d_mtl=d_mtl;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetNX(G4int nx){
  this->nx=nx;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetNY(G4int ny){
  this->ny=ny;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetNZ(G4int nz){
  this->nz=nz;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt){
  this->outerRadius_pmt=outerRadius_pmt;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::SetDefaults(){
  //Resets to default values
  //d_mtl=0.0635*cm;     //original values
  d_mtl=0.2*mm;
  
  //  scint_x = 17.8*cm;    //original values
  //  scint_y = 17.8*cm;
  //  scint_z = 22.6*cm;
  
  scint_x = 7.0*cm;  //take this as the radius of the scintillator
  scint_y = 7.0*cm;
  scint_z = 2.0*cm;
  
  //  nx = 2;   //original values
  //  ny = 2;
  //  nz = 3;
  
  nx = 1;
  ny = 1;
  nz = 1;

  // outerRadius_pmt = 2.3*cm; //original value
  outerRadius_pmt = 0.564*mm; //radius that gives an area equal to the 1x1mm PMT

  // sphereOn = true; //original value
  sphereOn = false;
  refl=1.0;
  
  nfibers=15;
  WLSslab=false;
  mainVolume=true;
  slab_z=2.5*mm;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/scint/detector/scintYieldFactor 1.");
  
  if(Pethylene_mt)Pethylene_mt->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  if(MPTPStyrene)MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);

  updated=true;
  
  /*  G4Tubs * SolidHadModule =
      new G4Tubs("HadModuleSolid", HadModuleRMin, HadModuleRMax, HadModuleLenght,
      HadModuleStartPhi,HadModuleDPhi);
      G4LogicalVolume * LogicalHadModule = 
      new G4LogicalVolume(SolidHadModule, FCALMaterials->Material("Copper"),
      "HadModuleLogical");
  */
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void scintDetectorConstruction::UpdateGeometry(){

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  updated=false;
}

void scintDetectorConstruction::SetMainScintYield(G4double y){
  scint_mt->AddConstProperty("SCINTILLATIONYIELD",y/MeV);  
}
 
void scintDetectorConstruction::SetWLSScintYield(G4double y){
  MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",y/MeV); 
}





