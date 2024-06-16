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
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
     
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", "Air");
  world_mat->SetMaterialPropertiesTable(vacuum_mt);
  G4Box* solidWorld =    
    new G4Box("WorldBox",                       //its name
      10*cm, 10*cm,10*cm);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "LogicalWorld");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

G4double Pi=3.14159265349;
G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
rotationMatrix->rotateY(90*deg); 
G4Material* scin_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
G4double A;
G4int Z;
//Define Oxygen
A = 16.0 * g/mole;
Z = 8;
G4Element* elO = new G4Element ("Oxygen", "O", Z, A); 
//Define Hydrogen 
A = 1.01 * g/mole;
Z = 1;
G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);
//Define Carbon
A = 12.01 * g/mole;
Z = 6;
G4Element* elC = new G4Element ("Carbon", "C", Z, A);
//Define Fluorine
A = 19.01 * g/mole;
Z = 9;
G4Element* elF = new G4Element ("Fluorine", "F", Z, A);
//Define PMMA (C502H8)
//NIST reference 
G4Material* PMMAcover1 = new G4Material("PMMA", 1.19*g/cm3, 3);
PMMAcover1 -> AddElement(elC, 5);
PMMAcover1 -> AddElement(elO, 2);
PMMAcover1 -> AddElement(elH, 8);
//fluorinated PMMA
G4Material* PMMAcover2 = new G4Material("fPMMA", 1.43*g/cm3, 4);
PMMAcover2 -> AddElement(elC, 5);
PMMAcover2 -> AddElement(elO, 2);
PMMAcover2 -> AddElement(elH, 7);
PMMAcover2 -> AddElement(elF, 1);
std::vector<G4double> wls_Energy = { 1.75 * eV, 2.00 * eV};
std::vector<G4double> rIndexPstyrene = { 1.59, 1.59};
std::vector<G4double> scintilFast    = { 1.0, 1.0};
std::vector<G4double> absorption1    = { 1 * m, 1 * m};
G4MaterialPropertiesTable *fMPTPStyrene = new G4MaterialPropertiesTable();
fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene);
fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1);
fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast);
fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10. * ns);
scin_mat->SetMaterialPropertiesTable(fMPTPStyrene);
scin_mat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
G4MaterialPropertiesTable *Cover1 = new G4MaterialPropertiesTable();
Cover1->AddConstProperty("SCINTILLATIONYIELD", 0. / keV);
std::vector<G4double> rcover1= { 1.490, 1.490};
Cover1->AddProperty("RINDEX", wls_Energy, rcover1);
Cover1->AddProperty("ABSLENGTH", wls_Energy, absorption1);
PMMAcover1->SetMaterialPropertiesTable(Cover1);
G4MaterialPropertiesTable *Cover2 = new G4MaterialPropertiesTable();
Cover2->AddConstProperty("SCINTILLATIONYIELD", 0. / keV);
std::vector<G4double> rcover2= { 1.420, 1.420};
Cover2->AddProperty("RINDEX", wls_Energy, rcover2);
Cover2->AddProperty("ABSLENGTH", wls_Energy, absorption1);
PMMAcover2->SetMaterialPropertiesTable(Cover2);
G4Tubs* solidEnv1 = new G4Tubs("SolidTub01", 0*cm, 0.050*cm, 2.38*cm, 0*cm, 2*Pi*rad); 
G4Tubs* solidEnv2 = new G4Tubs("SolidTub02", 0*cm, 0.052*cm, 2.38*cm, 0*cm, 2*Pi*rad);
G4Tubs* solidEnv3 = new G4Tubs("SolidTub03", 0*cm, 0.055*cm, 2.38*cm, 0*cm, 2*Pi*rad);
G4SubtractionSolid* Cover1Solid = new G4SubtractionSolid("Cover1", solidEnv2, solidEnv1);
G4SubtractionSolid* Cover2Solid = new G4SubtractionSolid("Cover2", solidEnv3, solidEnv2);
G4LogicalVolume* Xlayer =                         
new G4LogicalVolume(solidEnv1, scin_mat, "LogicSolid1");              
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector((-2.38+(0.136*j))*cm,(-2.38+(0.272*i))*cm, 0*cm), Xlayer,"Xlayers", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }   
G4LogicalVolume* XCover1 =                         
new G4LogicalVolume(Cover1Solid, PMMAcover1, "LogicSolid2");                
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector((-2.38+(0.136*j))*cm,(-2.38+(0.272*i))*cm,0*cm), XCover1,"XCover1", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }
  
G4LogicalVolume* XCover2 =                           
new G4LogicalVolume(Cover2Solid, PMMAcover2, "LogicSolid3");                     
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                        new G4PVPlacement(0,G4ThreeVector((-2.38+(0.136*j))*cm,(-2.38+(0.272*i))*cm,0*cm), XCover2,"XCover2", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }

G4LogicalVolume* Ylayer =                         
new G4LogicalVolume(solidEnv1, scin_mat, "LogicSolid4");                  
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(rotationMatrix,G4ThreeVector(0*cm,(-2.516+(0.272*i))*cm,(-2.38+(0.136*j))*cm), Ylayer,"Ylayers", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }   
G4LogicalVolume* YCover1 =                         
new G4LogicalVolume(Cover1Solid, PMMAcover1, "LogicSolid5");                      
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(rotationMatrix,G4ThreeVector(0*cm,(-2.516+(0.272*i))*cm,(-2.38+(0.136*j))*cm), YCover1,"YCover1", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }
  
G4LogicalVolume* YCover2 =                           
new G4LogicalVolume(Cover2Solid, PMMAcover2, "LogicSolid6");                     
   for (int i = 1; i < 18; i++){
   for (int j = 1; j < 35; j++){
                       new G4PVPlacement(rotationMatrix,G4ThreeVector(0*cm,(-2.516+(0.272*i))*cm,(-2.38+(0.136*j))*cm), YCover2,"YCover2", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }    
      

G4Box* solidEnv4 = new G4Box("Pixeled1", 0.068*cm, 0.068*cm, 0.00118*cm);    
G4LogicalVolume* Pixeledminusz =                         
new G4LogicalVolume(solidEnv4, world_mat, "LogicSolid7");                       
   for (int i = 1; i < 35; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector((-2.38+(0.136*j))*cm,(-2.38+(0.136*i))*cm,-2.39*cm), Pixeledminusz,"Pixeledminusz", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }   

G4LogicalVolume* Pixeledplusz =                         
new G4LogicalVolume(solidEnv4, world_mat, "LogicSolid8");                       
   for (int i = 1; i < 35; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector((-2.38+(0.136*j))*cm,(-2.38+(0.136*i))*cm,2.39*cm), Pixeledplusz,"Pixeledplusz", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }
G4Box* solidEnv5 = new G4Box("Pixeled2", 0.00118*cm, 0.068*cm,0.068*cm);     //its size
G4LogicalVolume* Pixeledminusx =                         
new G4LogicalVolume(solidEnv5, world_mat, "LogicSolid9");                        
   for (int i = 1; i < 35; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector(-2.39*cm,(-2.38+(0.136*i))*cm, (-2.38+(0.136*j))*cm), Pixeledminusx,"Pixeledminusx", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }   

G4LogicalVolume* Pixeledplusx =                         
new G4LogicalVolume(solidEnv5, world_mat, "LogicSolid10");                        
   for (int i = 1; i < 35; i++){
   for (int j = 1; j < 35; j++){
                         new G4PVPlacement(0,G4ThreeVector(2.39*cm,(-2.38+(0.136*i))*cm, (-2.38+(0.136*j))*cm), Pixeledplusx,"Pixeledplusx", logicWorld, false, j+(i-1)*34, checkOverlaps);          //                            
  }
  }            
         
  //always return the physical World
  //  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
