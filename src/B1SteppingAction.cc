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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "math.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
//G4Track
G4VPhysicalVolume* Volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
G4Track* track = step->GetTrack();
const G4DynamicParticle* dynpar = track->GetDynamicParticle();
G4ParticleDefinition* particle = dynpar->GetDefinition();
G4double KineticEnergy = dynpar->GetKineticEnergy();
G4ThreeVector Momentum = dynpar->GetMomentum();
G4double MomenX = dynpar->GetMomentum().x();
G4double MomenY = dynpar->GetMomentum().y();
G4double MomenZ = dynpar->GetMomentum().z();
G4double MomentumMagnitude = sqrt(pow(MomenX,2)+pow(MomenY,2)+pow(MomenZ,2));
G4StepPoint* PrePoint = step->GetPreStepPoint(); 
G4StepPoint* PostPoint = step->GetPostStepPoint();
G4ThreeVector Position = PostPoint->GetPosition();
G4double PosX = PostPoint->GetPosition().x()/cm;
G4double PosY = PostPoint->GetPosition().y()/cm;
G4double PosZ = PostPoint->GetPosition().z()/cm;
G4int TrackID = track->GetTrackID();
G4String ProcessName = PostPoint->GetProcessDefinedStep()->GetProcessName();
G4String ParticleName = particle->GetParticleName();
G4double steplen = step->GetStepLength()/cm;
G4double depositedEnergy = step->GetTotalEnergyDeposit();
G4TouchableHistory* theTouchable = (G4TouchableHistory*)(PrePoint->GetTouchable());
G4int CopyNumber = theTouchable->GetVolume()->GetCopyNo();
if (Volume->GetName() == "Pixeledminusz"){
if (ProcessName == "Transportation"){
if (ParticleName == "opticalphoton"){
std::ofstream CheckFile1;
CheckFile1.open("SONTRAC_pixeledminusz.dat", std::ios::app);
CheckFile1 << TrackID << " " << PosX << " " << PosY << " " << PosZ << " " << KineticEnergy << " " << ParticleName << " " << Volume->GetName() << " " << CopyNumber << " " << ProcessName << G4endl;
CheckFile1.close();
}
}
}
if (Volume->GetName() == "Pixeledplusz"){
if (ProcessName == "Transportation"){
if (ParticleName == "opticalphoton"){
std::ofstream CheckFile2;
CheckFile2.open("SONTRAC_pixeledplusz.dat", std::ios::app);
CheckFile2 << TrackID << " " << PosX << " " << PosY << " " << PosZ << " " << KineticEnergy << " " << ParticleName << " " << Volume->GetName() << " " << CopyNumber << " " << ProcessName << G4endl;
CheckFile2.close();
}
}
}
if (Volume->GetName() == "Pixeledplusx"){
if (ProcessName == "Transportation"){
if (ParticleName == "opticalphoton"){
std::ofstream CheckFile3;
CheckFile3.open("SONTRAC_pixeledplusx.dat", std::ios::app);
CheckFile3 << TrackID << " " << PosX << " " << PosY << " " << PosZ << " " << KineticEnergy << " " << ParticleName << " " << Volume->GetName() << " " << CopyNumber << " " << ProcessName << G4endl;
CheckFile3.close();
}
}
}
if (Volume->GetName() == "Pixeledminusx"){
if (ProcessName == "Transportation"){
if (ParticleName == "opticalphoton"){
std::ofstream CheckFile4;
CheckFile4.open("SONTRAC_pixeledminusx.dat", std::ios::app);
CheckFile4 << TrackID << " " << PosX << " " << PosY << " " << PosZ << " " << KineticEnergy << " " << ParticleName << " " << Volume->GetName() << " " << CopyNumber << " " << ProcessName << G4endl;
CheckFile4.close();
}
}
}
if (Volume->GetName() == "Xlayers" || Volume->GetName() == "Ylayers" || Volume->GetName() == "XCover1" || Volume->GetName() == "XCover2" || Volume->GetName() == "YCover1" || Volume->GetName() == "YCover2"){
if (ParticleName != "neutron"){
if (ParticleName != "opticalphoton"){  
std::ofstream CheckFile5;
CheckFile5.open("100_MeV_neutron_beam_all_secondaries.dat", std::ios::app);
CheckFile5 << TrackID << " " << PosX << " " << PosY << " " << PosZ << " " << KineticEnergy << " " << ParticleName << " " << Volume->GetName() << " " << CopyNumber << " " << ProcessName <<  " " << depositedEnergy << G4endl;
CheckFile5.close();
}
}
}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

