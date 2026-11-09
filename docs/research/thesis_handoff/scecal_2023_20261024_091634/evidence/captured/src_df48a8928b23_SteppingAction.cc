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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
SteppingAction* SteppingAction::fgInstance =0;

SteppingAction* SteppingAction::Instance()
{
	return fgInstance;
}

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event) 
  : G4UserSteppingAction(),
    fVolume(0),
    fDetector(det), fEventAction_Step(event)                                         
{
	fgInstance = this;
	kineticEn=0;
	volume1="none";
	volume2="none";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
	fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	if (aStep->GetTrack()->GetTrackID()==1 && aStep->GetTrack()->GetCurrentStepNumber()==1) {
		preEnergy = aStep->GetTrack()->GetKineticEnergy();
	 // fEventAction_Step ->SetEnergy(preEnergy);
	}
	G4int particleID = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
	G4String name = aStep->GetTrack()->GetDefinition()->GetParticleName();
	G4int stepNo = aStep->GetTrack()->GetCurrentStepNumber();
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4int copyNo = preStepPoint->GetPhysicalVolume()->GetCopyNo();
	
	// get volume of the current step
	G4LogicalVolume* preVolume 
		= preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4String volumeName = preVolume->GetName();
	
	G4LogicalVolume* postVolume 
		= preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4String postVolumeName = postVolume->GetName();
	
	G4double eleEnergy = preStepPoint->GetKineticEnergy();

	if(name=="e-"&&volumeName=="RealCrystal") fEventAction_Step->fStepTag++;

	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep==0) return;
	// check if we are in scoring volume
	// collect energy and track length step by step
	if(volumeName=="Absorber") {  fEventAction_Step->AddAbsorberEnDep(copyNo,edep); }
	if(volumeName=="RealCrystal") {  fEventAction_Step->AddCrystalEnDep(copyNo,edep); }
	if(volumeName=="PCBCu" || volumeName=="PCBFR4") {  fEventAction_Step->AddPCBEnDep(copyNo,edep); }
	//  if(fEventAction_Step->fStepTag==1&&name=="e-"&&volumeName=="Crystal") {  fEventAction_Step->AddElectronEn(copyNo,eleEnergy); }
}
 
void SteppingAction::Reset()
{
  //fEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
