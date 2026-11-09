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
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4GeneralParticleSource.hh"
#include "HistoManager.hh"
#include <TFile.h>
#include <TH2D.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction* det, HistoManager* histo)
 : G4VUserPrimaryGeneratorAction(), fGParticleSource(), fDetector(det), fHistoManager_Particle(histo)
{
// const char* filename = "pythia_event.data";
//  HEPEvt = new G4HEPEvtInterface();
	fGParticleSource  = new G4GeneralParticleSource();
	useHEPEvt = false;
	a = 11;
	messenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	// delete HEPEvt;
	delete fGParticleSource;
	delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// set random particle  position
	//  
	// create vertex
	//
	fHistoManager_Particle->fParticleInfo.reset();
	if (useHEPEvt) { 
		//HEPEvt->GeneratePrimaryVertex(anEvent);
		//G4cout<<" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<G4endl;
	} else {
		if (fHitMap2D) {
	    double x_mm = 0.0;
	    double y_mm = 0.0;
	    fHitMap2D->GetRandom2(x_mm, y_mm); // ヒストグラムの形状に従って乱数取得

	    // 粒子の入射位置（X, Y, Z）を設定 (Z位置 -60mm は実験配置に合わせて変更可)
	    fGParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x_mm * mm, y_mm * mm, -60.0 * mm));
	  }

	  fGParticleSource->GeneratePrimaryVertex(anEvent);
	  //G4cout<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ "<<G4endl;
	  G4ThreeVector direction = fGParticleSource->GetParticleMomentumDirection();
	  G4ThreeVector position = fGParticleSource->GetParticlePosition();
	  G4int PDG = fGParticleSource->GetParticleDefinition()->GetPDGEncoding();
	  G4double mass = fGParticleSource->GetParticleDefinition()->GetPDGMass();
	  G4double energy = fGParticleSource->GetParticleEnergy();
	  
	  fHistoManager_Particle->fParticleInfo.fPrimaryPDG = PDG;
	  fHistoManager_Particle->fParticleInfo.fPrimaryEnergy = energy;
	  fHistoManager_Particle->fParticleInfo.fPrimaryPosX = position.x();
	  fHistoManager_Particle->fParticleInfo.fPrimaryPosY = position.y();
	  fHistoManager_Particle->fParticleInfo.fPrimaryPosZ = position.z();
	  fHistoManager_Particle->fParticleInfo.fPrimaryDirection.push_back(direction.x());
	  fHistoManager_Particle->fParticleInfo.fPrimaryDirection.push_back(direction.y());
	  fHistoManager_Particle->fParticleInfo.fPrimaryDirection.push_back(direction.z());
	  //G4cout<<angle<<" : "<<position.x()<<", "<<position.y()<<", "<<position.z()<<G4endl;
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::LoadHitMap(G4String filename)
{
  if (fHitMap2D) {
    delete fHitMap2D; // すでに読み込まれていればメモリ解放
    fHitMap2D = nullptr;
  }

  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (file && !file->IsZombie()) {
    fHitMap2D = (TH2D*)file->Get("h_hitmap_smooth"); // ヒストグラム名
    if (fHitMap2D) {
      fHitMap2D->SetDirectory(0);
      G4cout << "Successfully loaded HitMap from: " << filename << G4endl;
    }
    file->Close();
    delete file;
  } else {
    G4cerr << "Error: Cannot open " << filename << G4endl;
  }
}
