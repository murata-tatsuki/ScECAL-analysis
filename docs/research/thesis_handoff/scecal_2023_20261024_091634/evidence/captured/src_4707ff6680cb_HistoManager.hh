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
// $Id: HistoManager.hh 68017 2013-03-13 13:29:53Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef HistoManager_h
#define HistoManager_h 

#include "globals.hh"
#include <vector>
//#include "g4root.hh"
#include "G4GenericAnalysisManager.hh"
#include <G4ThreeVector.hh>
#include <TTree.h>
#include <TFile.h>

//const G4int kMAXTrack=5000;//should be checked!!!
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ParticleInfo
{
    public:
        G4int 	 fPrimaryPDG;
        G4double fPrimaryEnergy;
        G4double fPrimaryPosX;
        G4double fPrimaryPosY;
        G4double fPrimaryPosZ;
        std::vector<G4double> 	fPrimaryDirection;
        std::vector<G4int>	  	fCrystalID;
        std::vector<G4int>	  	fCellID;
        std::vector<G4double> 	fCrystalE;
        std::vector<G4double> 	fCrystalX;
        std::vector<G4double> 	fCrystalY;
        std::vector<G4double> 	fCrystalZ;
        std::vector<G4int> 	  	fAbsorberID;
        std::vector<G4double> 	fAbsorberE;
        std::vector<G4int> 		fPCBID;
        std::vector<G4double> 	fPCBE;
        //std::vector<G4int> 		fFilmID;
        //std::vector<G4double> 	fFilmEnDep;
	public:
		void reset();
};

class HistoManager
{
    public:
        HistoManager(char* foutname);
        ~HistoManager();
        void save();
        void book();
        ParticleInfo fParticleInfo;
    private:
        G4String fFileName;
        G4String fOutFileName;
    public:  
        TFile* fRootFile;
        TTree* fNtuple;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

