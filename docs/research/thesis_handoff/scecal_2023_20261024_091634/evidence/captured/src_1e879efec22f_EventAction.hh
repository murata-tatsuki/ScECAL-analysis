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
/// \file DBDecay/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 68017 2013-03-13 13:29:53Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "DetectorConstruction.hh"
class EventMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
    public:
        EventAction(HistoManager*);
        ~EventAction();

    public:
        virtual void BeginOfEventAction(const G4Event*);
        virtual void   EndOfEventAction(const G4Event*);
        G4int fStepTag;

        void SetPrintModulo(G4int val)   {fPrintModulo = val;};
        void AddDecayChain(G4String val) {fDecayChain += val;};

		double* CellID2Position(int cellID);
		int CopyNo2CellID(int copyNo);
        void SetEnergy(G4double ke);
        void AddParticle(G4ThreeVector direction, G4ThreeVector position);
        void AddCrystalEnDep(G4int copyNo,G4double edep);
		void AddAbsorberEnDep(G4int copyNo, G4double edep);
        void AddPCBEnDep(G4int copyNo,G4double edep);
	public:
		static constexpr const int decodeID[6][36] =
		{  // decodeID[chipNo][chanNo] = copyNo (in a layer)
			0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13,
			55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,

			24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 
			33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,

			149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 
			108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,

			150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 
			159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
			
			191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181,
			139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
			
			137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 
			93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, -1, -1, -1, -1, -1, -1
		}; // chip0_chan1 -> 42 
    private:
        G4double     fEventEdep;
        G4int        fPrintModulo;
        G4String     fDecayChain;                   
        
		HistoManager* fHistoManager_Event;
        G4GeneralParticleSource * fGParticleSource;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


