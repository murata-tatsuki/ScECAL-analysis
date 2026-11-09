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
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 68030 2015-03-13 13:51:27Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
//#include "EventMessenger.hh"

#include "G4Event.hh"
#include <iomanip>
#include <HistoManager.hh>
#include "G4GenericAnalysisManager.hh"
#include "TTree.h"
#include "RunAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(HistoManager* histo)
:G4UserEventAction(),
 fEventEdep(0),fPrintModulo(10000),fDecayChain(),fHistoManager_Event(histo)
{
  	fGParticleSource  = new G4GeneralParticleSource();
  	//fHistoManager_Event = new HistoManager();
  	//fEventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  	delete fGParticleSource;
  	//delete fHistoManager_Event;
  	//delete fEventMessenger;
}

double* EventAction::CellID2Position(int cellID)
{
	int layerNo = cellID / (1e5);
	int chipNo 	= (cellID%100000) / (1e4);
	int chanNo  = (cellID%100);
	int layerOrderNo = decodeID[chipNo][chanNo];// layerOrderNo: 0-209 
	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	int _xID = layerOrderNo % rowNu; // 0-41
	int _yID = layerOrderNo / rowNu; // 0-4
	static double _position[3] = { 0 };
	double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
	double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;
	if (layerNo % 2 == 0) {
		_position[0] = -y0; 
		_position[1] = -x0;
		_position[2] = 1 + layerNo / 2 * 19.9;
	} else {
		_position[0] = -x0; 
		_position[1] = -y0;
		_position[2] = 12.2 + (layerNo - 1) / 2 * 19.9;
		//_position[2] = 12.95+(layerNo-1)/2*19.9;
	}
	return _position;
}

int EventAction::CopyNo2CellID(int copyNo) // copyNo: 0-6719
{
	int layerNo, decodeIDNow, chipNo, chanNo;
	layerNo  = copyNo / 210;
	decodeIDNow = copyNo % 210;
	bool isFound = false;
	for(int i_chip = 0 ; i_chip < 6 ; i_chip++) {
		for(int i_chan = 0 ; i_chan < 36; i_chan++) {
			if(i_chip==5 && i_chan>29) continue;
			if(decodeID[i_chip][i_chan]==decodeIDNow) {
				chipNo = i_chip;
				chanNo = i_chan;
				isFound = true;
				break;
			} else	
				continue;
		}
		if(isFound) break;
	}
	return (layerNo*1e5 + chipNo*1e4 + chanNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    fStepTag=0;
//  G4cout <<
//  		"..................6666666666666........................"
//  		<< G4endl;
  	fDecayChain = " ";
//  fHistoManager_Event->fParticleInfo.reset();
// G4cout<<"begin of event"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //  G4cout<<" >>>>>>>>>>>>>>>>> "<< fHistoManager_Event ->fParticleInfo.fPrimaryEnergy <<" <<<<<<<<<<<<<<<"<<G4endl;
//  G4cout <<
//  		"..................777777777777777........................"
//  		<< G4endl;
 	G4int evtNb = evt->GetEventID(); 
 //printing survey
 //
	if (evtNb%200 == 0) 
   	G4cout << "\n end of event " << std::setw(6) << evtNb 
          << " :" + fDecayChain << G4endl;
// G4cout<<"end of event "<<fHistoManager_Event->fParticleInfo.nTrack<<" "<<fHistoManager_Event->fParticleInfo.fTrackTime[0]<<G4endl;
 	fHistoManager_Event->fNtuple->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::SetEnergy(G4double ke)
{

	fHistoManager_Event->fParticleInfo.fPrimaryEnergy = ke;
//    G4cout<<" >>>>>>>>>>>>>>>>> "<<fHistoManager_Event->fParticleInfo.fPrimaryEnergy <<" <<<<<<<<<<<<<<<"<<G4endl;
}


void EventAction::AddParticle(G4ThreeVector direction, G4ThreeVector position)
{
	fHistoManager_Event->fParticleInfo.fPrimaryDirection.push_back(direction.x());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddCrystalEnDep(G4int copyNo,G4double edep){
	for(size_t i_copyNo=0; i_copyNo!=(fHistoManager_Event->fParticleInfo.fCrystalID.size()); ++i_copyNo){
		if(fHistoManager_Event->fParticleInfo.fCrystalID.at(i_copyNo)==copyNo){
			fHistoManager_Event->fParticleInfo.fCrystalE.at(i_copyNo) =  fHistoManager_Event->fParticleInfo.fCrystalE.at(i_copyNo) + edep;
			return;
		}
	}
	fHistoManager_Event->fParticleInfo.fCrystalID.push_back(copyNo);
	int cellID = CopyNo2CellID(copyNo);
	fHistoManager_Event->fParticleInfo.fCellID.push_back(cellID);
	double * pos = CellID2Position(cellID);
	fHistoManager_Event->fParticleInfo.fCrystalX.push_back(*(pos));
	fHistoManager_Event->fParticleInfo.fCrystalY.push_back(*(pos+1));
	fHistoManager_Event->fParticleInfo.fCrystalZ.push_back(*(pos+2));
	fHistoManager_Event->fParticleInfo.fCrystalE.push_back(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddAbsorberEnDep(G4int copyNo, G4double edep){
	for(size_t i_copyNo=0; i_copyNo!=(fHistoManager_Event->fParticleInfo.fAbsorberID.size()); ++i_copyNo){
		if(fHistoManager_Event->fParticleInfo.fAbsorberID.at(i_copyNo)==copyNo){
			fHistoManager_Event->fParticleInfo.fAbsorberE.at(i_copyNo) =  fHistoManager_Event->fParticleInfo.fAbsorberE.at(i_copyNo) + edep;
			return;
		}
	}
	fHistoManager_Event->fParticleInfo.fAbsorberID.push_back(copyNo);
	fHistoManager_Event->fParticleInfo.fAbsorberE.push_back(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddPCBEnDep(G4int copyNo,G4double edep){
	for(size_t i_copyNo=0; i_copyNo!=(fHistoManager_Event->fParticleInfo.fPCBID.size()); ++i_copyNo){
		if(fHistoManager_Event->fParticleInfo.fPCBID.at(i_copyNo)==copyNo){
			fHistoManager_Event->fParticleInfo.fPCBE.at(i_copyNo) =  fHistoManager_Event->fParticleInfo.fPCBE.at(i_copyNo) + edep;
			return;
		}
	}
	fHistoManager_Event->fParticleInfo.fPCBID.push_back(copyNo);
	fHistoManager_Event->fParticleInfo.fPCBE.push_back(edep);
}


