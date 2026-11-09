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

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

//...oooOOO0OOOooo......................................................
PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction * mpga)
:myAction(mpga)   //  ??????????
{
	mydetDirectory = new G4UIdirectory("/mydet/");
	mydetDirectory->SetGuidance("CEPCScECAL detector control commands.");

	genCmd = new G4UIcmdWithAString("/mydet/generator", this);
	genCmd->SetGuidance("Select primary generator.");
	genCmd->SetGuidance("Available generators : PYTHIA, fGParticleSource");
	genCmd->SetParameterName("generator", true);
	genCmd->SetDefaultValue("fGParticleSource");
	genCmd->SetCandidates("PYTHIA fGParticleSource");

	hitmapCmd = new G4UIcmdWithAString("/mydet/hitmapFile", this);
  hitmapCmd->SetGuidance("Set ROOT file path for beam hitmap.");
  hitmapCmd->SetParameterName("fileName", false); // 必須引数
  hitmapCmd->SetDefaultValue("beam_hitmap.root");
}

//...oooOOO0OOOooo.......................................................
PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete genCmd;
	delete hitmapCmd;
	delete mydetDirectory;
}

//...oooOOO0OOOooo.......................................................
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command==genCmd)
	{ myAction->SetHEPEvtGenerator(newValue=="PYTHIA");}
	else if( command==hitmapCmd ) // ★【追加】コマンド受信時に LoadHitMap を呼び出す
  { 
    myAction->LoadHitMap(newValue);
  }

}

//...oooOOO0OOOooo.......................................................
G4String PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
	G4String cv;
	if (command == genCmd) {
		if(myAction->GetHEPEvtGenerator()) { 
			cv = "PYTHIA"; 
		} else {
			cv = "fGParticleSource";
		}
	}
	return cv;
}





