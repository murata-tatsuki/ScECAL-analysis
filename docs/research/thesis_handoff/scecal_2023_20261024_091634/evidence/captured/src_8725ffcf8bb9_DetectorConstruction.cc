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
/// \file radioactivedecay/rdecay01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "SteppingAction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // define a material
    G4Material* Vacuum =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");// vacuum  

    G4Material* Kapton =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON"); 

    G4Material* PSD =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); 

    G4Material* Aluminum =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"); 

	G4Material* Air = 
		G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");	

	G4Material* SiO2 = 
		G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");	

    G4Material* tungsten =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

    G4Material* Copper =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");

    G4Element* elW = 
        G4NistManager::Instance()->FindOrBuildElement("W");

    G4Element* elCu = 
        G4NistManager::Instance()->FindOrBuildElement("Cu");

    G4Element* elH = 
        G4NistManager::Instance()->FindOrBuildElement("H");

    G4Element* elO = 
        G4NistManager::Instance()->FindOrBuildElement("O");

    G4Element* elC = 
        G4NistManager::Instance()->FindOrBuildElement("C");

    G4double density;
	G4String name;
    G4int ncomponents;
	G4int natoms;
    G4double fractionmass;

	//Epoxy (for FR4 )
	//density = 1.2*g/cm3;
	//G4Material* Epoxy = new G4Material("Epoxy" , density, ncomponents=2);
	//Epoxy->AddElement(elH, natoms=2);
	//Epoxy->AddElement(elC, natoms=2);
	density = 1.25*g/cm3;
	G4Material* Epoxy = new G4Material("Epoxy" , density, ncomponents=3);
	Epoxy->AddElement(elC, natoms=11);
	Epoxy->AddElement(elH, natoms=12);
	Epoxy->AddElement(elO, natoms=3);

	//FR4 (Glass + Epoxy)
	density = 1.86*g/cm3;
	G4Material* FR4 = new G4Material("FR4"  , density, ncomponents=2);
	//FR4->AddMaterial(SiO2, fractionmass=0.528);
	//FR4->AddMaterial(Epoxy, fractionmass=0.472);
	FR4->AddMaterial(SiO2, fractionmass=0.5);
	FR4->AddMaterial(Epoxy, fractionmass=0.5);

	//density = 14.98*g.323; // for 25Cu:75W
	//density = 16.32*g/cm3; // for 15Cu:85W
	//density = 16.45*g/cm3; // for 15Cu:85W //yzhou
    //G4Material* CuW = new G4Material(name="CuW", density, ncomponents=2);
    //CuW->AddElement(elCu,fractionmass=15*perCent);
    //CuW->AddElement(elW,fractionmass=85*perCent);
	//density = 16.32*g/cm3; // from xrf v2 
	density = 16.17*g/cm3; // from xrf v2 
    G4Material* CuW = new G4Material(name="CuW", density, ncomponents=2);
    CuW->AddElement(elCu,fractionmass=15*perCent);
    CuW->AddElement(elW,fractionmass=85*perCent);
	
	//jiaxuan : add
	density = 14.68*g/cm3; // from xrf v1
	G4Material* CuW2 = new G4Material(name="CuW2", density, ncomponents=2);
    CuW2->AddElement(elCu,fractionmass=19.05*perCent);
    CuW2->AddElement(elW,fractionmass=80.95*perCent);

	// ESR film
 	density = 1.38*g/cm3;
	G4Material* PET = new G4Material("PET", density, ncomponents=3);
	PET->AddElement(elC,natoms=10);
	PET->AddElement(elH,natoms=8);
	PET->AddElement(elO,natoms=4);
    //******************************************************
    // Option to switch checking of volumes Overlaps
    G4bool checkOverlaps = true;
	
	//G4bool overlaps = myPhysicalVolume->CheckOverlaps(); 
	//if(overlaps)	{G4cout << "Warning: Overlaps detected in the geometry!" << G4endl;}
    
	// World : Full sphere shape
    G4double solidWorld_rmax = 100*cm;
    G4Orb* solidWorld 
			= new G4Orb("World", solidWorld_rmax);             
    G4LogicalVolume* logicWorld 
			= new G4LogicalVolume(solidWorld, Air, "World");   
    G4VPhysicalVolume* physiWorld
			= new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);       

	//******************************************************
    //Absorber & Crystal(Scintillator) & PCB
    G4int LayerNo = 30;
    G4int crystalShortSideNo = 42;
    G4int crystalLongSideNo = 5;
	G4double skin_ESR = 0.08*mm;
	G4double skin_kapton = 0.05*mm;
    
	G4double crystalL = (45+2*skin_ESR)*mm;
    G4double crystalW = (5+2*skin_ESR)*mm;
    G4double crystalH = (2+2*skin_ESR)*mm;
    G4double crystalGapL = 0.24*mm;
    G4double crystalGapW = 0.14*mm;
    
	G4double absorberXY = 230*mm;
	G4double PCBFR4XY = absorberXY;
	G4double PCBCuXY = absorberXY;
	G4double absorberZ = 3.2*mm;
    G4double PCBCuZ = 0.24*mm;
    G4double PCBFR4Z = 1.56*mm;
    
	G4double crystal_absorberGap_1 = 0.9*mm;// in a superlayer
	G4double absorber_crystalGap_2 = 0.9*mm;// between superlayers
	G4double airgap = 0.64*mm; // jiaxuan : measurment : 2.85mm scint+ESR+kapton+air
    //G4double crystal_pcbGap = 0*mm;
    G4double pcb_absorberGap = 1.2*mm;

	//******************************************************
    G4double absorberPositionXY = 0*mm;
	G4double kaptonXY = absorberXY;
	G4double kaptonPositionXY = 0*mm;
   	G4double PCBPositionXY = 0*mm;
    
	G4double absorberPositionZ1 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ/2.);
    G4double absorberPositionZ2 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + pcb_absorberGap + PCBFR4Z + PCBCuZ + crystalH + airgap + skin_kapton + crystal_absorberGap_1 + absorberZ/2.);//order of superlayer
    
	G4double crystalPositionL 	= (crystalLongSideNo - 1) / 2. * (crystalL + crystalGapL);		// now : (4/2) * (45.4mm)
    G4double crystalPositionW	= (crystalShortSideNo - 1) / 2. * (crystalW + crystalGapW); 	// now : (41/2) * (5.3mm)
    G4double crystalPositionZ1	= (skin_kapton + airgap + crystalH / 2.);
    G4double crystalPositionZ2 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + pcb_absorberGap + PCBFR4Z + PCBCuZ + crystalH/2.);
    
	G4double kaptonPositionZ1 = skin_kapton / 2. ;
	G4double kaptonPositionZ2 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + pcb_absorberGap + PCBFR4Z + PCBCuZ + crystalH + airgap + skin_kapton/2);

    G4double PCBFR4PositionZ1 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z/2.);
    G4double PCBFR4PositionZ2 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + pcb_absorberGap + PCBFR4Z/2.);
    G4double PCBCuPositionZ1 = (skin_kapton + airgap + crystalH + PCBCuZ/2.);
    G4double PCBCuPositionZ2 = (skin_kapton + airgap + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + pcb_absorberGap + PCBFR4Z + PCBCuZ/2.);
   

	G4double superlayerGapZ = (skin_kapton + crystalH + PCBCuZ + PCBFR4Z + pcb_absorberGap + absorberZ + airgap )*2 + crystal_absorberGap_1 + absorber_crystalGap_2;

    //******************************************************
    //******************************************************
 	// Kapton Film
  	G4Box* solidKapton = new G4Box("Kapton",0.5*kaptonXY, 0.5*kaptonXY, 0.5*skin_kapton);
	G4LogicalVolume* logicKapton = new G4LogicalVolume(solidKapton,Kapton,"Kapton");
	for(G4int i_Layer=0; i_Layer < LayerNo; i_Layer++)
	{
		G4double XNow = kaptonPositionXY;
		G4double YNow = kaptonPositionXY;
		G4double ZNow = i_Layer%2==0 ? (kaptonPositionZ1+i_Layer/2*superlayerGapZ) : (kaptonPositionZ2+(i_Layer-1)/2*superlayerGapZ);
		G4VPhysicalVolume* physiKapton 
			= new G4PVPlacement(0, G4ThreeVector(XNow,YNow, ZNow), logicKapton, "Kapton", logicWorld, false, i_Layer, checkOverlaps);

	}


  	//Scintillator && ESR Film
    G4Box* solidCrystal
			= new G4Box("Crystal", 0.5*crystalL, 0.5*crystalW, 0.5*crystalH);    
	G4Box* solidRealCrystal 
			= new G4Box("RealCrystal", 0.5*crystalL-skin_ESR, 0.5*crystalW-skin_ESR, 0.5*crystalH-skin_ESR);
	G4SubtractionSolid* solidfilm
			= new G4SubtractionSolid("ESR",solidCrystal,solidRealCrystal);

    G4LogicalVolume* logicRealCrystal
			= new G4LogicalVolume(solidRealCrystal, PSD, "RealCrystal"); 
	G4LogicalVolume* logicFilm
			= new G4LogicalVolume(solidfilm, PET, "Film");

    for(G4int i_Layer=0; i_Layer!=LayerNo; ++i_Layer)
	{
        G4RotationMatrix rm;
        rm.rotateZ(90*deg);
        
		if(i_Layer%2==0)
		{
            //for(G4int i_Hor=0; i_Hor!=1; ++i_Hor)
            for(G4int i_Hor=0; i_Hor!=crystalLongSideNo; ++i_Hor)
			{
				//for(G4int i_Ver=0; i_Ver!=20; ++i_Ver)
                for(G4int i_Ver=0; i_Ver!=crystalShortSideNo; ++i_Ver)
				{
					// Real Crystal && Film
					G4double XNow = - (i_Hor * (crystalL + crystalGapL) - crystalPositionL) ;
					G4double YNow = - (i_Ver * (crystalW + crystalGapW) - crystalPositionW) ;
					G4double ZNow = (i_Layer/2*superlayerGapZ + crystalPositionZ1);
					G4int copyNo = i_Layer * crystalShortSideNo * crystalLongSideNo + i_Hor * crystalShortSideNo + i_Ver; 
					
					G4VPhysicalVolume* physicalCrystal
							= new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZNow), logicRealCrystal, "RealCrystal", logicWorld, false, copyNo, checkOverlaps);
					G4VPhysicalVolume* physicalFilm
							= new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZNow), logicFilm, "Film", logicWorld, false,	copyNo, checkOverlaps);    
				}
            }
        }
		if(i_Layer%2==1)
		{
            //for(G4int i_Ver=0; i_Ver!=2; ++i_Ver)
            for(G4int i_Ver=0; i_Ver!=crystalLongSideNo; ++i_Ver)
			{
                //for(G4int i_Hor=0; i_Hor!=30; ++i_Hor)
                for(G4int i_Hor=0; i_Hor!=crystalShortSideNo; ++i_Hor)
				{
					// Real Crystal && Film
					G4double XNow = - (i_Hor * (crystalW + crystalGapW) - crystalPositionW) ;
					G4double YNow = - (i_Ver * (crystalL + crystalGapL) - crystalPositionL) ;
					G4double ZNow = (i_Layer/2*superlayerGapZ + crystalPositionZ2);
					G4int copyNo = i_Layer * crystalShortSideNo * crystalLongSideNo + i_Ver * crystalShortSideNo + i_Hor; 

					G4VPhysicalVolume* physicalCrystal
							= new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(XNow,YNow,ZNow)), logicRealCrystal, "RealCrystal", logicWorld, false, copyNo, checkOverlaps);		
					G4VPhysicalVolume* physicalFilm
							= new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(XNow,YNow,ZNow)), logicFilm,	"Film",	logicWorld, false, copyNo, checkOverlaps);	

                }
            }
        }
    }

   //******************************************************
    //Absorber
    G4Box* solidAbsorber 
			= new G4Box("Absorber", 0.5*absorberXY, 0.5*absorberXY, 0.5*absorberZ); 
	G4LogicalVolume* logicAbsorber 
			= new G4LogicalVolume(solidAbsorber, CuW, "Absorber"); 
	// jiaxuan : add
	//G4LogicalVolume* logicAbsorber2 
	//		= new G4LogicalVolume(solidAbsorber, CuW2, "Absorber2");                                          // its name 

    for(G4int i_Layer=0; i_Layer!=LayerNo; ++i_Layer)
	{
		G4double XNow = absorberPositionXY;
		G4double YNow = absorberPositionXY;
		G4double ZNow = i_Layer%2==0 ? (absorberPositionZ1+i_Layer/2*superlayerGapZ) : (absorberPositionZ2+(i_Layer-1)/2*superlayerGapZ) ;
		//jiaxuan : add
		//if(i_Layer==14 || i_Layer==16 || i_Layer==18)
        new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZNow), logicAbsorber, "Absorber", logicWorld,  false, i_Layer, checkOverlaps);  
        //new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZNow), logicAbsorber2, "Absorber", logicWorld,  false, i_Layer, checkOverlaps);  
	}

    //******************************************************
    //PCB
    G4Box* solidPCBCu 
			= new G4Box("PCBCu", 0.5*PCBCuXY, 0.5*PCBCuXY, 0.5*PCBCuZ);   
	G4Box* solidPCBFR4 
			= new G4Box("PCBFR4", 0.5*PCBFR4XY, 0.5*PCBFR4XY, 0.5*PCBFR4Z);        
    G4LogicalVolume* logicPCBCu 
			= new G4LogicalVolume(solidPCBCu, Copper, "PCBCu");                                          // its name 
    G4LogicalVolume* logicPCBFR4 
			= new G4LogicalVolume(solidPCBFR4, FR4, "PCBFR4");                                          // its name 

    //  G4VPhysicalVolume* physiPCB;
    for(G4int i_Layer=0; i_Layer!=LayerNo; ++i_Layer)
	{
		G4double XNow = PCBPositionXY;
		G4double YNow = PCBPositionXY;
		G4double ZCuNow = i_Layer%2==0 ? (PCBCuPositionZ1+i_Layer/2*superlayerGapZ) : (PCBCuPositionZ2+(i_Layer-1)/2*superlayerGapZ);
		G4double ZFR4Now = i_Layer%2==0 ? (PCBFR4PositionZ1+i_Layer/2*superlayerGapZ) : (PCBFR4PositionZ2+(i_Layer-1)/2*superlayerGapZ);
		new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZCuNow), logicPCBCu, "PCBCu", logicWorld, false, i_Layer, checkOverlaps); 	
		new G4PVPlacement(0, G4ThreeVector(XNow,YNow,ZFR4Now), logicPCBFR4, "PCBFR4", logicWorld, false, i_Layer, checkOverlaps);            
    }

// **********************************************************************
// Additional material1 : Shell ECAL Al alloy 1mm
	G4Box *solidAdd1 
			= new G4Box("Shell",115*mm,115*mm,0.5*mm);
	G4LogicalVolume *logicAdd1
			= new G4LogicalVolume(solidAdd1,Aluminum,"Shell");
	new G4PVPlacement(0,G4ThreeVector(0,0,-2*mm),logicAdd1,"Shell",logicWorld,false,0,checkOverlaps);

// Additional material2 : Scintillator 41mm
	G4Box *solidAdd2 
			= new G4Box("Scint",115*mm,115*mm,20.5*mm);
	G4LogicalVolume *logicAdd2
			= new G4LogicalVolume(solidAdd2,PSD,"Scint");
	new G4PVPlacement(0,G4ThreeVector(0,0,-30*mm),logicAdd2,"Scint",logicWorld,false,0,checkOverlaps);

// **********************************************************************
// **********************************************************************
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.9,0.0,0.0));
    //G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    visAttributes ->SetVisibility(false);
    logicWorld ->SetVisAttributes(visAttributes);
    //logicAbsorber ->SetVisAttributes(visAttributes);
    //logicPCB ->SetVisAttributes(visAttributes);
	//logicFilm ->SetVisAttributes(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour(0.9,0.0,0.0));//red
    logicAbsorber ->SetVisAttributes(visAttributes);
    //logicAbsorber2 ->SetVisAttributes(visAttributes);
    logicAdd1 ->SetVisAttributes(visAttributes);
    //(0.9,0.9,0.0)yellow

    visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0));//green
    logicPCBCu ->SetVisAttributes(visAttributes);
    logicPCBFR4 ->SetVisAttributes(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));//white
    logicRealCrystal ->SetVisAttributes(visAttributes);
    logicAdd2 ->SetVisAttributes(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour::Yellow());
    logicFilm ->SetVisAttributes(visAttributes);


    // Set scoring volume to stepping action(where we will account energy deposit)
    SteppingAction* steppingAction = SteppingAction::Instance();
    steppingAction->SetVolume(logicRealCrystal);

    //always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
