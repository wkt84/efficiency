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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.cc     *
//    *                                      *
//    ****************************************
//
#include "G4CSGSolid.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4TransportationManager.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "BrachyMaterial.hh"
#include "BrachyFactoryLeipzig.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"
#include "BrachyPhantomROGeometry.hh"
#include "BrachyPhantomSD.hh"
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName)
: detectorChoice(0), phantomSD(0), phantomROGeometry(0), factory(0),
  World(0), WorldLog(0), WorldPhys(0),
  phantomAbsorberMaterial(0), Detector(0), DetectorLog(0), DetectorPhys(0),GRDLog(0),GRDPhys(0)
{
  // Define half size of the phantom along the x, y, z axis
  phantomSizeX = 15.*cm ;
  phantomSizeY = 15.*cm;
  phantomSizeZ = 15.*cm;

  // Define the number of voxels the phantom is subdivided along the
  // three axis 
  // Each voxel is 1 mm wide
  numberOfVoxelsAlongX = 300;
  numberOfVoxelsAlongY = 300;
  numberOfVoxelsAlongZ = 300;

  ComputeDimVoxel();
  
  // Define the sizes of the World volume contaning the phantom
  worldSizeX = 1.0*m;
  worldSizeY = 1.0*m;
  worldSizeZ = 1.0*m;

  sensitiveDetectorName = SDName;

  // Define the messenger of the Detector component
  // It is possible to modify geometrical parameters through UI
  detectorMessenger = new BrachyDetectorMessenger(this);

  // Define the Iridium source as default source modelled in the geometry
  factory = new BrachyFactoryIr();

  // BrachyMaterial defined the all the materials necessary
  // for the experimental set-up 
  pMaterial = new BrachyMaterial();
}

BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 
  delete pMaterial;
  delete factory;
  delete detectorMessenger;
  if (phantomROGeometry) delete phantomROGeometry;
}

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{
  pMaterial -> DefineMaterials();
 
  // Model the phantom (water box)  
  ConstructPhantom();

  // Model the source in the phantom
  factory -> CreateSource(WorldPhys); //Build the source inside the phantom
 
  // Define the sensitive volume: phantom
  ConstructSensitiveDetector();

  return WorldPhys;
}

void BrachyDetectorConstruction::SwitchBrachytherapicSeed()
{
  // Change the source in the water phantom
  factory -> CleanSource();
  delete factory;

  switch(detectorChoice)
  { 
    case 1:
      factory = new BrachyFactoryI();
      break;
    case 2:
      factory = new BrachyFactoryLeipzig();
      break;
    case 3:
      factory = new BrachyFactoryIr();
      break;   
    default:
      factory = new BrachyFactoryIr();
      break;
  }

  factory -> CreateSource(WorldPhys);

  // Notify run manager that the new geometry has been built
  G4RunManager::GetRunManager() -> DefineWorldVolume( WorldPhys );
}

void BrachyDetectorConstruction::SelectBrachytherapicSeed(G4String val)
{
  if(val == "Iodium") 
  {
   detectorChoice = 1;
  }
  else
  {
    if(val=="Leipzig")
    {
      detectorChoice = 2;
    }
    else
    {
      if(val=="Iridium")
      {
        detectorChoice = 3;
      }
    }
  }

  G4cout << "Now the source is " << val << G4endl;
}

void BrachyDetectorConstruction::SetBuildupDepth(G4double val)
{
}

void BrachyDetectorConstruction::ConstructPhantom()
{
  // Model the water phantom 
  
  // Define the light blue color
  G4Colour  lblue   (0.0, 0.0, .75);

//  G4Material* Water = pMaterial -> GetMat("Water");
	G4Material* air = pMaterial -> GetMat("Air");
	G4Material* water = pMaterial -> GetMat("Water");
	G4Material* Al = pMaterial -> GetMat("Alminium");
	G4Material* mylar = pMaterial -> GetMat("Mylar");

  ComputeDimVoxel();

  // World volume
  World = new G4Box("World",worldSizeX,worldSizeY,worldSizeZ);
  WorldLog = new G4LogicalVolume(World,air,"WorldLog",0,0,0);
  WorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                "WorldPhys",WorldLog,0,false,0);

	// Alminium cap
	G4Tubs* Outer = new G4Tubs("Outer", 38.3*mm, 39.3*mm, 97*mm/2, 0.*deg, 360.*deg);
	G4LogicalVolume *OuterLog = new G4LogicalVolume(Outer, Al, "OuterLog", 0, 0, 0);
	G4VPhysicalVolume *OuterPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -49.5*mm), OuterLog, "OuterPhys", WorldLog, false, 0);

	G4Tubs* Upper = new G4Tubs("Upper", 0.*mm, 39.3*mm, 1.*mm/2, 0.*deg, 360.*deg);
	G4LogicalVolume *UpperLog = new G4LogicalVolume(Upper, Al, "UpperLog", 0, 0, 0);
	G4VPhysicalVolume *UpperPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*mm), UpperLog, "UpperPhys", WorldLog, false, 0);

	G4Tubs* Inner = new G4Tubs("Inner", 32.5*mm, 33.3*mm, 47.*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *InnerLog = new G4LogicalVolume(Inner, Al, "InnerLog", 0, 0, 0);
	G4VPhysicalVolume *InnerPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -47*mm-4.*mm), InnerLog, "InnerPhys", WorldLog, false, 0);

	G4Tubs* AlCap = new G4Tubs("AlCap", 0.*mm, 33.3*mm, 0.03*mm/2, 0.*deg, 360.*deg);
	G4LogicalVolume *AlCapLog = new G4LogicalVolume(AlCap, Al, "AlCapLog", 0, 0, 0);
	G4VPhysicalVolume *AlCapPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -4.*mm+0.045*mm),AlCapLog, "AlCapPhys", WorldLog, false, 0);

	G4Tubs* MylarCap = new G4Tubs("MylarCap", 0.*mm, 33.3*mm, 0.03*mm/2, 0.*deg, 360.*deg);
	G4LogicalVolume *MylarCapLog = new G4LogicalVolume(MylarCap, mylar, "MylarCapLog", 0, 0, 0);
	G4VPhysicalVolume *MylarCapPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -4.*mm+0.015*mm),MylarCapLog, "MylarCapPhys", WorldLog, false, 0);

  WorldLog -> SetVisAttributes (G4VisAttributes::Invisible);

  // Visualization attributes of the phantom
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(lblue);
  G4VisAttributes* MylarVisAtt = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2));
  simpleBoxVisAtt -> SetVisibility(true);
	MylarVisAtt -> SetVisibility(true);
  OuterLog -> SetVisAttributes(simpleBoxVisAtt);
  UpperLog -> SetVisAttributes(simpleBoxVisAtt);
  InnerLog -> SetVisAttributes(simpleBoxVisAtt);
  AlCapLog -> SetVisAttributes(simpleBoxVisAtt);
  MylarCapLog -> SetVisAttributes(MylarVisAtt);
}

void  BrachyDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4Colour  lgreen  (0.0, .75, 0.0);
	G4Colour  dgreen  (0.0, .90, 0.0);
  G4Material* air = pMaterial -> GetMat("Air") ;
	G4Material* germa = pMaterial -> GetMat("Germanium");

	// r of top op Germanium
	G4double r = 4.;

	// Dead Layer of Ge-Li
	G4double dead = 0.7;


	// Cryo hole
	G4Tubs *center = new G4Tubs("center", 0, 4.4*mm, 18.6*mm/2, 0.*deg, 360.*deg);
	G4Sphere *top = new G4Sphere("top", 0, 4.4*mm, 0.*deg, 360.*deg, 0.*deg, 90.*deg);

	G4UnionSolid *uni = new G4UnionSolid("center+top", center, top, 0, G4ThreeVector(0,0,18.6*mm/2));


	// Outer Crystal
	G4Tubs *tube = new G4Tubs("tube", 0, 59.*mm/2, (36.7-r)*mm/2, 0.*deg, 360.*deg);

	G4Torus *torus = new G4Torus("torus", 0., r*mm, (59./2.-r)*mm, 0.*deg, 360.*deg);
	G4Tubs *toptube = new G4Tubs("toptube", 0., (59./2.-r)*mm, r*mm, 0.*deg, 360.*deg);
	
	G4UnionSolid *uni2 = new G4UnionSolid("torus+toptube", toptube, torus, 0, G4ThreeVector(0,0,0));

	G4UnionSolid *uni3 = new G4UnionSolid("tube+uni2", tube, uni2, 0, G4ThreeVector(0,0,(36.7-r)*mm/2));
	G4SubtractionSolid *Crystal_Outer = new G4SubtractionSolid("Crystal_Outer", uni3, uni, 0, G4ThreeVector(0,0,-(36.7-r)*mm/2.+18.6*mm/2.));


	// Inner Crystal
	G4Tubs *intube = new G4Tubs("intube", 0., 59.*mm/2-dead*mm, (36.7-r)*mm/2, 0.*deg, 360.*deg);
	G4Torus *intorus = new G4Torus("intorus", 0., (r-dead)*mm, (59./2.-r)*mm, 0.*deg, 360.*deg);
	G4Tubs *intoptube = new G4Tubs("intoptube", 0., (59./2.-r)*mm, (r-dead)*mm, 0.*deg, 360.*deg);

	G4UnionSolid *uni_in1 = new G4UnionSolid("intorus+intoptube", intoptube, intorus, 0, G4ThreeVector(0,0,0));
	G4UnionSolid *uni_in2 = new G4UnionSolid("intube+uni_in1", intube, uni_in1, 0, G4ThreeVector(0,0,(36.7-r)*mm/2.));
	G4SubtractionSolid *Crystal = new G4SubtractionSolid("Crystal", uni_in2, uni, 0, G4ThreeVector(0,0,-(36.7-r)*mm/2.+18.6*mm/2.));


	// Dead Layer
	G4SubtractionSolid *DeadLayer = new G4SubtractionSolid("DeadLayer", Crystal_Outer, Crystal, 0, G4ThreeVector(0,0,0));


  GRDLog = new G4LogicalVolume(Crystal,germa,"GRDLog",0,0,0);
  GRDPhys = new G4PVPlacement(0, G4ThreeVector(0,0,-(36.7-r)*mm/2-r*mm-4.*mm), GRDLog, "GRDPhys",WorldLog,false,0);

	G4LogicalVolume *DeadLog = new G4LogicalVolume(DeadLayer, germa, "DeadLog", 0, 0, 0);
	G4VPhysicalVolume *DeadPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, -(36.7-r)*mm/2-r*mm-4.*mm), DeadLog, "DeadPhys", WorldLog, false, 0);


  // Visualization attributes of the phantom
  G4VisAttributes* simpleDetVisAtt = new G4VisAttributes(lgreen);
  simpleDetVisAtt -> SetVisibility(true);
  simpleDetVisAtt -> SetForceWireframe(true);
  GRDLog -> SetVisAttributes(simpleDetVisAtt);

	G4VisAttributes* DeadVisAtt = new G4VisAttributes(dgreen);
	DeadVisAtt -> SetVisibility(true);
	DeadVisAtt -> SetForceWireframe(true);
	DeadLog -> SetVisAttributes(DeadVisAtt);

  phantomSD = new BrachyPhantomSD(sensitiveDetectorName);
  pSDManager->AddNewDetector(phantomSD);
  GRDLog->SetSensitiveDetector(phantomSD);

}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the phantom is a water box whose size is: " << G4endl
         << phantomSizeX *2./cm
         << " cm * "
         << phantomSizeY *2./cm
         << " cm * "
         << phantomSizeZ *2./cm
         << " cm" << G4endl
         << "number of Voxel: "
         << numberOfVoxelsAlongX <<G4endl
         << "Voxel size: "
         << dimVoxel * 2/mm
         << "mm" << G4endl 
         << "The phantom is made of "
         << phantomAbsorberMaterial -> GetName() <<G4endl
         << "the source is at the center of the phantom" << G4endl
         << "-------------------------------------------------------------------------"
         << G4endl;
}

void BrachyDetectorConstruction::SetPhantomMaterial(G4String materialChoice)
{
  // It is possible to change the material of the phantom
  // interactively

  // Search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
  {
    phantomAbsorberMaterial = pttoMaterial;
    PrintDetectorParameters();
  } 
  else
    G4cout << "WARNING: material '" << materialChoice
           << "' not available!" << G4endl;            
}
