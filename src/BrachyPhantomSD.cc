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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************
//    *                              *  
//    *    BrachyPhantomSD.cc       *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomSD.cc,v 1.14 2009-02-23 17:34:26 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//
#include "BrachyPhantomSD.hh"
#include "AnalysisManager.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

BrachyPhantomSD::BrachyPhantomSD(G4String name):G4VSensitiveDetector(name)
{
}

BrachyPhantomSD::~BrachyPhantomSD()
{
  
}

void BrachyPhantomSD::Initialize(G4HCofThisEvent*)
{
	edep = 0;
}

G4bool BrachyPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4double de = aStep->GetTotalEnergyDeposit();
	if(de <= 0)
		return true;
	edep += de;

	return true;
}

void BrachyPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
	AnalysisManager* analysis = AnalysisManager::getInstance();   

	// Fill the ntuple with position and energy deposit in the phantom
	if(edep > 0){
		analysis -> FillHistogramWithEnergy(edep/keV);
		analysis->count(edep);
		analysis->sumsquare(edep*edep);
//		G4cout << edep/keV << G4endl;
	}
}

void BrachyPhantomSD::clear()
{
} 

void BrachyPhantomSD::DrawAll()
{
}

void BrachyPhantomSD::PrintAll()
{
}



