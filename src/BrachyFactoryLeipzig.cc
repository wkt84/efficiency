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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyFactoryLeipzig.cc  *
//    *                             *
//    *******************************
//
// $Id: BrachyFactoryLeipzig.cc,v 1.6 2006-06-29 15:48:31 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//

#include "globals.hh"
#include "BrachyFactoryLeipzig.hh"
#include"BrachyPrimaryGeneratorActionIr.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh" 
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstructionLeipzig.hh"

BrachyFactoryLeipzig:: BrachyFactoryLeipzig()
{
 leipzigSource = new  BrachyDetectorConstructionLeipzig();
 iridiumPrimaryParticle = new BrachyPrimaryGeneratorActionIr();
}

BrachyFactoryLeipzig:: ~BrachyFactoryLeipzig()
{
  delete leipzigSource;
}

void BrachyFactoryLeipzig::CreatePrimaryGeneratorAction(G4Event* anEvent)
{
 iridiumPrimaryParticle -> GeneratePrimaries(anEvent);                                 
}

void BrachyFactoryLeipzig::CreateSource(G4VPhysicalVolume* mother)
{
  leipzigSource -> ConstructLeipzig(mother);
}

void BrachyFactoryLeipzig::CleanSource()
{
  ;
}
