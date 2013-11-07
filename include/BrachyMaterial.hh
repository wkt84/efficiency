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
// $Id: BrachyMaterial.hh,v 1.7 2006-06-29 15:47:42 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//
//    **********************************
//    *                                *
//    *      BrachyMaterial.hh          *
//    *                                *
//    **********************************
//
//Code developed by: Susanna Guatelli
//
// This class manages the elements and materials
//
#ifndef BrachyMaterial_H
#define BrachyMaterial_H 1

#include "globals.hh"
class G4Material;

class BrachyMaterial
{ 
public:
  BrachyMaterial();
  ~ BrachyMaterial();

public:
  void  DefineMaterials();
  G4Material* GetMat(G4String); //returns the material

private:
  G4Material* matW;
	G4Material* matGe;
	G4Material* matAl;
  G4Material* matplexiglass;
  G4Material* matPb;
  G4Material* matir192;
  G4Material* Titanium;
  G4Material* matAir;
	G4Material* matMylar;
  G4Material* matH2O;
  G4Material* soft;
  G4Material* matsteel;
  G4Material* gold;
  G4Material* matI; 
  G4Material* ceramic;
  G4Material* Vacuum; 
  G4Material* bone;
  G4Material* muscle;
  G4Material* PMMA;
  G4Material* glass;
  G4Material* platinum;
};
#endif
