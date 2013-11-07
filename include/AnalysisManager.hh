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
// $Id: BrachyAnalysisManager.hh,v 1.14 2008-06-05 13:45:39 cirrone Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $

// Code review: MG Pia, 14/05/2207
// Contact: Geant4-INFN Genova group, MariaGrazia.Pia@ge.infn.it

// The class BrachyAnalysisManager creates and manages histograms and ntuples
//

#ifndef G4ANALYSISMANAGER_HH
#define G4ANALYSISMANAGER_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

class TH1F;

class AnalysisManager
{

public:
  ~AnalysisManager();
  static AnalysisManager* getInstance();

  void book();
	void PrimaryParticleEnergySpectrum(G4double);
  void FillHistogramWithEnergy(G4double);
  void finish();
	void count(G4double);
	void sumsquare(G4double);
	G4double uncertainty();

private:
  AnalysisManager();
  static AnalysisManager* instance;

	TH1F* hist; 
	TH1F* h_ini;

	G4int number;
	G4double hit;
	G4double square;

public:
	G4int printnumber() const {return number;}
	G4double printsquare() const {return square;}
	G4double printhit() const {return hit;}

#ifdef G4ANALYSIS_USE
  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory     *tupFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram2D *h1;
  AIDA::IHistogram1D *h2;
  AIDA::IHistogram1D *h3;
  AIDA::ITuple *ntuple;
#endif

};

#endif



