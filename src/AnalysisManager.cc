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
//    *    AnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: AnalysisManager.cc,v 1.22 2009-11-12 10:32:59 pandola Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//

#include <stdlib.h>
#include <fstream>
#include "AnalysisManager.hh"

#include "G4ios.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"

/*
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/ITuple.h"
*/

AnalysisManager* AnalysisManager::instance = 0;

AnalysisManager::AnalysisManager() : 
hist(0), h_ini(0)
{
  gROOT->Reset();
  // Instantiate the factories
  // The factories manage the analysis objects
//  aFact = AIDA_createAnalysisFactory();

//  AIDA::ITreeFactory *treeFact = aFact -> createTreeFactory(); 
  
  // Definition of the output file
//  G4String fileName = "brachytherapy.root";
//  theTree = treeFact -> create(fileName,"ROOT",false, true);

//  delete treeFact;
}

AnalysisManager::~AnalysisManager() 
{ 
  delete hist;
  hist = 0;

	delete h_ini;
	h_ini = 0;
}

AnalysisManager* AnalysisManager::getInstance()
{
  if (instance == 0) instance = new AnalysisManager;
  return instance;
}

void AnalysisManager::book() 
{ 
  // Instantiate the histogram and ntuple factories
  hist = new TH1F("h1", "Deposit Energy", 500, 0, 5);
	hist->GetXaxis()->SetTitle("Deposit Energy (keV)");

	h_ini = new TH1F("h2", "Initial Spectrum", 600, 0, 1200);
	h_ini->GetXaxis()->SetTitle("Energy (keV)");

	number = 0;
	hit = 0;
	square = 0;
/*  
  // Creating a 2D histogram
  // Energy deposit in the plane containing the source
  h1 = histFact -> createHistogram2D("10","Energy, pos", //histoID,histo name
				     300 ,-150.,150.,    //bins'number,xmin,xmax 
                                     300,-150.,150.);    //bins'number,ymin,ymax 
  //creating a 1D histograms
  // Histogram containing the initial energy (MeV) of the photons delivered by the radioactive core
  h2 = histFact -> createHistogram1D("20","Initial Energy", //histoID, histo name 
				     1000,0.,1.);            //bins' number, xmin, xmax
   
  // Histogram containing the energy deposit in the plane containing the source, along the axis 
  // perpendicular to the source main axis
  h3 = histFact -> createHistogram1D("30","Energy deposit  Distribution", 
				     300,-150.,150.); //bins' number, xmin, xmax

  //defining the ntuple columns' name 
  std::string columnNames = "double energy, x , y , z ";
  std::string options = "";
  
  //creating a ntuple
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);
  // check for non-zero ...
  if (ntuple) G4cout<<"The Ntuple is non-zero"<<G4endl;
*/
}
 
/*
void AnalysisManager::FillNtupleWithEnergy(G4double xx,
                                                 G4double yy, 
                                                 G4double zz,
                                                 G4double en)
{
  if (ntuple == 0) 
   {
     G4cout << "AAAAAAAGH..... The Ntuple is 0" << G4endl;
     return;
    }
 
  // Fill the ntuple
  
  G4int indexX = ntuple -> findColumn( "x" );
  G4int indexY = ntuple -> findColumn( "y" );
  G4int indexZ = ntuple -> findColumn( "z" );
  G4int indexEnergy = ntuple -> findColumn( "energy" );

  ntuple -> fill(indexEnergy, en);// method: fill ( int column, double value )
  ntuple -> fill(indexX, xx);
  ntuple -> fill(indexY, yy);
  ntuple -> fill(indexZ, zz);

  ntuple->addRow();
}
*/
void AnalysisManager::FillHistogramWithEnergy(G4double kE)
{
  // 2DHistogram: energy deposit in a voxel which center is fixed in position (x,z)  
  hist -> Fill(kE);
}

void AnalysisManager::PrimaryParticleEnergySpectrum(G4double primaryParticleEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  h_ini -> Fill(primaryParticleEnergy);
  return;
}
/*
void AnalysisManager::DoseDistribution(G4double x,G4double energy)
{
  // 1DHistogram: energy spectrum of primary particles  
  h3 -> fill(x, energy);
}
*/

void AnalysisManager::count(G4double h)
{
	number += 1;
	hit += h;
}

void AnalysisManager::sumsquare(G4double a)
{
	square += a;
}

G4double AnalysisManager::uncertainty()
{
	G4double abs, rel;
	abs = std::sqrt((square-hit*hit/number)/(number*(number-1)));
	rel = abs/(hit/number)*100;

	return rel;
}

void AnalysisManager::finish() 
{  
	TFile *f=new TFile("hist.root","RECREATE", "Geant4 ROOT Analysis");
	hist->Write();
	h_ini->Write();
	f->Close();
	delete f;
}
