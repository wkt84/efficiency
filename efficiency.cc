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
//    *******************************
//    *                             *
//    *    Brachy.cc                *
//    *                             *
//    *******************************
//
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "BrachyFactoryIr.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "BrachyDetectorConstruction.hh"
#include "BrachyPhysicsList.hh"
#include "BrachyPhantomSD.hh"
#include "BrachyPrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include "BrachyRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"

#ifdef G4ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "G4ScoringManager.hh"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc ,char ** argv)

{
  G4ScoringManager::GetScoringManager(); // instanciate the interactive scoring manager

  G4RunManager* pRunManager = new G4RunManager;

  G4String sensitiveDetectorName = "Phantom";

  // Initialize the physics component
  pRunManager -> SetUserInitialization(new BrachyPhysicsList);

  // Initialize the detector component
  BrachyDetectorConstruction  *pDetectorConstruction = new  BrachyDetectorConstruction(sensitiveDetectorName);
  pRunManager -> SetUserInitialization(pDetectorConstruction);

  // Initialize the primary particles
  BrachyPrimaryGeneratorAction* primary = new BrachyPrimaryGeneratorAction();
  pRunManager -> SetUserAction(primary);

  BrachyRunAction *pRunAction = new BrachyRunAction();
  pRunManager -> SetUserAction(pRunAction);

  //// Initialize the Visualization component 
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // Output environment variables concerning the analysis:
#ifdef G4ANALYSIS_USE
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " Using AIDA 3.2.1 analysis " << G4endl;
# else
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " G4ANALYSIS_USE environment variable not set, NO ANALYSIS " 
	 << G4endl;
#endif
  
  //Initialize G4 kernel
  pRunManager -> Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  if (argc == 1)   // Define UI session for interactive mode.
    { 
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      G4cout << " UI session starts ..." << G4endl;
      UImanager -> ApplyCommand("/control/execute VisualisationMacro.mac");    
      ui -> SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }  

  // Close the output file containing histograms and n-tuple with the 
  // results of the simulation.
  // Finish the analysis.
#ifdef G4ANALYSIS_USE
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
  analysis -> finish();
#endif
  
  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete pRunManager;

  return 0;
}
