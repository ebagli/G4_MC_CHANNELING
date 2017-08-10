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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4MPI_USE.hh"

#include "G4MTRunManager.hh"
#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "UserActionInitialization.hh"
#include "DetectorConstruction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "QGSP_BIC.hh"

#include "G4EmStandardPhysics_option4_channeling.hh"
#include "G4EmStandardPhysicsSS_channeling.hh"
#include "PhysicsList.hh"
#include "G4ChannelingPhysics.hh"
#include "G4GenericBiasingPhysics.hh"


#ifdef G4MPI_USE
#include "G4MPImanager.hh"
#include "G4MPIsession.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    CLHEP::MTwistEngine randomEngine;
    G4Random::setTheEngine(&randomEngine);

#ifdef G4MPI_USE
    // --------------------------------------------------------------------
    // MPI session
    // --------------------------------------------------------------------
    // At first, G4MPImanager/G4MPIsession should be created.
    G4MPImanager* g4MPI = new G4MPImanager(argc, argv);
    
    // MPI session (G4MPIsession) instead of G4UIterminal
    // Terminal availability depends on your MPI implementation.
    G4MPIsession* session = g4MPI-> GetMPIsession();
    
    // LAM/MPI users can use G4tcsh.
    G4String prompt = "[40;01;33m";
    prompt += "G4MPI";
    prompt += "[40;31m(%s)[40;36m[%/][00;30m:";
    session-> SetPrompt(prompt);
#endif

    // Construct the default run manager
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores() - 2);

    // Activate UI-command base scorer
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);

    // Set mandatory initialization classes
    G4VModularPhysicsList* physlist= new QGSP_BIC();
    G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
    physlist->RegisterPhysics(new G4ChannelingPhysics());
    physlist->ReplacePhysics(new G4EmStandardPhysics_option4_channeling());
    //physlist->ReplacePhysics(new G4EmStandardPhysicsSS_channeling());
    biasingPhysics->PhysicsBiasAllCharged();
    physlist->RegisterPhysics(biasingPhysics);
    
    runManager->SetUserInitialization(physlist);
    runManager->SetUserInitialization(new UserActionInitialization());
    runManager->SetUserInitialization(new DetectorConstruction());
    
    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();  
    
    if(argc!=1) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    else {
        // Define visualization and UI terminal for interactive mode
#ifdef G4VIS_USE
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
#endif
        
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
        ui->SessionStart();
        delete ui;
#endif
        
#ifdef G4VIS_USE
        delete visManager;
#endif
    }
    
#ifdef G4MPI_USE
    delete g4MPI;
#endif

    // Job termination
    delete runManager;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
