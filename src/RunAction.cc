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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Analysis.hh"
#include "G4GenericMessenger.hh"

#include "G4MPI_USE.hh"
#ifdef G4MPI_USE
#include "G4MPImanager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():
G4UserRunAction(),
fFileName("out"){
    G4RunManager::GetRunManager()->SetPrintProgress(1);
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    
    //** Create directories **//
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(1);
    
    //** Creating ntuple **//
    analysisManager->CreateNtuple("sim", "Channeling");
    analysisManager->CreateNtupleDColumn("angXin");
    analysisManager->CreateNtupleDColumn("angYin");
    analysisManager->CreateNtupleDColumn("posXin");
    analysisManager->CreateNtupleDColumn("posYin");
    analysisManager->CreateNtupleDColumn("angXout");
    analysisManager->CreateNtupleDColumn("angYout");
    analysisManager->CreateNtupleDColumn("efx");
    analysisManager->CreateNtupleDColumn("efy");
    analysisManager->CreateNtupleDColumn("nud");
    analysisManager->CreateNtupleDColumn("eld");
    analysisManager->CreateNtupleDColumn("sx");
    analysisManager->CreateNtupleDColumn("sy");
    analysisManager->CreateNtupleDColumn("sz");
    analysisManager->CreateNtupleDColumn("energy");
    analysisManager->FinishNtuple();
    
    // -- Define messengers:
    fChangeFileName =
    new G4GenericMessenger(this, "/filename/","Change File Name" );
    
    //G4GenericMessenger::Command& changeFileName =
    fChangeFileName->DeclareProperty("set",
                                     fFileName,
                                     "Set user name." );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction(){
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*){
#ifdef G4MPI_USE
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int rank = G4MPImanager::GetManager()->GetRank();
    std::ostringstream fileName;
    fileName << fFileName << rank;
    analysisManager->OpenFile(fileName.str());
#else
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fFileName);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
