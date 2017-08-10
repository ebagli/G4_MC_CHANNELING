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

#include "CrystalDetector.hh"
#include "CrystalDetectorHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"
#include "G4ChannelingTrackData.hh"

CrystalDetector::CrystalDetector(G4String name):
G4VSensitiveDetector(name){
    G4String HCname;
    collectionName.insert(HCname="collection");
    fHCID = -1;
    fChannelingID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetector::~CrystalDetector(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetector::Initialize(G4HCofThisEvent*HCE){
    fHitsCollection =
        new CrystalDetectorHitsCollection(SensitiveDetectorName,
                                                        collectionName[0]);
    if(fHCID<0){
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool CrystalDetector::ProcessHits(G4Step*aStep,
                                                  G4TouchableHistory*){

    if(aStep->GetTrack()->GetWeight()!=1){
        //while(!getchar());
        //G4cout << aStep->GetTrack()->GetWeight() << G4endl;
    }
    if(aStep->GetTrack()->GetTrackID()>1) return true;
    
    G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
    if((postStepPoint->GetStepStatus() == fGeomBoundary)) return true;

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    if((preStepPoint->GetStepStatus() == fGeomBoundary)) return true;
    
    if(fChannelingID == -1){
        fChannelingID = G4PhysicsModelCatalog::Register("channeling");
    }
    if(fChannelingID == -1){
        return true;
    }
    G4ChannelingTrackData* trackdata =
    (G4ChannelingTrackData*)(aStep->GetTrack()->GetAuxiliaryTrackInformation(fChannelingID));

    G4TouchableHistory* theTouchable =
        (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0);
    G4int copyNo = thePhysical->GetCopyNo();

    G4ThreeVector worldPos = preStepPoint->GetPosition();
    
    CrystalDetectorHit* aHit =
        new CrystalDetectorHit(copyNo);
    aHit->SetLayerID(copyNo);
    aHit->SetEFX(trackdata->GetEFX());
    aHit->SetEFY(trackdata->GetEFY());
    aHit->SetNud(trackdata->GetNuD());
    aHit->SetEld(trackdata->GetElD());
    G4ThreeVector pos = trackdata->GetPosCh();
    aHit->SetWorldPos(worldPos);
    aHit->SetStep(aStep->GetStepLength());
    aHit->SetKinE(aStep->GetTrack()->GetKineticEnergy());
    
    fHitsCollection->insert(aHit);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetector::EndOfEvent(G4HCofThisEvent* /*HCE*/){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
