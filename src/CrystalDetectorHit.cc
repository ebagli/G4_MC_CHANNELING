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

#include "CrystalDetectorHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#ifdef G4MULTITHREADED
G4ThreadLocal G4Allocator<CrystalDetectorHit>*
CrystalDetectorHitAllocator = 0;
#else
G4Allocator<CrystalDetectorHit>
CrystalDetectorHitAllocator;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetectorHit::CrystalDetectorHit(){
    fLayerID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetectorHit::CrystalDetectorHit(G4int z){
    fLayerID = z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetectorHit::~CrystalDetectorHit(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetectorHit::CrystalDetectorHit(
                                       const CrystalDetectorHit &right): G4VHit(){
    fLayerID = right.fLayerID;
    fWorldPos = right.fWorldPos;
    fEFX = right.fEFX;
    fEFY = right.fEFY;
    fNud = right.fNud;
    fEld = right.fEld;
    fStep = right.fStep;
    fKinE = right.fKinE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const CrystalDetectorHit&
CrystalDetectorHit::operator=(
                              const CrystalDetectorHit &right){
    fLayerID = right.fLayerID;
    fWorldPos = right.fWorldPos;
    fEFX = right.fEFX;
    fEFY = right.fEFY;
    fNud = right.fNud;
    fEld = right.fEld;
    fStep = right.fStep;
    fKinE = right.fKinE;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int CrystalDetectorHit::operator==
(const CrystalDetectorHit &/*right*/) const{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetectorHit::Draw(){
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fWorldPos);
        circle.SetScreenSize(2);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,1.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::map<G4String,G4AttDef>*
CrystalDetectorHit::GetAttDefs() const{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store =
    G4AttDefStore::GetInstance("CrystalDetectorHit",isNew);
    if (isNew) {
        G4String ID("ID");
        (*store)[ID] = G4AttDef(ID,"ID","Physics","","G4int");
        
        G4String Pos("Pos");
        (*store)[Pos] =
        G4AttDef(Pos, "Position","Physics","G4BestUnit","G4ThreeVector");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AttValue>*
CrystalDetectorHit::CreateAttValues() const{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
    
    values->push_back(G4AttValue("ID",
                                 G4UIcommand::ConvertToString(fLayerID),
                                 ""));
    
    values->push_back(G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));
    
    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetectorHit::Print(){
    G4cout << fLayerID <<
    "," << fWorldPos.x() <<
    "," << fWorldPos.z() <<
    "," << fWorldPos.y() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
