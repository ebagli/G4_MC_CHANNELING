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

#include "G4ChannelingMaterialData.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ChannelingECHARM.hh"
#include "G4LogicalCrystalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Box.hh"

G4ChannelingMaterialData::G4ChannelingMaterialData(const G4String& name):
G4VMaterialExtension(name),
bIsBent(false){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMaterialData::SetFilename(const G4String& fileName){
    G4String filePot = fileName + "_pot.txt";
    G4String fileEFX = fileName + "_efx.txt";
    G4String fileEFY = fileName + "_efy.txt";
    G4String fileAtD = fileName + "_atd.txt";
    G4String fileElD = fileName + "_eld.txt";

    fPotential = new G4ChannelingECHARM(filePot,CLHEP::eV);
    fElectricFieldX =  new G4ChannelingECHARM(fileEFX,CLHEP::eV/CLHEP::m);
    fElectricFieldY =  new G4ChannelingECHARM(fileEFY,CLHEP::eV/CLHEP::m);
    fNucleiDensity =   new G4ChannelingECHARM(fileAtD,1.);
    fElectronDensity = new G4ChannelingECHARM(fileElD,1.);

    G4cout <<  filePot << G4endl;
    G4cout <<  fileEFX << G4endl;
    G4cout <<  fileEFY << G4endl;
    G4cout <<  fileAtD << G4endl;
    G4cout <<  fileElD << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMaterialData::SetFilenameElement(const G4String& fileName,std::string elementName){
    G4String filePot = fileName + "_pot.txt";
    G4String fileEFX = fileName + "_efx.txt";
    G4String fileEFY = fileName + "_efy.txt";
    G4String fileAtD = fileName + "_atd.txt";
    G4String fileElD = fileName + "_eld.txt";
    
    fPotentialElement[elementName] = new G4ChannelingECHARM(filePot,CLHEP::eV);
    fElectricFieldXElement[elementName] =  new G4ChannelingECHARM(fileEFX,CLHEP::eV/CLHEP::m);
    fElectricFieldYElement[elementName] =  new G4ChannelingECHARM(fileEFY,CLHEP::eV/CLHEP::m);
    fNucleiDensityElement[elementName] =   new G4ChannelingECHARM(fileAtD,1.);
    fElectronDensityElement[elementName] = new G4ChannelingECHARM(fileElD,1.);
    
    G4cout <<  filePot << G4endl;
    G4cout <<  fileEFX << G4endl;
    G4cout <<  fileEFY << G4endl;
    G4cout <<  fileAtD << G4endl;
    G4cout <<  fileElD << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMaterialData::SetBR(G4double val){
    fVectorR = new G4PhysicsLinearVector(0,DBL_MAX,2);
    fVectorR->PutValue(0,val);
    fVectorR->PutValue(1,val);
    bIsBent = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMaterialData::SetBR(const G4String& filename){
    std::ifstream vFileIn;
    int points;
    float maximum;
    vFileIn.open(filename);
    vFileIn >> points >> maximum;
    
    fVectorR = new G4PhysicsLinearVector(0,maximum * CLHEP::millimeter,points);
    double vTempX;
    for(G4int i0=0;i0<points; i0++){
        vFileIn >> vTempX;
        fVectorR->PutValue(i0,vTempX * CLHEP::meter);
    }
    bIsBent = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ChannelingMaterialData::~G4ChannelingMaterialData(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
