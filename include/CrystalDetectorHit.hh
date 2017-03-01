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
// --------------------------------------------------------------
//

#ifndef CrystalDetectorHit_h
#define CrystalDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class CrystalDetectorHit : public G4VHit
{
public:
    CrystalDetectorHit();
    CrystalDetectorHit(G4int z);
    virtual ~CrystalDetectorHit();
    CrystalDetectorHit(
                        const CrystalDetectorHit &right);
    const CrystalDetectorHit& operator=(
                        const CrystalDetectorHit &right);
    int operator==(const CrystalDetectorHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
private:
    G4int fLayerID;
    G4double fNud;
    G4double fEld;
    G4double fEFX;
    G4double fEFY;
    G4double fStep;
    G4double fKinE;
    G4ThreeVector fWorldPos;
    
public:
    inline void SetLayerID(G4int z) { fLayerID = z; }
    inline G4int GetLayerID() const { return fLayerID; }
    inline void SetNud(G4double z) { fNud = z; }
    inline G4double GetNud() const { return fNud; }
    inline void SetEld(G4double z) { fEld = z; }
    inline G4double GetEld() const { return fEld; }
    inline void SetEFX(G4double z) { fEFX = z; }
    inline G4double GetEFX() const { return fEFX; }
    inline void SetEFY(G4double z) { fEFY = z; }
    inline G4double GetEFY() const { return fEFY; }
    inline void SetStep(G4double z) { fStep = z; }
    inline G4double GetStep() const { return fStep; }
    inline void SetKinE(G4double z) { fKinE = z; }
    inline G4double GetKinE() const { return fKinE; }
    inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
};

typedef G4THitsCollection<CrystalDetectorHit>
    CrystalDetectorHitsCollection;

#ifdef G4MULTITHREADED
extern G4ThreadLocal G4Allocator<CrystalDetectorHit>*
    CrystalDetectorHitAllocator;
#else
extern G4Allocator<CrystalDetectorHit>
    CrystalDetectorHitAllocator;
#endif

inline void* CrystalDetectorHit::operator new(size_t)
{
#ifdef G4MULTITHREADED
    if(!CrystalDetectorHitAllocator)
        CrystalDetectorHitAllocator =
        new G4Allocator<CrystalDetectorHit>;
    return (void *) CrystalDetectorHitAllocator->MallocSingle();
#else
    void* aHit;
    aHit = (void*)CrystalDetectorHitAllocator.MallocSingle();
    return aHit;
#endif
}

inline void CrystalDetectorHit::operator delete(void* aHit)
{
#ifdef G4MULTITHREADED
    CrystalDetectorHitAllocator->FreeSingle((
                    CrystalDetectorHit*) aHit);
#else
   CrystalDetectorHitAllocator.FreeSingle((
                    CrystalDetectorHit*) aHit);
#endif
}

#endif


