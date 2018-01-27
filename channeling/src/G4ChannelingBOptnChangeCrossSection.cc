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
#include "G4ChannelingBOptnChangeCrossSection.hh"
#include "G4InteractionLawPhysical.hh"
#include "G4BiasingProcessInterface.hh"


G4ChannelingBOptnChangeCrossSection::G4ChannelingBOptnChangeCrossSection(G4String name)
  : G4VBiasingOperation( name  ),
    fInteractionOccured( false )
{
  fBiasedExponentialLaw = new G4InteractionLawPhysical("LawForOperation"+name);
}

G4ChannelingBOptnChangeCrossSection::~G4ChannelingBOptnChangeCrossSection()
{
  if ( fBiasedExponentialLaw ) delete fBiasedExponentialLaw;
}

const G4VBiasingInteractionLaw* G4ChannelingBOptnChangeCrossSection::ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& )
{
  return fBiasedExponentialLaw;
}

void G4ChannelingBOptnChangeCrossSection::SetBiasedCrossSection(G4double xst)
{
  fBiasedExponentialLaw->SetPhysicalCrossSection( xst );
}

G4double G4ChannelingBOptnChangeCrossSection::GetBiasedCrossSection() const
{
  return fBiasedExponentialLaw->GetPhysicalCrossSection();
}

void G4ChannelingBOptnChangeCrossSection::Sample()
{
  fInteractionOccured = false;
  fBiasedExponentialLaw->Sample();
}

void G4ChannelingBOptnChangeCrossSection::UpdateForStep( G4double truePathLength )
{
  fBiasedExponentialLaw->UpdateForStep( truePathLength );
}

G4VParticleChange*  G4ChannelingBOptnChangeCrossSection::ApplyFinalStateBiasing(const G4BiasingProcessInterface* callingProcess,
                                                                                const G4Track* track,
                                                                                const G4Step* step,
                                                                                G4bool&  forceBiasedFinalState)
{
    //G4cout << "HERE" << G4endl;
    forceBiasedFinalState = true;
    return callingProcess->GetWrappedProcess()->PostStepDoIt( *track, *step);
}

