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

#include "G4Channeling.hh"

#include "Randomize.hh"

#include "G4ChannelingTrackData.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4ChargeState.hh"
#include "G4LambdacPlus.hh"
#include "G4AntiLambdacPlus.hh"
#include "G4XibMinus.hh"
#include "G4AntiXibMinus.hh"
#include "G4XiMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiOmegaMinus.hh"
#include "G4SigmaPlus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

G4Channeling::G4Channeling():
G4VDiscreteProcess("channeling"),
fChannelingID(-1),
fTimeStepMin(0.),
fTimeStepMax(0.),
fTransverseVariationMax(2.E-2 * CLHEP::angstrom),
k010(G4ThreeVector(0.,1.,0.)){
    fChannelingID = G4PhysicsModelCatalog::GetIndex("channeling");
    if(fChannelingID == -1){
        fChannelingID = G4PhysicsModelCatalog::Register("channeling");
    }
    fSpin = G4ThreeVector(0.,0.,0.);
#ifdef bSaveTrajectoryToFileForChanneling
    outfile.open("bSaveTrajectoryToFileForChanneling.txt", std::ios_base::app);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Channeling::~G4Channeling(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Channeling::G4Channeling(G4Channeling& right):G4VDiscreteProcess(right){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingTrackData* G4Channeling::GetTrackData(const G4Track& aTrack){
    G4ChannelingTrackData* trackdata =
        (G4ChannelingTrackData*)(aTrack.GetAuxiliaryTrackInformation(fChannelingID));
    if(trackdata == nullptr){
        trackdata = new G4ChannelingTrackData();
        aTrack.SetAuxiliaryTrackInformation(fChannelingID,trackdata);
    }
    return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(const G4Track& aTrack,
                         G4ThreeVector& pos,
                         G4ThreeVector& out){
    out = G4ThreeVector((GetMatData(aTrack)->GetEFX()->GetEC(pos)),
                        (GetMatData(aTrack)->GetEFY()->GetEC(pos)),
                        0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::PosToLattice(G4StepPoint* step,G4ThreeVector& pos){
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetTouchable());
    G4Box* box = (G4Box*) step->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();

    pos -= theTouchable->GetTranslation();
    pos = ((*theTouchable->GetRotation()).inverse())(pos);
    pos += G4ThreeVector(0.,0.,box->GetZHalfLength());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::UpdateParameters(const G4Track& aTrack){

    G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
    
    G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
    G4StepPoint* preStepPoint = aTrack.GetStep()->GetPreStepPoint();
    
    G4double mass = aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass();
    G4double gamma = aTrack.GetTotalEnergy()/mass;
    G4double beta = aTrack.GetVelocity()/CLHEP::c_light;
    G4double m_e = CLHEP::electron_mass_c2;
    G4double Tmax = 2. * m_e * gamma * gamma * beta * beta;
    Tmax /= (1. + 2.*gamma*m_e/mass + m_e*m_e/mass/mass);
    G4double T = aTrack.GetStep()->GetTotalEnergyDeposit() - aTrack.GetStep()->GetNonIonizingEnergyDeposit();
    G4double mommod = aTrack.GetStep()->GetPreStepPoint()->GetMomentum().mag();
    G4double theta_sc_el = sqrt(2.*m_e*T*(1.-T/Tmax))/mommod;
   
    /*
    G4cout << 2.*m_e*T/CLHEP::eV << G4endl;
    G4cout << T/Tmax << G4endl;
    G4cout << mommod/CLHEP::eV << G4endl;
    G4cout << aTrack.GetStep()->GetTotalEnergyDeposit()/CLHEP::eV << G4endl;
    G4cout << aTrack.GetStep()->GetNonIonizingEnergyDeposit()/CLHEP::eV << G4endl;

    while(!getchar());
    */
    
    G4ThreeVector posPost = postStepPoint->GetPosition();
    aLCV->RotateToLattice(posPost);
    G4ThreeVector posPre = preStepPoint->GetPosition();
    aLCV->RotateToLattice(posPre);
    
    G4double integrationLimit = fabs(posPost.z() - posPre.z());
    
    //G4cout << aTrack.GetWeight() << G4endl;
    //while(!getchar());
    if(integrationLimit > 0.){
        //----------------------------------------
        // Check if the crystal is bent
        //----------------------------------------
        G4bool isBent = GetMatData(aTrack)->IsBent();

        //----------------------------------------
        // Get the momentum in the world reference
        // frame and rotate to the solid reference frame
        //----------------------------------------
        G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
        G4ThreeVector momWorld = aTrack.GetStep()->GetPreStepPoint()->GetMomentum();
        G4ThreeVector mom = (*theTouchable->GetRotation())(momWorld);
        
        //----------------------------------------
        // Get the momentum in the solid reference
        // frame and rotate to the crystal reference frame
        //----------------------------------------
        aLCV->RotateToLattice(mom);

        //----------------------------------------
        // Get the momentum in the crystal reference
        // frame and rotate to the reference frame
        // solidal to the bent planes
        //----------------------------------------
        if(isBent){
            PosToLattice(preStepPoint,posPre);
            G4ThreeVector axis010 = (*theTouchable->GetRotation())(k010);
            mom.rotate(axis010,-posPre.z()/GetMatData(aTrack)->GetBR(posPre).x());
        }

        //----------------------------------------
        // Take the position stored in the track data.
        // If the particle enters the crystal,
        // the position in the channel is randomly
        // generated using a uniform distribution
        //----------------------------------------
        G4ThreeVector pos;
        if(GetTrackData(aTrack)->GetPosCh().x() == DBL_MAX){
            G4double posX = G4UniformRand() * GetMatData(aTrack)->GetPot()->GetIntSp(0);
            G4double posY = G4UniformRand() * GetMatData(aTrack)->GetPot()->GetIntSp(1);
            pos = G4ThreeVector(posX,posY,0.);
        }
        else{
            pos = GetTrackData(aTrack)->GetPosCh();
        }

        G4double step=0., stepTot=0.;
        G4double nud =0., eld    =0.;
        G4double efx =0., efy    =0.;
        G4double nud_temp =0., eld_temp    =0.;

        //G4double beta = aTrack.GetVelocity()/CLHEP::c_light;
        G4double Z = GetParticleDefinition(aTrack)->GetPDGCharge();
        
        const G4double oneSixth = 1./6.;
        G4ThreeVector posk1,posk2,posk3,posk4,posk5,posk6;
        G4ThreeVector momk1,momk2,momk3,momk4,momk5,momk6;
        G4ThreeVector pos_temp, efxy;

        //----------------------------------------
        // Initialize variables for polarization modification
        //----------------------------------------
        const G4DynamicParticle*  pParticle  = aTrack.GetDynamicParticle() ;
        const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
        
        fSpin                   = aTrack.GetPolarization();
        //G4double mass
        G4double magMoment, spin;
        G4double muB;
        G4double g_BMT = 2.;
        G4double anomaly = 0.;
        G4double omegac = 0.;
        //G4double gamma = 0.;
        G4double charge = 0.;
        G4double spinmag = fSpin.mag();
        G4double elMoment = 0.;
        G4double d = 0.05;
        G4double eta = 0.;
        
        if(spinmag!=0.){
            //mass      = pParticleDef->GetPDGMass() ;
            charge    = pParticleDef->GetPDGCharge();
            magMoment = pParticleDef->GetPDGMagneticMoment();
            spin      = pParticleDef->GetPDGSpin();
            //gamma     = aTrack.GetTotalEnergy()/mass;
            eta       = 0.;
        
            omegac    = (CLHEP::eplus/mass)*CLHEP::c_light;
            muB       = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/(mass/CLHEP::c_squared);
            elMoment  = 2.* d * muB * spin;
            G4double muN = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/(938.2720813*CLHEP::MeV/CLHEP::c_squared);

            
            if(pParticleDef == G4LambdacPlus::Definition() || pParticleDef == G4AntiLambdacPlus::Definition()){
                anomaly = +0.3;
                //elMoment = 1.977E-3 * muB;
                magMoment = 2.*(anomaly + 1.)*spin*muB;
            }
            else if(pParticleDef == G4XibMinus::Definition() || pParticleDef == G4AntiXibMinus::Definition()){
                anomaly = +1.38;
                //elMoment = 0.;
                magMoment = 2.*(anomaly + 1.)*spin*muB;
            }
            else if(pParticleDef == G4SigmaPlus::Definition() || pParticleDef == G4AntiSigmaPlus::Definition()){
                //anomaly = +1.38;
                //elMoment = 0.;
                magMoment = 2.458 * muN;
                g_BMT = std::abs(magMoment)/muB/spin;
                anomaly   = -(g_BMT - 2.)/2.;

            }
            else{
                if(spin != 0. && muB != 0.){
                    g_BMT = std::abs(magMoment)/muB/spin;
                }
                anomaly   = (g_BMT + 2.)/2.;
            }
            
            
            if(spin != 0. && muB != 0.){
                //eta = 0.5 * d / spin;
                eta = d;
            }
            else{
                eta = 0.;
            }
            
            if(pParticleDef->GetLeptonNumber() > 0. || pParticleDef->GetBaryonNumber() > 0. ){
                anomaly = -anomaly;
            }
            
            
            /*
            G4cout << "Expected precession (to be multiplied by channeling deflection angle): " << gamma*anomaly << G4endl;
            G4cout << "Expected precession (to be multiplied by channeling deflection angle in rad): " << gamma*anomaly*180./3.14 << G4endl;
            G4cout << "Magnetic Dipole Moment [unit of mu_N]: " << magMoment / muN << G4endl;
            G4cout << "Electric Dipole Moment [unit of mu_N]: " << elMoment / muN << G4endl;
            G4cout << "Anomaly: " << anomaly << G4endl;
            G4cout << "g: " << g_BMT << G4endl;
            G4cout << "eta: " << eta << G4endl;
            G4cout << "Mass: " << mass/CLHEP::GeV << G4endl;
            G4cout << "Spin: " << spin << G4endl;
            G4cout << "Spin: " << pParticleDef->GetPDGIsospin() << G4endl;
            G4cout << "muB: " << muB/muN << G4endl;
            while(!getchar());
            */
            
        }

        // Scattering from energy deposition in electrons
        
         if(theta_sc_el!=0.0){
            G4double eld_sc_el = GetTrackData(aTrack)->GetElD();
            G4double rot_angle = G4UniformRand()*CLHEP::twopi;
            G4double rot_mod   = G4RandGauss::shoot(0.,theta_sc_el*eld_sc_el);
            mom.rotateX(cos(rot_angle)*rot_mod);
            mom.rotateY(sin(rot_angle)*rot_mod);
            /*
            G4cout << rot_angle << G4endl;
            G4cout << rot_mod << G4endl;
            G4cout << theta_sc_el << G4endl;
            G4cout << eld_sc_el << G4endl;
            while(!getchar());
            */
        }
        
        //
        
        do{
            //----------------------------------------
            // Limit the variable step length for the
            // integration via the selected algorithm
            // and update variables for the integration
            //----------------------------------------

            UpdateIntegrationStep(aTrack,mom,step);
            if(step + stepTot > integrationLimit){
                step = integrationLimit - stepTot;
            }

            //----------------------------------------
            // Function integration algorithm
            // 4th Order Runge-Kutta
            //----------------------------------------
            
            GetEF(aTrack,pos,efxy);
            posk1 = step / mom.z() * mom;
            momk1 = step / beta * Z * efxy;
            if(isBent) momk1.setX(momk1.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos)).x());
            
            GetEF(aTrack,pos_temp = pos + posk1 * 0.5,efxy);
            posk2 = step / mom.z() * (mom + momk1 * 0.5);
            momk2 = step / beta * Z * efxy;
            if(isBent) momk2.setX(momk2.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());

            GetEF(aTrack,pos_temp = pos + posk2 * 0.5,efxy);
            posk3 = step / mom.z() * (mom + momk2 * 0.5);
            momk3 = step / beta * Z * efxy;
            if(isBent) momk3.setX(momk3.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());
            
            GetEF(aTrack,pos_temp = pos + posk3,efxy);
            posk4 = step / mom.z() * (mom + momk3);
            momk4 = step / beta * Z * efxy;
            if(isBent) momk4.setX(momk4.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());

            pos = pos + oneSixth * (posk1 + 2.*posk2 + 2.*posk3 + posk4);
            mom = mom + oneSixth * (momk1 + 2.*momk2 + 2.*momk3 + momk4);
            
            //----------------------------------------
            // Function integration algorithm
            // 2th Order Velocity-Verlet
            //----------------------------------------
           
            /*
            GetEF(aTrack,pos,efxy);
            posk1 = pos + (step * 0.5 / mom.z()) * mom;
            //momk1 = mom + step * 0.5 / betaZ * efxy;
            momk1 = mom;
            if(isBent) momk1.setX(momk1.x() - step * 0.5 * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos)).x());
            
            GetEF(aTrack,posk1,efxy);
            pos = pos + (step / momk1.z()) * momk1;
            //mom = mom + step / betaZ * efxy;
            mom = mom;
            if(isBent) mom.setX(mom.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(posk1)).x());
            */

            //----------------------------------------
            // Update the total step and the electron
            // and nuclei density experienced by
            // the particle during its motion
            //----------------------------------------

            stepTot += step;

            nud_temp = GetMatData(aTrack)->GetNuD()->GetEC(pos);
            eld_temp = GetMatData(aTrack)->GetElD()->GetEC(pos);
            
            if(nud_temp < 0.) {nud_temp = 0.;}
            if(eld_temp < 0.) {eld_temp = 0.;}

            nud += (step * nud_temp);
            eld += (step * eld_temp);

            efx += (step * GetMatData(aTrack)->GetEFX()->GetEC(pos));
            efy += (step * GetMatData(aTrack)->GetEFY()->GetEC(pos));
            
            //----------------------------------------
            // Spin Precession
            //----------------------------------------
            if (spinmag != 0.) {
                G4ThreeVector mu = mom;
                mu *= (1./mom.mag());
                G4ThreeVector ef;
                GetEF(aTrack,pos,ef);
                ef /= CLHEP::c_light;

                //----------------------------------------
                // Get the momentum in the reference frame
                // solidal to the bent planes and rotate
                // to the reference frame
                //----------------------------------------
                if(GetMatData(aTrack)->IsBent()){
                    G4ThreeVector axis010 = (*theTouchable->GetRotation())(k010);
                    mu.rotate(axis010,(posPre.z()+stepTot)/GetMatData(aTrack)->GetBR(posPre).x());
                    ef.rotate(axis010,(posPre.z()+stepTot)/GetMatData(aTrack)->GetBR(posPre).x());
                }
            
                //----------------------------------------
                // Get the momentum in the crystal reference
                // frame and rotate to the solid reference frame
                //----------------------------------------
            
                aLCV->RotateToSolid(mu);
                aLCV->RotateToSolid(ef);
            
                //----------------------------------------
                // Get the momentum in the solid reference
                // frame and rotate to the world reference frame
                //----------------------------------------
                mu = ((*theTouchable->GetRotation()).inverse())(mu);
                ef = ((*theTouchable->GetRotation()).inverse())(ef);
            
                G4ThreeVector BField = G4ThreeVector(0.,0.,0.);
            
                G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * mu);
                G4double ucb = (anomaly+1./gamma)/beta;
                G4double uce = anomaly + 1./(gamma+1.);
                G4double ude = beta*gamma/(1.+gamma)*(ef*mu);
            
                G4double pcharge;
                if (charge == 0.) pcharge = 1.;
                else pcharge = charge;
            
                G4ThreeVector dSpin(0.,0.,0.);
                dSpin =
                pcharge*omegac*( ucb*(fSpin.cross(BField))-udb*(fSpin.cross(mu))
                                - uce*(mu*(fSpin*ef) - ef*(fSpin*mu))
                                + eta/2.*(fSpin.cross(ef) - ude*(fSpin.cross(mu))
                                          + (mu*(fSpin*BField) - BField*(fSpin*mu)) ) );
                fSpin += dSpin*step;
                fSpin.setMag(spinmag);
            }
            //----------------------------------------


        } while(stepTot<integrationLimit);
        
        nud /= stepTot;
        eld /= stepTot;

        if(nud < 1.E-3) {nud = 1.E-3;}
        if(eld < 1.E-3) {eld = 1.E-3;}
        
        GetTrackData(aTrack)->SetNuD(nud);
        GetTrackData(aTrack)->SetElD(eld);

        GetTrackData(aTrack)->SetEFX(efx);
        GetTrackData(aTrack)->SetEFY(efy);
        
        GetTrackData(aTrack)->SetMomCh(mom);
        GetTrackData(aTrack)->SetPosCh(pos);
        
        
#ifdef bSaveTrajectoryToFileForChanneling
        G4int eID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        
        outfile << pos.x()/CLHEP::angstrom << "," << pos.y()/CLHEP::angstrom << "," << pos.z()/CLHEP::angstrom
        << "," << mom.x()/CLHEP::keV << "," << mom.y()/CLHEP::keV << "," << mom.z()/CLHEP::keV
        << "," << nud << "," << eld << "," << efx/stepTot
        << "," << GetMatData(aTrack)->GetPot()->GetEC(pos) << "," << eID << "," << stepTot/CLHEP::mm
        << "," << fSpin.x() << "," << fSpin.y()<< "," << fSpin.z() << G4endl;
#endif

        //G4cout << pos.x() << " " << pos.z() << " " << nud << " " << eld << G4endl;
        return true;
    }
    else{
        return false;
    }
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::
UpdateIntegrationStep(const G4Track& aTrack,
                      G4ThreeVector& mom,
                      G4double& step){
    
    if(mom.x() != 0.0 || mom.y() != 0.0){
        double xy2 = mom.x() * mom.x() + mom.y()*mom.y();
        
        if(xy2!=0.){
            step = std::fabs(fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy() / std::pow(xy2,0.5));
            if(step < fTimeStepMin) step = fTimeStepMin;
            else{
                fTimeStepMax = sqrt( fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy()
                                    / fabs(GetMatData(aTrack)->GetEFX()->GetMax()));
                
                if(step > fTimeStepMax) step = fTimeStepMax;
            }
        }
        else{
            step = fTimeStepMin;
        }
        
        return true;
    }
    else{
        step = fTimeStepMin;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Channeling::
GetMeanFreePath(const G4Track& aTrack,
                G4double, // previousStepSize
                G4ForceCondition* condition){
    
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step.
    // if it hasn't, return DBL_MAX
    //----------------------------------------
    
    *condition = Forced;
    
    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();
    
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true){
        G4double osc_per = GetOscillationPeriod(aTrack);
        fTimeStepMin = osc_per * 2.E-4;
        return osc_per * 0.01;
    }
    else{
        GetTrackData(aTrack)->Reset();
        return DBL_MAX;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4Channeling::
PostStepDoIt(const G4Track& aTrack,
             const G4Step&){
    
    //----------------------------------------
    // check if the volume has a lattice
    // and if the particle is in channeling.
    // If it is so, the particle is forced
    // to follow the channeling plane
    // direction. If the particle has
    // dechanneled or exited the crystal,
    // the outgoing angle is evaluated
    //----------------------------------------
    
    aParticleChange.Initialize(aTrack);
    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();
    
    
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true){

        G4bool bModifiedTraj = UpdateParameters(aTrack);

        if(bModifiedTraj==true){
            //----------------------------------------
            // Get the momentum in the reference frame
            // solidal to the bent planes and rotate
            // to the reference frame
            //----------------------------------------
            G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
            G4ThreeVector momCh = GetTrackData(aTrack)->GetMomCh();

            G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
            G4TouchableHistory* theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());

            if(GetMatData(aTrack)->IsBent()){
                G4ThreeVector posPost = postStepPoint->GetPosition();
                PosToLattice(postStepPoint,posPost);
                G4ThreeVector axis010 = (*theTouchable->GetRotation())(k010);
                momCh.rotate(axis010,posPost.z()/GetMatData(aTrack)->GetBR(posPost).x());
            }
            
            //----------------------------------------
            // Get the momentum in the crystal reference
            // frame and rotate to the solid reference frame
            //----------------------------------------

            aLCV->RotateToSolid(momCh);
            
            //----------------------------------------
            // Get the momentum in the solid reference
            // frame and rotate to the world reference frame
            //----------------------------------------
            G4ThreeVector mom = ((*theTouchable->GetRotation()).inverse())(momCh);

            aParticleChange.ProposeMomentumDirection(mom.unit());
            aParticleChange.ProposePolarization(fSpin);
        }
    }
    else{
        // if the volume has no lattice it resets the density factors
        GetTrackData(aTrack)->Reset();
    }
    //G4cout << "Weight " << aParticleChange.GetParentWeight() << G4endl;

    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
