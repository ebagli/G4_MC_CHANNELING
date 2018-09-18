# Geant4 Channeling Simulation
Enrico Bagli - INFN (Italy) bagli@fe.infn.it

## 1. Introduction
The application simulates channeling in a straight or bent crystal via Geant4.
Charged particles impinging on a crystal can be captured into the channeling regime provided that their trajectories are aligned with crystalline planes or axes within the critical angle for channeling. Under such circumstance, particle motion is confined within the potential well between atomic planes or axes.

Channeling is simulated by including [DYNECHARM++](http://www.sciencedirect.com/science/article/pii/S0168583X1300308X) and the results of the [ECHARM software](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.81.026708) into the [Geant4 channeling package](https://link.springer.com/article/10.1140/epjc/s10052-014-2996-y). DYNECHARM++ allows the tracking of a relativistic charged particle inside a crystalline medium via the numerical integration of the classical equations of motion. The continuum potential approximation proposed by Lindhard is used. ECHARM is a software for the calculation of the potential and related physical quantities experienced by a particle traversing an aligned periodic complex atomic structure. The software uses classical physics equations and the expansion of periodic functions as a Fourier series. The Geant4 channeling package is an extension of the Geant4 toolkit for the simulation of the channeling effect.


## 2. Experimental Setup
The setup is composed by a crystal and two pairs of detectors to reveal the incoming and the outgoing angle of the particle from the crystal.

## 3. Installation
The application requires [Geant4 10.3](www.geant4.org).
Geant4 is a toolkit for the simulation of the passage of particles through matter.
Its areas of application include high energy, nuclear and accelerator physics, as well as studies in medical and space science.
The three main reference papers for Geant4 are:
- [Nuclear Instruments and Methods in Physics Research A 506 (2003) 250-303](http://www.sciencedirect.com/science/article/pii/S0168900203013688)
- [IEEE Transactions on Nuclear Science 53 No. 1 (2006) 270-278](http://ieeexplore.ieee.org/xpls/abs_all.jsp?isnumber=33833&arnumber=1610988&count=33&index=7)
- [Nuclear Instruments and Methods in Physics Research A 835 (2016) 186-225](http://www.sciencedirect.com/science/article/pii/S0168900216306957)

The source code can be found also on [GitHub](https://github.com/Geant4/geant4/tree/geant4-10.3-release).

After the installation of the Geant4 10.3 libraries, the code can be compiled via [CMake](www.cmake.org).

## 4. Usage
The application can be launched via the `channeling` executable and the macros in the `mac` directory, e.g. `./channeling mac/H8_PL01.mac`.


## 5. Output
The output of the application is a [ROOT](https://root.cern.ch) file with a tree that contains the particle incoming and outgoing angle from the crystal, the average electrid field and nuclei/electron densities experienced by the particle. Only the primary particles are detected.

- `angXin` - particle horizontal incoming angle at the crystal in µrad
- `angYin` - particle vertical incoming angle at the crystal in µrad
- `posXin` - particle horizontal hit position on the crystal surface in mm
- `posYin` - particle vertical hit position on the crystal surface in mm
- `angXout` - particle horizontal outgoing angle from the crystal in µrad
- `angYout` - particle horizontal outgoing angle from the crystal in µrad
- `efx` - average horizontal electric field experienced under channeling by the particle
- `efy` - average vertical electric field experienced under channeling by the particle
- `nud` - average nuclei density experienced under channeling by the particle
- `eld` - average electron density experienced under channeling by the particle
- `sx` - X spin component of the outgoing channeled particle
- `sy` - Y spin component of the outgoing channeled particle
- `sz` - Z spin component of the outgoing channeled particle
- `energy` - the channeled particle energy in MeV

## 6. Macros
The `.mac` files contain macro command that modify dynamically the setup.
The primary beam properties can be configured via [General Particle Source - GPS - commands](https://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch02s07.html) of Geant4.

The detector properties can be changed via the macro commands:
- `/mydet/setDetMaterial G4_Si` - set the detector material (`G4_Si`)
- `/mydet/setSize 50. 50. 0.6 mm` - set the detector sizes (`x y z [unit]`)
- `/mydet/setDistance1 -10087. mm` - set the detector 1 distance from the crystal (`z [unit]`)
- `/mydet/setDistance2  -332 mm` - set the detector 2 distance from the crystal (`z [unit]`)
- `/mydet/setDistance3  +493. mm` - set the detector 3 distance from the crystal (`z [unit]`)
- `/mydet/setDistance4  +10945. mm` - set the detector 4 distance from the crystal (`z [unit]`)

The crystal properties can be changed via the macro commands:
- `/xtal/setMaterial G4_Si` - set the crystal material (`G4_Si`)
- `/xtal/setSize 50. 50. 1. mm` - set the crystal sizes (`x y z [unit]`)
- `/xtal/setAngle 0. 0. 0. rad` - set the crystal rotation angle with respect to the incoming beam on the `z` direction, i.e. the `y` angle rotate the crystal on the horizontal plane, the `x` angle rotate the crystal on the vertical plane.
- `/mydet/setBR 10. 0. 0. m` - set the crystal bending radius (`x y z [unit]`). The `x` bending radius bend the crystal on the horizontal plane, the `y` and `z` bending radii do not modify the setup in this version of the application.
- `/mydet/setEC data/Si220pl` - set the filename for the potential, electric field, nuclei and electron density evaluated through the continuum approximation by Lindhard. The `pl` suffix means that the planar characteristics are loaded, `ax` that the axial characteristics are loaded.
