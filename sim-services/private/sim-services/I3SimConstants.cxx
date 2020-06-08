
#include <sim-services/I3SimConstants.h>

namespace I3SimConstants {

/** Electromagnetic and hadronic shower shape parameters from icecube/201210001 */
ShowerParameters::ShowerParameters(I3Particle::ParticleType type, double E, double density)
    : a(0.), b(0.), emScale(1.), emScaleSigma(0.)
{
    // protect against extremely low energies
    // NB: while Equation 4.11 of Leif's Masters' thesis is written in terms
    // of log10, we use the natural log here and divide the energy-scaling
    // coefficients (beta) below by ln(10) to compensate
    const double logE = std::max(0., std::log(E));
    const double Lrad = 0.358*(I3Units::g/I3Units::cm3)/density;
    
    const bool isElectron =
    (type==I3Particle::EMinus) ||
    (type==I3Particle::EPlus) ||
    (type==I3Particle::Brems) ||
    (type==I3Particle::DeltaE) ||
    (type==I3Particle::PairProd) ||
    (type==I3Particle::Gamma) ||
    (type==I3Particle::Pi0) ; // Pi0 decays to 2 gammas and produce EM showers

    bool isHadron =
    (type==I3Particle::Hadrons) ||
    (type==I3Particle::Neutron) ||
    (type==I3Particle::PiPlus) ||
    (type==I3Particle::PiMinus) ||
    (type==I3Particle::K0_Long) ||
    (type==I3Particle::KPlus) ||
    (type==I3Particle::KMinus) ||
    (type==I3Particle::PPlus) ||
    (type==I3Particle::PMinus) ||
    (type==I3Particle::K0_Short) ||
    (type==I3Particle::Eta) ||
    (type==I3Particle::Lambda) ||
    (type==I3Particle::SigmaPlus) ||
    (type==I3Particle::Sigma0) ||
    (type==I3Particle::SigmaMinus) ||
    (type==I3Particle::Xi0) ||
    (type==I3Particle::XiMinus) ||
    (type==I3Particle::OmegaMinus) ||
    (type==I3Particle::NeutronBar) ||
    (type==I3Particle::LambdaBar) ||
    (type==I3Particle::SigmaMinusBar) ||
    (type==I3Particle::Sigma0Bar) ||
    (type==I3Particle::SigmaPlusBar) ||
    (type==I3Particle::Xi0Bar) ||
    (type==I3Particle::XiPlusBar) ||
    (type==I3Particle::OmegaPlusBar) ||
    (type==I3Particle::DPlus) ||
    (type==I3Particle::DMinus) ||
    (type==I3Particle::D0) ||
    (type==I3Particle::D0Bar) ||
    (type==I3Particle::DsPlus) ||
    (type==I3Particle::DsMinusBar) ||
    (type==I3Particle::LambdacPlus) ||
    (type==I3Particle::WPlus) ||
    (type==I3Particle::WMinus) ||
    (type==I3Particle::Z0) ||
    (type==I3Particle::NuclInt);
    
    const bool isMuon =
    (type==I3Particle::MuMinus) ||
    (type==I3Particle::MuPlus);

    const bool isTau =
    (type==I3Particle::TauMinus) ||
    (type==I3Particle::TauPlus);
    
    if ((!isHadron) && (!isElectron) && (!isMuon) && (!isTau))
    {
        // if we don't know it but it has a pdg code,
        // it is probably a hadron..
        isHadron = true;
    }
    
    if (isElectron) {
        
        switch(type){
            default:
            case I3Particle::DeltaE:
            case I3Particle::Brems:
            case I3Particle::PairProd:
            case I3Particle::EMinus:
                a=2.01849+0.63176*logE;
                b=Lrad/0.63207;
                break;
            case I3Particle::EPlus:  // e+
                a=2.00035+0.63190*logE; 
                b=Lrad/0.63008;
                break;
            case I3Particle::Gamma:
            case I3Particle::Pi0:   // gamma, pi0
                a=2.83923+0.58209*logE; 
                b=Lrad/0.64526;
                break;  
        }   
    
    } else if (isHadron) {
        
        double E0(0.), m(0.), f0(1.), rms0(0.), gamma(0.);
        
        switch(type){
            default:
            case I3Particle::NuclInt:
            case I3Particle::Hadrons:
            case I3Particle::PiPlus:
                a=1.58357292+0.41886807*logE; 
                b=Lrad/0.33833116;
                E0=0.18791678;
                m =0.16267529;
                f0=0.30974123;
                rms0 =0.95899551;
                gamma=1.35589541;
                break;
            case I3Particle::PiMinus:
                a=1.69176636+0.40803489*logE; 
                b=Lrad/0.34108075;
                E0=0.19826506;
                m =0.16218006;
                f0=0.31859323;
                rms0 =0.94033488;
                gamma=1.35070162;
                break;
            case I3Particle::K0_Long:
                a=1.95948974+0.34934666*logE;
                b=Lrad/0.34535151;
                E0=0.21687243;
                m =0.16861530;
                f0=0.27724987;
                rms0 =1.00318874;
                gamma=1.37528605;
                break;
            case I3Particle::PPlus:
                a=1.47495778+0.40450398*logE;
                b=Lrad/0.35226706;
                E0=0.29579368;
                m =0.19373018;
                f0=0.02455403;
                rms0 =1.01619344;
                gamma=1.45477346;
                break;
            case I3Particle::Neutron:
                a=1.57739060+0.40631102*logE; 
                b=Lrad/0.35269455;
                E0=0.66725124;
                m =0.19263595;
                f0=0.17559033;
                rms0 =1.01414337;
                gamma=1.45086895;
                break;
            case I3Particle::PMinus:
                a=1.92249171+0.33701751*logE;
                b=Lrad/0.34969748;
                E0=0.29579368;
                m =0.19373018;
                f0=0.02455403;
                rms0 =1.01094637;
                gamma=1.50438415;
                break;
            }
            
            double e=std::max(2.71828183, E);
            emScale = 1.-pow(e/E0, -m)*(1.-f0);
            emScaleSigma = emScale*rms0*pow(log(e), -gamma);
        
    } else {
        log_fatal_stream("Particle type " << type << " is not a shower");
    }
    
    if (E < 1.*I3Units::GeV) b=0.; // this sets the cascade length to 0

}

}