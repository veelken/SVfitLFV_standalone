#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneLikelihoodLFV.h"

#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
#include "TauAnalysis/SVfitStandalone/interface/LikelihoodFunctions.h"

using namespace svFitStandalone;

/// global function pointer for minuit or VEGAS
const SVfitStandaloneLikelihoodLFV* SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV = 0;
/// indicate first iteration for integration or fit cycle for debugging
static bool FIRST = true;

SVfitStandaloneLikelihoodLFV::SVfitStandaloneLikelihoodLFV(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, bool verbose) 
  : SVfitStandaloneLikelihood(measuredTauLeptons, measuredMET, covMET, verbose) 
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV::SVfitStandaloneLikelihoodLFV>:" << std::endl;
  }
  int numPromptLeptons = 0;
  int numMeasuredTauLeptons = 0;
  for ( std::vector<MeasuredTauLepton>::const_iterator measuredTauLepton = measuredTauLeptons.begin();
	measuredTauLepton != measuredTauLeptons.end(); ++measuredTauLepton ) {
    if ( measuredTauLepton->decayType() == svFitStandalone::kHadDecay ||
	 measuredTauLepton->decayType() == svFitStandalone::kLepDecay ) {
      measuredTauLepton_ = (*measuredTauLepton);
      ++numMeasuredTauLeptons;
    } else if ( measuredTauLepton->decayType() == svFitStandalone::kPrompt ) { 
      promptLepton_ = (*measuredTauLepton);
      ++numPromptLeptons;
    }
  }
  if ( !(numPromptLeptons == 1 && numMeasuredTauLeptons == 1) ) {
    std::cout << " >> ERROR : there must be exactly one measured tau lepton + one prompt lepton" << std::endl;
    errorCode_ |= LeptonNumber;
  }
  measuredTauLeptons_.clear();
  measuredTauLeptons_.push_back(measuredTauLepton_);
  measuredTauLeptons_.push_back(promptLepton_);
  // set global function pointer to this
  gSVfitStandaloneLikelihoodLFV = this;
}

const double*
SVfitStandaloneLikelihoodLFV::transformint(double* xPrime, const double* x, const double mtest, const int par) const
{
  bool skipLOG = false;
  double nunuMass, labframeXFrac, labframePhi;
  // used for determination of xFrac for second lepton
  double vmm = (measuredTauLepton_.p4() + promptLepton_.p4()).mass();
  labframeXFrac = square(vmm/mtest); // visible energy fraction x in labframe
  if ( verbose_ && FIRST ) {
    skipLOG = ((pow(vmm/mtest, 2)/labframeXFrac) > 1.);
    if ( !skipLOG ) {
      std::cout << "Boundary Check: labframeXFrac[" << measuredTauLepton_.decayType() << "]" << std::endl; 
    }
  }
  xPrime[kMaxNLLParams] = labframeXFrac;
  if ( measuredTauLepton_.decayType() == kLepDecay ) {
    nunuMass = x[0]; // nunu inv mass (can be const 0 for had tau decays)
    if ( verbose_ && FIRST && !skipLOG ) { 
      std::cout << "Boundary Check: nunuMass     [" << measuredTauLepton_.decayType() << "] -> 0" << std::endl; 
    }
    labframePhi = x[1]; // phi in labframe
    if ( verbose_ && FIRST && !skipLOG ) { 
      std::cout << "Boundary Check: labframePhi  [" << measuredTauLepton_.decayType() << "] -> 1" << std::endl; 
    }
  } else {
    nunuMass = 0.;
    labframePhi = x[0];
    if ( verbose_ && FIRST && !skipLOG ) { 
      std::cout << "Boundary Check: labframePhi  [" << measuredTauLepton_.decayType() << "] -> 0" << std::endl; 
    }
  }     
  double labframeVisMom = measuredTauLepton_.momentum(); // visible momentum in lab-frame
  double labframeVisEn  = measuredTauLepton_.energy();   // visible energy in lab-frame
  double visMass        = measuredTauLepton_.mass();     // vis mass
  // add protection against zero mass for visMass. If visMass is lower than the electron mass, set it
  // to the electron mass
  if ( visMass < 5.1e-4 ) { 
    visMass = 5.1e-4; 
  }    
  // momentum of visible decay products in tau lepton restframe
  double restframeVisMom     = svFitStandalone::pVisRestFrame(visMass, nunuMass, svFitStandalone::tauLeptonMass);
  // tau lepton decay angle in tau lepton restframe (as function of the energy ratio of visible decay products/tau lepton energy)
  double restframeDecayAngle = svFitStandalone::gjAngleFromX(labframeXFrac, visMass, restframeVisMom, labframeVisEn, svFitStandalone::tauLeptonMass);
  // tau lepton decay angle in labframe
  double labframeDecayAngle  = svFitStandalone::gjAngleToLabFrame(restframeVisMom, restframeDecayAngle, labframeVisMom);
  // tau lepton momentum in labframe
  double labframeTauMom      = svFitStandalone::motherMomentumLabFrame(visMass, restframeVisMom, restframeDecayAngle, labframeVisMom, svFitStandalone::tauLeptonMass);
  Vector labframeTauDir      = svFitStandalone::motherDirection(measuredTauLepton_.direction(), labframeDecayAngle, labframePhi).unit();
  // tau lepton four vector in labframe
  LorentzVector fittedTauLepton = svFitStandalone::motherP4(labframeTauDir, labframeTauMom, svFitStandalone::tauLeptonMass);
  LorentzVector fittedDiTauSystem = fittedTauLepton + promptLepton_.p4();
  // fill branch-wise nll parameters
  xPrime[ kNuNuMass1   ] = nunuMass;
  xPrime[ kVisMass1    ] = visMass;
  xPrime[ kDecayAngle1 ] = restframeDecayAngle;
  xPrime[ kNuNuMass2   ] = 0.;
  xPrime[ kVisMass2    ] = 0.;
  xPrime[ kDecayAngle2 ] = 0.;
  // subtract the visible part from it. The remainder is the pure neutrino part. Minus the the remainder is the estimate of the fittedMET
  Vector fittedMET = fittedDiTauSystem.Vect() - (measuredTauLeptons_[0].p() + measuredTauLeptons_[1].p()); 
  // fill event-wise nll parameters
  xPrime[ kDMETx   ] = measuredMET_.x() - fittedMET.x(); 
  xPrime[ kDMETy   ] = measuredMET_.y() - fittedMET.y();
  xPrime[ kMTauTau ] = mtest;

  if ( verbose_ && FIRST ) {
    std::cout << " >> input values for transformed variables: " << std::endl;
    std::cout << "    MET[x] = " <<  fittedMET.x() << " (fitted)  " << measuredMET_.x() << " (measured) " << std::endl; 
    std::cout << "    MET[y] = " <<  fittedMET.y() << " (fitted)  " << measuredMET_.y() << " (measured) " << std::endl; 
    std::cout << "    fittedDiTauSystem: [" 
	      << " px = " << fittedDiTauSystem.px() 
	      << " py = " << fittedDiTauSystem.py() 
	      << " pz = " << fittedDiTauSystem.pz() 
	      << " En = " << fittedDiTauSystem.energy() 
	      << " ]" << std::endl; 
    std::cout << " >> nll parameters after transformation: " << std::endl;
    std::cout << "    x[kNuNuMass1  ] = " << xPrime[kNuNuMass1  ] << std::endl;
    std::cout << "    x[kVisMass1   ] = " << xPrime[kVisMass1   ] << std::endl;
    std::cout << "    x[kDecayAngle1] = " << xPrime[kDecayAngle1] << std::endl;
    std::cout << "    x[kDMETx      ] = " << xPrime[kDMETx      ] << std::endl;
    std::cout << "    x[kDMETy      ] = " << xPrime[kDMETy      ] << std::endl;
    std::cout << "    x[kMTauTau    ] = " << xPrime[kMTauTau    ] << std::endl;
  }
  return xPrime;
}

double
SVfitStandaloneLikelihoodLFV::probint(const double* x, const double mtest, const int par) const 
{
  // in case of initialization errors don't start to do anything
  if ( error() ) {
    return 0.; 
  }
  double phiPenalty = 0.;
  double xPrime[kMaxNLLParams + 2];
  const double* xPrime_ptr = transformint(xPrime, x, mtest, par);
  if ( xPrime_ptr ) {
    return prob(xPrime_ptr, phiPenalty);
  } else{
    return 0.;
  }
}

const double*
SVfitStandaloneLikelihoodLFV::transform(double* xPrime, const double* x) const
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:transform(double*, const double*)>:" << std::endl;
  }
  LorentzVector fittedDiTauSystem;
  // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
  // for each tau decay
  double nunuMass       = x[ kMNuNu ];       // nunu inv mass (can be const 0 for had tau decays) 
  double labframeXFrac  = x[ kXFrac ];       // visible energy fraction x in labframe
  double labframePhi    = x[ kPhi   ];       // phi in labframe 
  double labframeVisMom = measuredTauLepton_.momentum(); // visible momentum in lab-frame
  double labframeVisEn  = measuredTauLepton_.energy();   // visible energy in lab-frame
  double visMass        = measuredTauLepton_.mass();     // vis mass

  // add protection against zero mass for visMass. If visMass is lower than the electron mass, set it
  // to the electron mass
  if ( visMass < 5.1e-4 ) { 
    visMass = 5.1e-4; 
  }    
  // momentum of visible decay products in tau lepton restframe
  double restframeVisMom     = pVisRestFrame(visMass, nunuMass, tauLeptonMass);
  // tau lepton decay angle in tau lepton restframe (as function of the energy ratio of visible decay products/tau lepton energy)
  double restframeDecayAngle = gjAngleFromX(labframeXFrac, visMass, restframeVisMom, labframeVisEn, tauLeptonMass);
  // tau lepton decay angle in labframe
  double labframeDecayAngle  = gjAngleToLabFrame(restframeVisMom, restframeDecayAngle, labframeVisMom);
  // tau lepton momentum in labframe
  double labframeTauMom      = motherMomentumLabFrame(visMass, restframeVisMom, restframeDecayAngle, labframeVisMom, tauLeptonMass);
  Vector labframeTauDir      = motherDirection(measuredTauLepton_.direction(), labframeDecayAngle, labframePhi).unit();
  // tau lepton four vector in labframe
  fittedDiTauSystem += motherP4(labframeTauDir, labframeTauMom, tauLeptonMass);
  // fill branch-wise nll parameters
  xPrime[ kNuNuMass1    ] = nunuMass;
  xPrime[ kVisMass1     ] = visMass;
  xPrime[ kDecayAngle1  ] = restframeDecayAngle;
  xPrime[ kMaxNLLParams ] = labframeXFrac;
 
  Vector fittedMET = fittedDiTauSystem.Vect() - (measuredTauLepton_.p() + promptLepton_.p()); 
  // fill event-wise nll parameters
  xPrime[ kDMETx   ] = measuredMET_.x() - fittedMET.x(); 
  xPrime[ kDMETy   ] = measuredMET_.y() - fittedMET.y();
  xPrime[ kMTauTau ] = fittedDiTauSystem.mass();

  if ( verbose_ && FIRST ) {
    std::cout << " >> input values for transformed variables: " << std::endl;
    std::cout << "    MET[x] = " <<  fittedMET.x() << " (fitted)  " << measuredMET_.x() << " (measured) " << std::endl; 
    std::cout << "    MET[y] = " <<  fittedMET.y() << " (fitted)  " << measuredMET_.y() << " (measured) " << std::endl; 
    std::cout << "    fittedDiTauSystem: [" 
	      << " px = " << fittedDiTauSystem.px() 
	      << " py = " << fittedDiTauSystem.py() 
	      << " pz = " << fittedDiTauSystem.pz() 
	      << " En = " << fittedDiTauSystem.energy() 
	      << " ]" << std::endl; 
    std::cout << " >> nll parameters after transformation: " << std::endl;
    std::cout << "    x[kNuNuMass1  ] = " << xPrime[ kNuNuMass1   ] << std::endl;
    std::cout << "    x[kVisMass1   ] = " << xPrime[ kVisMass1    ] << std::endl;
    std::cout << "    x[kDecayAngle1] = " << xPrime[ kDecayAngle1 ] << std::endl;
    std::cout << "    x[kDMETx      ] = " << xPrime[ kDMETx       ] << std::endl;
    std::cout << "    x[kDMETy      ] = " << xPrime[ kDMETy       ] << std::endl;
    std::cout << "    x[kMTauTau    ] = " << xPrime[ kMTauTau     ] << std::endl;
  }
  return xPrime;
}

double
SVfitStandaloneLikelihoodLFV::prob(const double* x) const 
{
  // in case of initialization errors don't start to do anything
  if ( error() ) { 
    return 0.;
  }
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:prob(const double*)>:" << std::endl;
  }
  ++idxObjFunctionCall_;
  if ( verbose_ && FIRST ) {
    std::cout << " >> ixdObjFunctionCall : " << idxObjFunctionCall_ << std::endl;  
    std::cout << " >> fit parameters before transformation: " << std::endl;
    std::cout << "    x[kXFrac1] = " << x[ kXFrac ] << std::endl;
    std::cout << "    x[kMNuNu1] = " << x[ kMNuNu ] << std::endl;
    std::cout << "    x[kPhi1  ] = " << x[ kPhi   ] << std::endl;
  }
  // prevent kPhi in the fit parameters (kFitParams) from trespassing the 
  // +/-pi boundaries
  double phiPenalty = 0.;
  if ( addPhiPenalty_ ) {
    if ( TMath::Abs(x[kPhi]) > TMath::Pi() ) {
      phiPenalty += (TMath::Abs(x[kPhi]) - TMath::Pi())*(TMath::Abs(x[kPhi]) - TMath::Pi());
    }
  }
  // xPrime are the transformed variables from which to construct the nll
  // transform performs the transformation from the fit parameters x to the 
  // nll parameters xPrime. prob is the actual combined likelihood. The
  // phiPenalty prevents the fit to converge to unphysical values beyond
  // +/-phi 
  double xPrime[kMaxNLLParams + 2];
  return prob(transform(xPrime, x), phiPenalty);
}

double 
SVfitStandaloneLikelihoodLFV::prob(const double* xPrime, double phiPenalty) const
{
  if ( verbose_ && FIRST ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:prob(const double*, double)>:" << std::endl;
  }
  // start the combined likelihood construction from MET
  double prob = probMET(xPrime[kDMETx], xPrime[kDMETy], covDet_, invCovMET_, metPower_, (verbose_&& FIRST));
  if ( verbose_ && FIRST ) {
    std::cout << "probMET         = " << prob << std::endl;
  }
  // add likelihood for the tau decay branch
  switch ( measuredTauLepton_.decayType() ) {
  case kHadDecay :
    prob *= probTauToHadPhaseSpace(
	      xPrime[kDecayAngle1], 
	      xPrime[kNuNuMass1], 
	      xPrime[kVisMass1], 
	      xPrime[kMaxNLLParams], 
	      addSinTheta_, 
	      (verbose_&& FIRST));
    if ( verbose_ && FIRST ) {
      std::cout << " *probTauToHad  = " << prob << std::endl;
    }
    break;
  case kLepDecay :
    prob *= probTauToLepMatrixElement(
	      xPrime[kDecayAngle1], 
	      xPrime[kNuNuMass1], 
	      xPrime[kVisMass1], 
	      xPrime[kMaxNLLParams], 
	      addSinTheta_, 
	      (verbose_&& FIRST));
    if ( verbose_ && FIRST ) {
      std::cout << " *probTauToLep  = " << prob << std::endl;
    }
    break;
  }
  // add additional logM term if configured such 
  if ( addLogM_ ) {
    if ( xPrime[kMTauTau] > 0. ) {
      prob *= (1.0/xPrime[kMTauTau]);
    }
    if ( verbose_ && FIRST ) {
      std::cout << " *1/mtautau     = " << prob << std::endl;
    }
  }
  if ( addDelta_ ) {
    prob *= (2.0*xPrime[kMaxNLLParams]/xPrime[kMTauTau]);
    if ( verbose_ && FIRST ) {
      std::cout << " *deltaDeriv.   = " << prob << std::endl;
    }
  }
  // add additional phiPenalty in case kPhi in the fit parameters 
  // (kFitParams) trespassed the physical boundaries from +/-pi 
  if ( phiPenalty > 0. ) {
    prob *= TMath::Exp(-phiPenalty);
    if ( verbose_ && FIRST ) {
      std::cout << "* phiPenalty   = " << prob << std::endl;
    }
  }
  // set FIRST to false after the first complete evaluation of the likelihood 
  FIRST = false;
  return prob;
}

void
SVfitStandaloneLikelihoodLFV::results(std::vector<LorentzVector>& fittedTauLeptons, const double* x) const
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:results(std::vector<LorentzVector>&, const double*)>:" << std::endl;
  }
  // map to local variables to be more clear on the meaning of the individual parameters. 
  double nunuMass       = x[ kMNuNu ];                   // nunu inv mass (can be const 0 for had tau decays) 
  double labframeXFrac  = x[ kXFrac ];                   // visible energy fraction x in labframe
  double labframePhi    = x[ kPhi   ];                   // phi in labframe 
  double labframeVisMom = measuredTauLepton_.momentum(); // visible momentum in lab-frame
  double labframeVisEn  = measuredTauLepton_.energy();   // visible energy in lab-frame
  double visMass        = measuredTauLepton_.mass();     // vis mass
  
  // momentum of visible decay products in tau lepton restframe
  double restframeVisMom     = pVisRestFrame(visMass, nunuMass, tauLeptonMass);
  // tau lepton decay angle in tau lepton restframe (as function of the energy ratio of visible decay products/tau lepton energy)
  double restframeDecayAngle = gjAngleFromX(labframeXFrac, visMass, restframeVisMom, labframeVisEn, tauLeptonMass);
  // tau lepton decay angle in labframe
  double labframeDecayAngle  = gjAngleToLabFrame(restframeVisMom, restframeDecayAngle, labframeVisMom);
  // tau lepton momentum in labframe
  double labframeTauMom      = motherMomentumLabFrame(visMass, restframeVisMom, restframeDecayAngle, labframeVisMom, tauLeptonMass);
  Vector labframeTauDir      = motherDirection(measuredTauLepton_.direction(), labframeDecayAngle, labframePhi).unit();
  // tau lepton four vector in labframe
  if ( fittedTauLeptons.size() >= 1 ) fittedTauLeptons[0] = motherP4(labframeTauDir, labframeTauMom, tauLeptonMass);
  else fittedTauLeptons.push_back(motherP4(labframeTauDir, labframeTauMom, tauLeptonMass));
  if ( fittedTauLeptons.size() >= 2 ) fittedTauLeptons[1] = promptLepton_.p4();
  else fittedTauLeptons.push_back(promptLepton_.p4());

  return;
}
