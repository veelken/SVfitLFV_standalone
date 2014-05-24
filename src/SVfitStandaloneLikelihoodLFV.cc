#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneLikelihoodLFV.h"

#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
#include "TauAnalysis/SVfitStandalone/interface/LikelihoodFunctions.h"

using namespace svFitStandalone;

/// global function pointer for minuit or VEGAS
const SVfitStandaloneLikelihoodLFV* SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV = 0;
/// indicate first iteration for integration or fit cycle for debugging
static bool FIRST = true;

SVfitStandaloneLikelihoodLFV::SVfitStandaloneLikelihoodLFV(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, bool verbose) 
  : SVfitStandaloneLikelihood(measuredTauLeptons, measuredMET, covMET, verbose),    
    lutVisMassRes_(0),
    lutVisPtRes_(0)
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV::SVfitStandaloneLikelihoodLFV>:" << std::endl;
  }
  int numPromptLeptons = 0;
  int numMeasuredTauLeptons = 0;
  for ( std::vector<MeasuredTauLepton>::const_iterator measuredTauLepton = measuredTauLeptons.begin();
	measuredTauLepton != measuredTauLeptons.end(); ++measuredTauLepton ) {
    if ( measuredTauLepton->type() == svFitStandalone::kTauToHadDecay  ||
	 measuredTauLepton->type() == svFitStandalone::kTauToElecDecay ||
	 measuredTauLepton->type() == svFitStandalone::kTauToMuDecay   ) {
      measuredTauLepton_ = (*measuredTauLepton);
      ++numMeasuredTauLeptons;
    } else if ( measuredTauLepton->type() == svFitStandalone::kPrompt ) { 
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

void 
SVfitStandaloneLikelihoodLFV::shiftVisMassAndPt(bool value, const TH1* lutVisMassRes, const TH1* lutVisPtRes)
{
  std::cout << "<SVfitStandaloneLikelihoodLFV::shiftVisMassAndPt>:" << std::endl;
  if ( value ) std::cout << "enabling shiftVisMassAndPt: lutVisMassRes = " << lutVisMassRes << ", lutVisPtRes = " << lutVisPtRes << std::endl;
  else std::cout << "disabling shiftVisMassAndPt" << lutVisPtRes << std::endl;
  shiftVisMassAndPt_ = value;
  if ( shiftVisMassAndPt_ ) {
    lutVisMassRes_ = lutVisMassRes;
    lutVisPtRes_   = lutVisPtRes;
  }
}

const double*
SVfitStandaloneLikelihoodLFV::transform(double* xPrime, const double* x, bool fixToMtest, double mtest) const
{
#ifdef SVFIT_DEBUG 
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:transform(double*, const double*)>:" << std::endl;
  }
#endif 
  LorentzVector fittedDiTauSystem = promptLepton_.p4();
  // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
  // for each tau decay
  double nunuMass, labframeXFrac, labframePhi;
  double visMass_unshifted = measuredTauLepton_.mass();
  double visMass = visMass_unshifted; // visible momentum in lab-frame
  double labframeVisMom_unshifted = measuredTauLepton_.momentum(); 
  double labframeVisMom = labframeVisMom_unshifted; // visible momentum in lab-frame
  double labframeVisEn = measuredTauLepton_.energy(); // visible energy in lab-frame
  if ( measuredTauLepton_.type() == kTauToElecDecay || measuredTauLepton_.type() == kTauToMuDecay ) {
    labframeXFrac = x[kXFrac];
    nunuMass = x[kMNuNu];
    labframePhi = x[kPhi];
  } else {
    labframeXFrac = x[kXFrac];
    nunuMass = 0.;
    labframePhi = x[kPhi];
    if ( shiftVisMassAndPt_ ) {
      visMass = x[kVisMassShifted];
      double shiftInv = 1. + x[kRecTauPtDivGenTauPt];
      double shift = ( shiftInv > 1.e-1 ) ?
        (1./shiftInv) : 1.e+1;
      labframeVisMom *= shift;
      //visMass *= shift; // CV: take mass and momentum to be correlated
      //labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
      labframeVisEn *= shift;
    }
  }
#ifdef SVFIT_DEBUG 
  if ( verbose_ ) {
    std::cout << " labframeXFrac = " << labframeXFrac << std::endl;
    std::cout << " nunuMass = " << nunuMass << std::endl;
    std::cout << " labframePhi = " << labframePhi << std::endl;
    std::cout << " visMass = " << visMass << std::endl;
    std::cout << " labframeVisMom = " << labframeVisMom << std::endl;
    std::cout << " labframeVisEn = " << labframeVisEn << std::endl;    
  }
#endif 
  bool isValidSolution = true;
  // add protection against unphysical mass of visible tau decay products
  if ( visMass < electronMass || visMass > tauLeptonMass ) { 
    isValidSolution = false;
  }  
  // add protection against unphysical visible energy fractions
  //if ( !(labframeXFrac >= 0. && labframeXFrac <= 1.) ) {
  //  isValidSolution = false;
  //}
  // CV: do not spend time on unphysical solutions: returning 0 pointer will lead to 0 evaluation of prob
  if ( !isValidSolution ) {
    return 0;
  }
  double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
  double enTau_lab = labframeVisEn/labframeXFrac;
  if ( (enTau_lab*enTau_lab) < tauLeptonMass2 ) {
    enTau_lab = tauLeptonMass;
    isValidSolution = false;
  }  
  double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
  double gamma = enTau_lab/tauLeptonMass;
  double beta = TMath::Sqrt(1. - 1./(gamma*gamma));
  double pVis_parl_rf = -beta*gamma*labframeVisEn + gamma*TMath::Cos(gjAngle_lab)*labframeVisMom;
  double pVis_perp = labframeVisMom*TMath::Sin(gjAngle_lab);
  double gjAngle_rf = TMath::ATan2(pVis_perp, pVis_parl_rf);
  Vector p3Tau_unit = motherDirection(measuredTauLepton_.direction(), gjAngle_lab, labframePhi);
  LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
#ifdef SVFIT_DEBUG 
  if ( verbose_ ) {
    std::cout << " gjAngle_lab = " << gjAngle_lab << std::endl;
    std::cout << " enTau_lab = " << enTau_lab << std::endl;
    std::cout << " pTau_lab = " << pTau_lab << std::endl;
    std::cout << " gamma = " << gamma << std::endl;
    std::cout << " beta = " << beta << std::endl;
    std::cout << " pVis_parl_rf = " << pVis_parl_rf << std::endl;
    std::cout << " pVis_perp = " << pVis_perp << std::endl;
    std::cout << " gjAngle_rf = " << gjAngle_rf << std::endl;
    std::cout << " p3Tau_unit: Px = " << p3Tau_unit.X() << ", Py = " << p3Tau_unit.Y() << ", Pz = " << p3Tau_unit.X() << std::endl;
    std::cout << " p4Tau_lab: Pt = " << p4Tau_lab.pt() << ", eta = " << p4Tau_lab.eta() << ", phi = " << p4Tau_lab.phi() << ", mass = " << p4Tau_lab.mass() << std::endl;
  }
#endif 
  fittedDiTauSystem += p4Tau_lab;
  // fill branch-wise nll parameters
  xPrime[ kNuNuMass1            ] = nunuMass;
  xPrime[ kVisMass1             ] = visMass;
  xPrime[ kDecayAngle1          ] = gjAngle_rf;
  xPrime[ kDeltaVisMass1        ] = visMass_unshifted - visMass;
  xPrime[ kRecTauPtDivGenTauPt1 ] = ( labframeVisMom > 0. ) ? (labframeVisMom_unshifted/labframeVisMom) : 1.e+3;
  xPrime[ kMaxNLLParams         ] = labframeXFrac;
  xPrime[ kMaxNLLParams + 1     ] = isValidSolution;
 
  Vector fittedMET = fittedDiTauSystem.Vect() - (measuredTauLepton_.p() + promptLepton_.p());
  // fill event-wise nll parameters
  xPrime[ kDMETx ] = measuredMET_.x() - fittedMET.x(); 
  xPrime[ kDMETy ] = measuredMET_.y() - fittedMET.y();
  if ( fixToMtest ) xPrime[ kMTauTau ] = mtest;       // CV: evaluate delta-function derrivate in case of VEGAS integration for nominal test mass,
  else xPrime[ kMTauTau ] = fittedDiTauSystem.mass(); //     not for fitted mass, to improve numerical stability of integration (this is what the SVfit plugin version does)
#ifdef SVFIT_DEBUG 
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
    std::cout << "    x[kNuNuMass    ] = " << xPrime[ kNuNuMass1    ] << std::endl;
    std::cout << "    x[kVisMass     ] = " << xPrime[ kVisMass1     ] << std::endl;
    std::cout << "    x[kDecayAngle  ] = " << xPrime[ kDecayAngle1  ] << std::endl;
    std::cout << "    x[kDMETx       ] = " << xPrime[ kDMETx        ] << std::endl;
    std::cout << "    x[kDMETy       ] = " << xPrime[ kDMETy        ] << std::endl;
    std::cout << "    x[kMTauTau     ] = " << xPrime[ kMTauTau      ] << std::endl;
    std::cout << "    x[kMaxNLLParams] = " << xPrime[ kMaxNLLParams ] << std::endl;
  }
#endif 
  return xPrime;
}

double
SVfitStandaloneLikelihoodLFV::prob(const double* x, bool fixToMtest, double mtest) const 
{
  // in case of initialization errors don't start to do anything
  if ( error() ) { 
    return 0.;
  }
#ifdef SVFIT_DEBUG 
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:prob(const double*)>:" << std::endl;
  }
#endif 
  ++idxObjFunctionCall_;
#ifdef SVFIT_DEBUG 
  if ( verbose_ && FIRST ) {
    std::cout << " >> ixdObjFunctionCall : " << idxObjFunctionCall_ << std::endl;  
    std::cout << " >> fit parameters before transformation: " << std::endl;
    std::cout << "    x[kXFrac]               = " << x[ kXFrac               ] << std::endl;
    std::cout << "    x[kMNuNu]               = " << x[ kMNuNu               ] << std::endl;
    std::cout << "    x[kPhi  ]               = " << x[ kPhi                 ] << std::endl;
    std::cout << "    x[kVisMassShifted]      = " << x[ kVisMassShifted      ] << std::endl;
    std::cout << "    x[kRecTauPtDivGenTauPt] = " << x[ kRecTauPtDivGenTauPt ] << std::endl;
  }
#endif 
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
  const double* xPrime_ptr = transform(xPrime, x, fixToMtest, mtest);
  if ( xPrime_ptr ) {
    return prob(xPrime_ptr, phiPenalty);
  } else {
    return 0.;
  }
}

double 
SVfitStandaloneLikelihoodLFV::prob(const double* xPrime, double phiPenalty) const
{
#ifdef SVFIT_DEBUG 
  if ( verbose_ && FIRST ) {
    std::cout << "<SVfitStandaloneLikelihoodLFV:prob(const double*, double)>:" << std::endl;
  }
#endif 
  if ( requirePhysicalSolution_ && xPrime[ kMaxNLLParams + 1 ] < 0.5 ) return 0.;
  // start the combined likelihood construction from MET
  double prob = probMET(xPrime[kDMETx], xPrime[kDMETy], covDet_, invCovMET_, metPower_, (verbose_&& FIRST));
#ifdef SVFIT_DEBUG 
  if ( verbose_ && FIRST ) {
    std::cout << "probMET         = " << prob << std::endl;
  }
#endif 
  // add likelihood for the tau decay branch
  switch ( measuredTauLepton_.type() ) {
  case kTauToHadDecay :
    prob *= probTauToHadPhaseSpace(
	      xPrime[kDecayAngle1], 
	      xPrime[kNuNuMass1], 
	      xPrime[kVisMass1], 
	      xPrime[kMaxNLLParams], 
	      addSinTheta_, 
	      (verbose_&& FIRST));
    if ( shiftVisMassAndPt_ ) {
      prob *= probVisMassAndPtShift(
		xPrime[kDeltaVisMass1], 
		xPrime[kRecTauPtDivGenTauPt1], 
		lutVisMassRes_, lutVisPtRes_, 
		(verbose_&& FIRST));
      }
#ifdef SVFIT_DEBUG 
    if ( verbose_ && FIRST ) {
      std::cout << " *probTauToHad  = " << prob << std::endl;
    }
#endif 
    break;
  case kTauToElecDecay :
  case kTauToMuDecay :  
    prob *= probTauToLepMatrixElement(
	      xPrime[kDecayAngle1], 
	      xPrime[kNuNuMass1], 
	      xPrime[kVisMass1], 
	      xPrime[kMaxNLLParams], 
	      addSinTheta_, 
	      (verbose_&& FIRST));
#ifdef SVFIT_DEBUG 
    if ( verbose_ && FIRST ) {
      std::cout << " *probTauToLep  = " << prob << std::endl;
    }
#endif 
    break;
  }
  // add additional logM term if configured such 
  if ( addLogM_ ) {
    if ( xPrime[kMTauTau] > 0. ) {
      prob *= (1.0/xPrime[kMTauTau]);
    }
#ifdef SVFIT_DEBUG 
    if ( verbose_ && FIRST ) {
      std::cout << " *1/mtautau     = " << prob << std::endl;
    }
#endif 
  }
  if ( addDelta_ ) {
    prob *= (2.0*xPrime[kMaxNLLParams]/xPrime[kMTauTau]);
#ifdef SVFIT_DEBUG 
    if ( verbose_ && FIRST ) {
      std::cout << "xPrime[kMaxNLLParams] = " << xPrime[kMaxNLLParams] << std::endl;
      std::cout << "xPrime[kMTauTau] = " << xPrime[kMTauTau] << std::endl;
      std::cout << " *deltaDeriv.   = " << prob << std::endl;
    }
#endif 
  }
  // add additional phiPenalty in case kPhi in the fit parameters 
  // (kFitParams) trespassed the physical boundaries from +/-pi 
  if ( phiPenalty > 0. ) {
    prob *= TMath::Exp(-phiPenalty);
#ifdef SVFIT_DEBUG 
    if ( verbose_ && FIRST ) {
      std::cout << "* phiPenalty   = " << prob << std::endl;
    }
#endif 
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
  double nunuMass                 = x[ kMNuNu ];       // nunu inv mass (can be const 0 for had tau decays) 
  double labframeXFrac            = x[ kXFrac ];       // visible energy fraction x in labframe
  double labframePhi              = x[ kPhi   ];       // phi in labframe 
  double visMass                  = measuredTauLepton_.mass(); 
  double labframeVisMom_unshifted = measuredTauLepton_.momentum(); 
  double labframeVisMom           = labframeVisMom_unshifted; // visible momentum in lab-frame
  double labframeVisEn            = measuredTauLepton_.energy(); // visible energy in lab-frame    
  if ( measuredTauLepton_.type() == kTauToHadDecay && shiftVisMassAndPt_ ) {
    visMass = x[kVisMassShifted];
    double shiftInv = 1. + x[kRecTauPtDivGenTauPt];
    double shift = ( shiftInv > 1.e-1 ) ?
      (1./shiftInv) : 1.e+1;
    labframeVisMom *= shift;
    //visMass *= shift; // CV: take mass and momentum to be correlated
    //labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
    labframeVisEn *= shift;
  }
  if ( visMass < 5.1e-4 ) { 
    visMass = 5.1e-4; 
  } 
  bool isValidSolution = true;
  double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
  double enTau_lab = labframeVisEn/labframeXFrac;
  double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
  Vector p3Tau_unit = motherDirection(measuredTauLepton_.direction(), gjAngle_lab, labframePhi);
  LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
  // tau lepton four vector in labframe
  if ( fittedTauLeptons.size() >= 1 ) fittedTauLeptons[0] = promptLepton_.p4();
  else fittedTauLeptons.push_back(promptLepton_.p4());
  if ( fittedTauLeptons.size() >= 2 ) fittedTauLeptons[1] = p4Tau_lab;
  else fittedTauLeptons.push_back(p4Tau_lab);
}
