#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneAlgorithmLFV.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

#include <TGraphErrors.h>

namespace svFitStandalone
{
  void map_xVEGAS_LFV(const double* x, bool isLep, bool shiftVisMassAndPt, double mvis, double mtest, double* x_mapped)
  {
    double xFrac = TMath::Power(mvis/mtest, 2.);
    if ( isLep ) {
      x_mapped[kXFrac]               = xFrac;
      x_mapped[kMNuNu]               = x[0];
      x_mapped[kPhi]                 = x[1];
      x_mapped[kVisMassShifted]      = 0.;
      x_mapped[kRecTauPtDivGenTauPt] = 0.;
    } else {
      x_mapped[kXFrac]                 = xFrac;
      x_mapped[kMNuNu]                 = 0.;
      x_mapped[kPhi]                   = x[0];
      if ( shiftVisMassAndPt ) {
	x_mapped[kVisMassShifted]      = x[1];
	x_mapped[kRecTauPtDivGenTauPt] = x[2];
      }
    }
  }
  
  void map_xMarkovChain_LFV(const double* x, bool isLep, bool shiftVisMassAndPt, double* x_mapped)
  {
    if ( isLep ) {
      x_mapped[kXFrac]               = x[0];
      x_mapped[kMNuNu]               = x[1];
      x_mapped[kPhi]                 = x[2];
      x_mapped[kVisMassShifted]      = 0.;
      x_mapped[kRecTauPtDivGenTauPt] = 0.;
    } else {
      x_mapped[kXFrac]                 = x[0];
      x_mapped[kMNuNu]                 = 0.;
      x_mapped[kPhi]                   = x[1];
      if ( shiftVisMassAndPt ) {
	x_mapped[kVisMassShifted]      = x[2];
	x_mapped[kRecTauPtDivGenTauPt] = x[3];
      }
    }
  }
}

SVfitStandaloneAlgorithmLFV::SVfitStandaloneAlgorithmLFV(const std::vector<svFitStandalone::MeasuredTauLepton>& measuredTauLeptons, const svFitStandalone::Vector& measuredMET, const TMatrixD& covMET, 
							 unsigned int verbose) 
  : SVfitStandaloneAlgorithm(measuredTauLeptons, measuredMET, covMET, verbose),
    standaloneObjectiveFunctionAdapterVEGAS_LFV_(0),
    mcObjectiveFunctionAdapterLFV_(0),
    mcPtEtaPhiMassAdapterLFV_(0)
{ 
  // instantiate the combined likelihood
  nllLFV_ = new svFitStandalone::SVfitStandaloneLikelihoodLFV(measuredTauLeptons, measuredMET, covMET, (verbose_ > 2));
  nllStatus_ = nllLFV_->error();

  standaloneObjectiveFunctionAdapterVEGAS_LFV_ = new svFitStandalone::ObjectiveFunctionAdapterVEGAS_LFV();
}

SVfitStandaloneAlgorithmLFV::~SVfitStandaloneAlgorithmLFV() 
{
  delete nllLFV_;
  delete mcObjectiveFunctionAdapterLFV_;
  delete mcPtEtaPhiMassAdapterLFV_;
  delete standaloneObjectiveFunctionAdapterVEGAS_LFV_;
}

void
SVfitStandaloneAlgorithmLFV::setup()
{
  using namespace svFitStandalone;

  if ( verbose_ >= 1 ) {
    std::cout << "<SVfitStandaloneAlgorithmLFV::setup()>:" << std::endl;
  }
  if ( verbose_ >= 1 ) {
    std::cout << " --> upper limit of leg1::mNuNu will be set to "; 
    if ( nllLFV_->measuredTauLepton().type() == kHadDecay ) { 
      std::cout << "0";
    } else {
      std::cout << (svFitStandalone::tauLeptonMass - TMath::Min(nllLFV_->measuredTauLepton().mass(), 1.5));
    } 
    std::cout << std::endl;
  }
  // start values for xFrac
  minimizer_->SetLimitedVariable(
    kXFrac, 
    std::string(TString::Format("leg%i::xFrac", 1)).c_str(), 
    0.5, 0.1, 0., 1.);
  // start values for nunuMass
  if ( nllLFV_->measuredTauLepton().type() == kHadDecay ) { 
    minimizer_->SetFixedVariable(
      kMNuNu, 
      std::string(TString::Format("leg%i::mNuNu", 1)).c_str(), 
      0.); 
  } else { 
    minimizer_->SetLimitedVariable(
      kMNuNu, 
      std::string(TString::Format("leg%i::mNuNu", 1)).c_str(), 
      0.8, 0.10, 0., svFitStandalone::tauLeptonMass - TMath::Min(nllLFV_->measuredTauLepton().mass(), 1.5)); 
  }
  // start values for phi
  minimizer_->SetVariable(
    kPhi, 
    std::string(TString::Format("leg%i::phi", 1)).c_str(), 
    0.0, 0.25);
  // start values for Pt and mass of visible tau decay products (hadronic tau decays only)
  if ( nllLFV_->measuredTauLepton().type() == kHadDecay && shiftVisMassAndPt_ ) {
    minimizer_->SetLimitedVariable(
      kVisMassShifted, 
      std::string(TString::Format("leg%i::mVisShift", 1)).c_str(), 
      0.8, 0.10, svFitStandalone::chargedPionMass, svFitStandalone::tauLeptonMass);
    minimizer_->SetLimitedVariable(
      kRecTauPtDivGenTauPt, 
      std::string(TString::Format("leg%i::tauPtDivGenVisPt", 1)).c_str(), 
      0., 0.10, -1., +1.5);
  } else {
    minimizer_->SetFixedVariable(
      kVisMassShifted, 
      std::string(TString::Format("leg%i::mVisShift", 1)).c_str(), 
      nllLFV_->measuredTauLepton().mass());
    minimizer_->SetFixedVariable(
      kRecTauPtDivGenTauPt, 
      std::string(TString::Format("leg%i::tauPtDivGenVisPt", 1)).c_str(), 
      0.);
  }
  // set values corresponding to prompt lepton to constant values
  // (effectively elliminating all parameters corresponding to prompt lepton in the fit)
  minimizer_->SetFixedVariable(
    kMaxFitParams + kXFrac, 
    std::string(TString::Format("leg%i::xFrac", 2)).c_str(), 
    0.); 
  minimizer_->SetFixedVariable(
    kMaxFitParams + kMNuNu, 
    std::string(TString::Format("leg%i::mNuNu", 2)).c_str(), 
    0.); 
  minimizer_->SetFixedVariable(
    kMaxFitParams + kPhi, 
    std::string(TString::Format("leg%i::phi", 2)).c_str(), 
    0.); 
  minimizer_->SetFixedVariable(
    kMaxFitParams + kVisMassShifted, 
    std::string(TString::Format("leg%i::mVisShift", 2)).c_str(), 
    nllLFV_->measuredTauLepton().mass());
  minimizer_->SetFixedVariable(
    kMaxFitParams + kRecTauPtDivGenTauPt, 
    std::string(TString::Format("leg%i::tauPtDivGenVisPt", 2)).c_str(), 
    0.);
}

void
SVfitStandaloneAlgorithmLFV::fit()
{
  if ( verbose_ >= 1 ) {
    std::cout << "<SVfitStandaloneAlgorithmLFV::fit()>" << std::endl
	      << " dimension of fit    : " << svFitStandalone::kMaxFitParams << std::endl
	      << " maxObjFunctionCalls : " << maxObjFunctionCalls_ << std::endl; 
  }
  // clear minimizer
  minimizer_->Clear();
  // set verbosity level of minimizer
  minimizer_->SetPrintLevel(-1);
  // setup the function to be called and the dimension of the fit
  ROOT::Math::Functor toMinimize(standaloneObjectiveFunctionAdapterMINUIT_LFV_, svFitStandalone::kMaxFitParams);
  minimizer_->SetFunction(toMinimize); 
  setup();
  minimizer_->SetMaxFunctionCalls(maxObjFunctionCalls_);
  // set Minuit strategy = 2, in order to get reliable error estimates:
  // http://www-cdf.fnal.gov/physics/statistics/recommendations/minuit.html
  minimizer_->SetStrategy(2);
  // compute uncertainties for increase of objective function by 0.5 wrt. 
  // minimum (objective function is log-likelihood function)
  minimizer_->SetErrorDef(0.5);
  if ( verbose_ >= 1 ) {
    std::cout << "starting ROOT::Math::Minimizer::Minimize..." << std::endl;
    std::cout << " #freeParameters = " << minimizer_->NFree() << ","
  	      << " #constrainedParameters = " << (minimizer_->NDim() - minimizer_->NFree()) << std::endl;
  }
  // do the minimization
  nllLFV_->addDelta(false);
  nllLFV_->addSinTheta(true);
  minimizer_->Minimize();
  if ( verbose_ >= 2 ) { 
    minimizer_->PrintResults(); 
  };
  /* get Minimizer status code, check if solution is valid:
    
     0: Valid solution
     1: Covariance matrix was made positive definite
     2: Hesse matrix is invalid
     3: Estimated distance to minimum (EDM) is above maximum
     4: Reached maximum number of function calls before reaching convergence
     5: Any other failure
  */
  fitStatus_ = minimizer_->Status();
  if ( verbose_ >=1 ) { 
    std::cout << "--> fitStatus = " << fitStatus_ << std::endl; 
  }
  
  // and write out the result
  using svFitStandalone::kXFrac;
  using svFitStandalone::kMNuNu;
  using svFitStandalone::kPhi;
  using svFitStandalone::kMaxFitParams;
  // update di-tau system with final fit results
  nllLFV_->results(fittedTauLeptons_, minimizer_->X());
  // determine uncertainty of the fitted di-tau mass
  double x1RelErr = minimizer_->Errors()[kXFrac]/minimizer_->X()[kXFrac];
  // this gives a unified treatment for retrieving the result for integration mode and fit mode
  fittedDiTauSystem_ = fittedTauLeptons_[0] + fittedTauLeptons_[1];
  mass_ = fittedDiTauSystem().mass();
  massUncert_ = 0.50*x1RelErr*fittedDiTauSystem().mass();
  if ( verbose_ >= 2 ) {
    std::cout << ">> -------------------------------------------------------------" << std::endl;
    std::cout << ">> Resonance Record: " << std::endl;
    std::cout << ">> -------------------------------------------------------------" << std::endl;
    std::cout << ">> pt  (di-tau)    = " << fittedDiTauSystem().pt  () << std::endl;
    std::cout << ">> eta (di-tau)    = " << fittedDiTauSystem().eta () << std::endl;
    std::cout << ">> phi (di-tau)    = " << fittedDiTauSystem().phi () << std::endl;
    std::cout << ">> mass(di-tau)    = " << fittedDiTauSystem().mass() << std::endl;  
    std::cout << ">> massUncert      = " << massUncert_ << std::endl
	      << "   error[xFrac1]   = " << minimizer_->Errors()[kXFrac] << std::endl
	      << "   value[xFrac1]   = " << minimizer_->X()[kXFrac]      << std::endl
	      << "   error[xFrac2]   = " << minimizer_->Errors()[kMaxFitParams+kXFrac] << std::endl
	      << "   value[xFrac2]   = " << minimizer_->X()[kMaxFitParams+kXFrac]      << std::endl;
    for ( size_t leg = 0; leg < 2 ; ++leg ) {
      std::cout << ">> -------------------------------------------------------------" << std::endl;
      std::cout << ">> Leg " << leg+1 << " Record: " << std::endl;
      std::cout << ">> -------------------------------------------------------------" << std::endl;
      std::cout << ">> pt  (meas)      = " << nllLFV_->measuredTauLeptons()[leg].p4().pt () << std::endl;
      std::cout << ">> eta (meas)      = " << nllLFV_->measuredTauLeptons()[leg].p4().eta() << std::endl;
      std::cout << ">> phi (meas)      = " << nllLFV_->measuredTauLeptons()[leg].p4().phi() << std::endl; 
      std::cout << ">> pt  (fit )      = " << fittedTauLeptons()[leg].pt () << std::endl;
      std::cout << ">> eta (fit )      = " << fittedTauLeptons()[leg].eta() << std::endl;
      std::cout << ">> phi (fit )      = " << fittedTauLeptons()[leg].phi() << std::endl; 
    }
  }
}

void
SVfitStandaloneAlgorithmLFV::integrateVEGAS(const std::string& likelihoodFileName)
{
  using namespace svFitStandalone;
  
  if ( verbose_ >= 1 ) {
    std::cout << "<SVfitStandaloneAlgorithmLFV::integrateVEGAS>:" << std::endl;
    clock_->Start("<SVfitStandaloneAlgorithmLFV::integrateVEGAS>");
  }

  // number of parameters for fit
  int nDim = 0;
  isLep_ = false;
  const TH1* lutVisMassRes = 0;
  const TH1* lutVisPtRes = 0;
  for ( size_t idx = 0; idx < nll_->measuredTauLeptons().size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = nll_->measuredTauLeptons()[idx];
    if ( measuredTauLepton.type() == kHadDecay ) { 
      if ( measuredTauLepton.decayMode() == 0 ) {
	lutVisMassRes = lutVisMassResDM0_;
	lutVisPtRes = lutVisPtResDM0_;
      } else if ( measuredTauLepton.decayMode() == 1 ) {
	lutVisMassRes = lutVisMassResDM1_;
	lutVisPtRes = lutVisPtResDM1_;
      } else if ( measuredTauLepton.decayMode() == 10 ) {
	lutVisMassRes = lutVisMassResDM10_;
	lutVisPtRes = lutVisPtResDM10_;
      } 
      if ( shiftVisMassAndPt_ ) nDim = 4;
      else nDim = 2;
    } else if ( measuredTauLepton.type() == kLepDecay ) {
      isLep_ = true;
      nDim = 3;
    }
  }
  nDim -= 1; // xFrac for second tau is fixed by delta function for test mass

  /* --------------------------------------------------------------------------------------
     lower and upper bounds for integration. Boundaries are defined for each decay channel
     separately. The order is: 
     
     - hadronic {phihad1, (masshad1, pthad1)}
     - leptonic {nunuMass1, philep1}
     
     x0* defines the start value for the integration, xl* defines the lower integation bound, 
     xh* defines the upper integration bound in the following definitions. 
  */
  double* x0 = new double[nDim];
  double* xl = new double[nDim];
  double* xh = new double[nDim];
  if ( isLep_ ) {
    x0[0] = 0.8; 
    xl[0] = 0.0; 
    xh[0] = svFitStandalone::tauLeptonMass; 
    x0[1] = 0.0; 
    xl[1] = -TMath::Pi(); 
    xh[1] = +TMath::Pi();
  } else {
    x0[0] = 0.0; 
    xl[0] = -TMath::Pi(); 
    xh[0] = +TMath::Pi();
    if ( shiftVisMassAndPt_ ) {
      x0[1] = 0.8; 
      xl[1] = svFitStandalone::chargedPionMass; 
      xh[1] = svFitStandalone::tauLeptonMass; 
      x0[2] = 0.0; 
      xl[2] = -1.0; 
      xh[2] = +1.5;
    }
  }

  TH1* histogramMass = makeHistogram("SVfitStandaloneAlgorithm_histogramMass", measuredDiTauSystem().mass()/1.0125, 1.e+4, 1.025);
  TH1* histogramMass_density = (TH1*)histogramMass->Clone(Form("%s_density", histogramMass->GetName()));

  std::vector<double> xGraph;
  std::vector<double> xErrGraph;
  std::vector<double> yGraph;
  std::vector<double> yErrGraph;

  // integrator instance
  ROOT::Math::GSLMCIntegrator ig2("vegas", 0., 1.e-6, 10000);
  //ROOT::Math::GSLMCIntegrator ig2("vegas", 0., 1.e-6, 2000);
  ROOT::Math::Functor toIntegrate(standaloneObjectiveFunctionAdapterVEGAS_LFV_, &ObjectiveFunctionAdapterVEGAS_LFV::Eval, nDim); 
  standaloneObjectiveFunctionAdapterVEGAS_LFV_->SetIsLep(isLep_);
  standaloneObjectiveFunctionAdapterVEGAS_LFV_->SetShiftVisMassAndPt(shiftVisMassAndPt_);
  ig2.SetFunction(toIntegrate);
  nllLFV_->addDelta(true);
  nllLFV_->addSinTheta(false);
  nllLFV_->addPhiPenalty(false);
  nllLFV_->shiftVisMassAndPt(shiftVisMassAndPt_, lutVisMassRes, lutVisPtRes);
  int count = 0;
  double pMax = 0.;
  double mvis = measuredDiTauSystem().mass();
  double mtest = mvis*1.0125;
  bool skiphighmasstail = false;
  for ( int i = 0; i < 100 && (!skiphighmasstail); ++i ) {
    standaloneObjectiveFunctionAdapterVEGAS_LFV_->SetMvis(mvis);
    standaloneObjectiveFunctionAdapterVEGAS_LFV_->SetMtest(mtest);
    double p = ig2.Integral(xl, xh);
    double pErr = ig2.Error();
    if ( verbose_ >= 2 ) {
      std::cout << "--> scan idx = " << i << ": mtest = " << mtest << ", p = " << p << " +/- " << pErr << " (pMax = " << pMax << ")" << std::endl;
    }    
    if ( p > pMax ) {
      mass_ = mtest;
      pMax  = p;
      count = 0;
    } else {
      if ( p < (1.e-3*pMax) ) {
	++count;
	if ( count>= 5 ) {
	  skiphighmasstail = true;
	}
      } else {
	count = 0;
      }
    }
    double mtest_step = TMath::Max(2.5, 0.025*mtest);
    int bin = histogramMass->FindBin(mtest);
    histogramMass->SetBinContent(bin, p*mtest_step);
    histogramMass->SetBinError(bin, pErr*mtest_step);
    xGraph.push_back(mtest);
    xErrGraph.push_back(0.5*mtest_step);
    yGraph.push_back(p);
    yErrGraph.push_back(pErr);
    mtest += mtest_step;
  }
  //mass_ = extractValue(histogramMass, histogramMass_density);
  massUncert_ = extractUncertainty(histogramMass, histogramMass_density);
  if ( verbose_ >= 1 ) {
    std::cout << "--> mass  = " << mass_  << " +/- " << massUncert_ << std::endl;
    std::cout << "   (pMax = " << pMax << ", count = " << count << ")" << std::endl;
  }
  delete histogramMass;
  delete histogramMass_density;
  if ( likelihoodFileName != "" ) {
    size_t numPoints = xGraph.size();
    TGraphErrors* likelihoodGraph = new TGraphErrors(numPoints);
    likelihoodGraph->SetName("svFitLikelihoodGraph");
    for ( size_t iPoint = 0; iPoint < numPoints; ++iPoint ) {
      likelihoodGraph->SetPoint(iPoint, xGraph[iPoint], yGraph[iPoint]);
      likelihoodGraph->SetPointError(iPoint, xErrGraph[iPoint], yErrGraph[iPoint]);
    }
    TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
    likelihoodGraph->Write();
    delete likelihoodFile;
    delete likelihoodGraph;
  }

  delete x0;
  delete xl;
  delete xh;

  if ( verbose_ >= 1 ) {
    clock_->Show("<SVfitStandaloneAlgorithmLFV::integrateVEGAS>");
  }
}

void
SVfitStandaloneAlgorithmLFV::integrateMarkovChain()
{
  using namespace svFitStandalone;
  
  if ( verbose_ >= 1 ) {
    std::cout << "<SVfitStandaloneAlgorithmLFV::integrateMarkovChain>:" << std::endl;
    clock_->Start("<SVfitStandaloneAlgorithmLFV::integrateMarkovChain>");
  }
  if ( isInitialized2_ ) {
    mcPtEtaPhiMassAdapterLFV_->Reset();
  } else {
    // initialize    
    std::string initMode = "none";
    unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls2_);
    unsigned numIterSampling = maxObjFunctionCalls2_;
    unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.02*maxObjFunctionCalls2_);
    unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.06*maxObjFunctionCalls2_);
    double T0 = 15.;
    double alpha = 1.0 - 1.e+2/maxObjFunctionCalls2_;
    unsigned numChains = 1;
    unsigned numBatches = 1;
    unsigned L = 1;
    double epsilon0 = 1.e-2;
    double nu = 0.71;
    int verbose = -1;
    integrator2_ = new SVfitStandaloneMarkovChainIntegrator(
                         initMode, numIterBurnin, numIterSampling, numIterSimAnnealingPhase1, numIterSimAnnealingPhase2,
			 T0, alpha, numChains, numBatches, L, epsilon0, nu,
			 verbose);
    mcObjectiveFunctionAdapterLFV_ = new MCObjectiveFunctionAdapterLFV();
    integrator2_->setIntegrand(*mcObjectiveFunctionAdapterLFV_);
    integrator2_nDim_ = 0;
    mcPtEtaPhiMassAdapterLFV_ = new MCPtEtaPhiMassAdapterLFV();
    integrator2_->registerCallBackFunction(*mcPtEtaPhiMassAdapterLFV_);
    isInitialized2_ = true;    
  }

  // number of parameters for fit
  int nDim = 0;
  isLep_ = false;
  const TH1* lutVisMassRes = 0;
  const TH1* lutVisPtRes = 0;
  for ( size_t idx = 0; idx < nll_->measuredTauLeptons().size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = nll_->measuredTauLeptons()[idx];
    if ( measuredTauLepton.type() == kHadDecay ) { 
      if ( measuredTauLepton.decayMode() == 0 ) {
	lutVisMassRes = lutVisMassResDM0_;
	lutVisPtRes = lutVisPtResDM0_;
      } else if ( measuredTauLepton.decayMode() == 1 ) {
	lutVisMassRes = lutVisMassResDM1_;
	lutVisPtRes = lutVisPtResDM1_;
      } else if ( measuredTauLepton.decayMode() == 10 ) {
	lutVisMassRes = lutVisMassResDM10_;
	lutVisPtRes = lutVisPtResDM10_;
      }	
      if ( shiftVisMassAndPt_ ) nDim = 4;
      else nDim = 2;
    } else if ( measuredTauLepton.type() == kLepDecay ) {
      isLep_ = true;
      nDim = 3;
    }
  }
  
  if ( nDim != integrator2_nDim_ ) {
    mcObjectiveFunctionAdapterLFV_->SetNDim(nDim);    
    integrator2_->setIntegrand(*mcObjectiveFunctionAdapterLFV_);
    mcPtEtaPhiMassAdapterLFV_->SetNDim(nDim);
    integrator2_nDim_ = nDim;
  }
  mcObjectiveFunctionAdapterLFV_->SetIsLep(isLep_);
  mcObjectiveFunctionAdapterLFV_->SetShiftVisMassAndPt(shiftVisMassAndPt_);
  mcPtEtaPhiMassAdapterLFV_->SetIsLep(isLep_);
  mcPtEtaPhiMassAdapterLFV_->SetShiftVisMassAndPt(shiftVisMassAndPt_);
  
  /* --------------------------------------------------------------------------------------
     lower and upper bounds for integration. Boundaries are defined for each decay channel
     separately. The order is: 
     
     - hadronic {xhad, phihad1, (masshad1, pthad1)}
     - leptonic {xlep, nunuMass1, philep1}
     
     x0* defines the start value for the integration, xl* defines the lower integation bound, 
     xh* defines the upper integration bound in the following definitions. 
  */
  std::vector<double> x0(nDim);
  std::vector<double> xl(nDim);
  std::vector<double> xh(nDim);
  if ( isLep_ ) {
    x0[0] = 0.5; 
    xl[0] = 0.0; 
    xh[0] = 1.0; 
    x0[1] = 0.8; 
    xl[1] = 0.0; 
    xh[1] = svFitStandalone::tauLeptonMass; 
    x0[2] = 0.0; 
    xl[2] = -TMath::Pi(); 
    xh[2] = +TMath::Pi();
  } else {
    x0[0] = 0.5; 
    xl[0] = 0.0; 
    xh[0] = 1.0; 
    x0[1] = 0.0; 
    xl[1] = -TMath::Pi(); 
    xh[1] = +TMath::Pi();
    if ( shiftVisMassAndPt_ ) {
      x0[2] = 0.8; 
      xl[2] = svFitStandalone::chargedPionMass; 
      xh[2] = svFitStandalone::tauLeptonMass; 
      x0[3] = 0.0; 
      xl[3] = -1.0; 
      xh[3] = +1.5;
    }
  }
  for ( int i = 0; i < nDim; ++i ) {
    // transform startPosition into interval ]0..1[
    // expected by MarkovChainIntegrator class
    x0[i] = (x0[i] - xl[i])/(xh[i] - xl[i]);
    //std::cout << "x0[" << i << "] = " << x0[i] << std::endl;
  }
  integrator2_->initializeStartPosition_and_Momentum(x0);
  nllLFV_->addDelta(false);
  nllLFV_->addSinTheta(false);
  nllLFV_->addPhiPenalty(false);
  nllLFV_->shiftVisMassAndPt(shiftVisMassAndPt_, lutVisMassRes, lutVisPtRes);
  double integral = 0.;
  double integralErr = 0.;
  int errorFlag = 0;
  integrator2_->integrate(xl, xh, integral, integralErr, errorFlag);
  fitStatus_ = errorFlag;
  pt_ = mcPtEtaPhiMassAdapterLFV_->getPt();
  ptUncert_ = mcPtEtaPhiMassAdapterLFV_->getPtUncert();
  eta_ = mcPtEtaPhiMassAdapterLFV_->getEta();
  etaUncert_ = mcPtEtaPhiMassAdapterLFV_->getEtaUncert();
  phi_ = mcPtEtaPhiMassAdapterLFV_->getPhi();
  phiUncert_ = mcPtEtaPhiMassAdapterLFV_->getPhiUncert();
  mass_ = mcPtEtaPhiMassAdapterLFV_->getMass();
  massUncert_ = mcPtEtaPhiMassAdapterLFV_->getMassUncert();
  fittedDiTauSystem_ = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(pt_, eta_, phi_, mass_);
  if ( verbose_ >= 1 ) {
    std::cout << "--> Pt = " << pt_ << ", eta = " << eta_ << ", phi = " << phi_ << ", mass  = " << mass_ << std::endl;
    clock_->Show("<SVfitStandaloneAlgorithmLFV::integrateMarkovChain>");
  }
}

