#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneAlgorithmLFV.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

#include <TGraphErrors.h>

namespace svFitStandalone
{
  void map_xLFV(const double* x, int nDim, double* x_mapped)
  {
    if ( nDim == 2 ) {
      x_mapped[kXFrac] = x[0];
      x_mapped[kMNuNu] = 0.;
      x_mapped[kPhi]   = x[1];
    } else if ( nDim == 3 ) {
      x_mapped[kXFrac] = x[0];
      x_mapped[kMNuNu] = x[1];
      x_mapped[kPhi]   = x[2];
    } else assert(0);
    //std::cout << "<map_x>:" << std::endl;
    //for ( int i = 0; i < 3; ++i ) {
    //  std::cout << " x_mapped[" << i << "] = " << x_mapped[i] << std::endl;
    //}
  }
}

SVfitStandaloneAlgorithmLFV::SVfitStandaloneAlgorithmLFV(const std::vector<svFitStandalone::MeasuredTauLepton>& measuredTauLeptons, const svFitStandalone::Vector& measuredMET, const TMatrixD& covMET, 
							 unsigned int verbose) 
  : SVfitStandaloneAlgorithm(measuredTauLeptons, measuredMET, covMET, verbose),
    mcObjectiveFunctionAdapterLFV_(0),
    mcPtEtaPhiMassAdapterLFV_(0)
{ 
  // instantiate the combined likelihood
  nllLFV_ = new svFitStandalone::SVfitStandaloneLikelihoodLFV(measuredTauLeptons, measuredMET, covMET, (verbose_ > 2));
  nllStatus_ = nllLFV_->error();
}

SVfitStandaloneAlgorithmLFV::~SVfitStandaloneAlgorithmLFV() 
{
  delete nllLFV_;
  delete mcObjectiveFunctionAdapterLFV_;
  delete mcPtEtaPhiMassAdapterLFV_;
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
    if ( nllLFV_->measuredTauLepton().decayType() == kHadDecay ) { 
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
  if ( nllLFV_->measuredTauLepton().decayType() == kHadDecay ) { 
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
  ROOT::Math::Functor toMinimize(standaloneObjectiveFunctionAdapterLFV_, svFitStandalone::kMaxFitParams);
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

  double pi = 3.14159265;
  // number of hadrponic decays
  int khad = 0;
  for ( size_t idx = 0; idx < nllLFV_->measuredTauLeptons().size(); ++idx ) {
    if ( nllLFV_->measuredTauLeptons()[idx].decayType() == kHadDecay ) { 
      khad++; 
    }
  }
  // number of parameters for fit
  int par = svFitStandalone::kMaxFitParams - (khad + 1);
  if ( verbose_ >= 1 ) {
    std::cout << "par = " << par << std::endl;
  }
  /* --------------------------------------------------------------------------------------
     lower and upper bounds for integration. Boundaries are defined for each decay channel
     separately. The order is: 
     
     - 1dim : hadronic {phihad1}
     - 2dim : leptonic {nunuMass, philep}
     
     xl* defines the lower integation bound, xu* defines the upper integration bound in 
     the following definitions. 
  */
  double xl1[1] = { -pi };
  double xu1[1] = {  pi };
  double xl2[2] = { 0.0, -pi };
  double xu2[2] = { svFitStandalone::tauLeptonMass, pi };

  TH1* histogramMass = makeHistogram("SVfitStandaloneAlgorithmLFV_histogramMass", measuredDiTauSystem().mass()/1.0125, 1.e+4, 1.025);
  TH1* histogramMass_density = (TH1*)histogramMass->Clone(Form("%s_density", histogramMass->GetName()));

  std::vector<double> xGraph;
  std::vector<double> xErrGraph;
  std::vector<double> yGraph;
  std::vector<double> yErrGraph;

  // integrator instance
  //ROOT::Math::GSLMCIntegrator ig2("vegas", 0., 1.e-6, 2000);
  ROOT::Math::GSLMCIntegrator ig2("vegas", 0., 1.e-6, 10000);
  ROOT::Math::Functor toIntegrate(&standaloneObjectiveFunctionAdapterLFV_, &ObjectiveFunctionAdapterLFV::Eval, par); 
  standaloneObjectiveFunctionAdapterLFV_.SetPar(par);
  ig2.SetFunction(toIntegrate);
  nllLFV_->addDelta(true);
  nllLFV_->addSinTheta(false);
  nllLFV_->addPhiPenalty(false);
  int count = 0;
  double pMax = 0.;
  double mtest = measuredDiTauSystem().mass();
  bool skiphighmasstail = false;
  for ( int i = 0; i < 100 && (!skiphighmasstail); ++i ) {
    standaloneObjectiveFunctionAdapterLFV_.SetM(mtest);
    double p = -1.;
    if ( par == 1 ) {
      p = ig2.Integral(xl1, xu1);
    } else if ( par == 2 ) {
      p = ig2.Integral(xl2, xu2);
    } else {
      std::cout << " >> ERROR : there must be one measured lepton + one prompt lepton" << std::endl;
      assert(0);
    }
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
    mcPtEtaPhiMassAdapter_ = new MCPtEtaPhiMassAdapterLFV();
    integrator2_->registerCallBackFunction(*mcPtEtaPhiMassAdapterLFV_);
    isInitialized2_ = true;    
  }

  const double pi = TMath::Pi();
  // number of hadronic decays
  int khad = 0;
  for ( size_t idx = 0; idx < nllLFV_->measuredTauLeptons().size(); ++idx ) {
    if ( nllLFV_->measuredTauLeptons()[idx].decayType() == kHadDecay ) { 
      ++khad; 
    }
  }
  // number of parameters for fit
  int nDim = svFitStandalone::kMaxFitParams - khad;  
  if ( nDim != integrator2_nDim_ ) {
    mcObjectiveFunctionAdapterLFV_->SetNDim(nDim);
    integrator2_->setIntegrand(*mcObjectiveFunctionAdapterLFV_);
    mcPtEtaPhiMassAdapterLFV_->SetNDim(nDim);
    integrator2_nDim_ = nDim;
  }
  /* --------------------------------------------------------------------------------------
     lower and upper bounds for integration. Boundaries are defined for each decay channel
     separately. The order is: 
     
     - 2dim : hadronic {xhad, phihad}
     - 3dim : leptonic {xlep, nunuMass, philep}
     
     x0* defines the start value for the integration, xl* defines the lower integation bound, 
     xu* defines the upper integration bound in the following definitions. 
     ATTENTION: order matters here! In the semi-leptonic decay the lepton must go first in 
     the parametrization, as it is first in the definition of integral boundaries. This is
     the reason why the measuredLeptons are eventually re-ordered in the constructor of 
     this class before passing them on to SVfitStandaloneLikelihood.
  */
  double x02[2] = { 0.5, 0.0 };
  double xl2[2] = { 0.0, -pi };
  double xu2[2] = { 1.0,  pi };
  double x03[3] = { 0.5, 0.8, 0.0 };
  double xl3[3] = { 0.0, 0.0, -pi };
  double xu3[3] = { 1.0, svFitStandalone::tauLeptonMass, pi };
  xu3[1] = svFitStandalone::tauLeptonMass - TMath::Min(nllLFV_->measuredTauLepton().mass(), 1.6);
  x03[1] = 0.5*(xl3[1] + xu3[1]);
  std::vector<double> x0(nDim);
  std::vector<double> xl(nDim);
  std::vector<double> xu(nDim);
  for ( int i = 0; i < nDim; ++i ) {
    if ( nDim == 2 ) {
      x0[i] = x02[i];
      xl[i] = xl2[i];
      xu[i] = xu2[i];
    } else if ( nDim == 3 ) {
      x0[i] = x03[i];
      xl[i] = xl3[i];
      xu[i] = xu3[i];
    } else {
      std::cerr << "<SVfitStandaloneAlgorithmLFV>:"
		<< "Exactly one measured lepton + one prompt lepton required --> ABORTING !!\n";
      assert(0);
    }
    // transform startPosition into interval ]0..1[
    // expected by MarkovChainIntegrator class
    x0[i] = (x0[i] - xl[i])/(xu[i] - xl[i]);
    //std::cout << "x0[" << i << "] = " << x0[i] << std::endl;
  }
  integrator2_->initializeStartPosition_and_Momentum(x0);
  nllLFV_->addDelta(false);
  nllLFV_->addSinTheta(false);
  nllLFV_->addPhiPenalty(false);
  double integral = 0.;
  double integralErr = 0.;
  int errorFlag = 0;
  integrator2_->integrate(xl, xu, integral, integralErr, errorFlag);
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
