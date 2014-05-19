#ifndef TauAnalysis_SVfitStandalone_SVfitStandaloneAlgorithmLFV_h
#define TauAnalysis_SVfitStandalone_SVfitStandaloneAlgorithmLFV_h

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneLikelihoodLFV.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneMarkovChainIntegrator.h"
#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"

#include <TMath.h>
#include <TArrayF.h>
#include <TString.h>
#include <TH1.h>
#include <TBenchmark.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

using svFitStandalone::Vector;
using svFitStandalone::LorentzVector;
using svFitStandalone::MeasuredTauLepton;

/**
   \class   SVfitStandaloneAlgorithmLFV.h "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithmLFV.h"
   
   \brief   Function interface to minuit.
   
   This class is an interface, which is used as global function pointer of the combined likelihood as defined in src/SVfitStandaloneLikelihood.cc
   to VEGAS or minuit. It is a member of the of the SVfitStandaloneAlgorithm class defined below and is used in SVfitStandalone::fit(), or 
   SVfitStandalone::integrate(), where it is passed on to a ROOT::Math::Functor. The parameters x correspond to the array of fit/integration 
   paramters as defined in interface/SVfitStandaloneLikelihood.h of this package. In the fit mode these are made known to minuit in the function
   SVfitStandaloneAlgorithm::setup. In the integration mode the mapping is done internally in the SVfitStandaloneLikelihood::tansformint. This
   has to be in sync. with the definition of the integration boundaries in SVfitStandaloneAlgorithm::integrate. 
*/

namespace svFitStandalone
{
  class ObjectiveFunctionAdapterMINUIT_LFV
  {
  public:
    // for minuit fit
    double operator()(const double* x) const // NOTE: return value = -log(likelihood)
    {
      double prob = SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV->prob(x);
      double nll;
      if ( prob > 0. ) nll = -TMath::Log(prob);
      else nll = std::numeric_limits<float>::max();
      return nll;
    }
  };
  // for VEGAS integration
  void map_xVEGAS_LFV(const double*, bool, bool, double, double, double*);
  class ObjectiveFunctionAdapterVEGAS_LFV
  {
  public:
    double Eval(const double* x) const // NOTE: return value = likelihood, **not** -log(likelihood)
    {
      map_xVEGAS_LFV(x, isLep_, shiftVisMassAndPt_, mvis_, mtest_, x_mapped_);  
      double prob = SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV->prob(x_mapped_, true, mtest_);      
      if ( TMath::IsNaN(prob) ) prob = 0.;
      return prob;
    }
    void SetIsLep(bool isLep) { isLep_ = isLep; }
    void SetShiftVisMassAndPt(bool shiftVisMassAndPt) { shiftVisMassAndPt_ = shiftVisMassAndPt; }
    void SetMvis(double mvis) { mvis_ = mvis; }
    void SetMtest(double mtest) { mtest_ = mtest; }
  private:
    mutable double x_mapped_[5];
    bool isLep_;
    bool shiftVisMassAndPt_;
    double mvis_;  // mass of visible tau decay products
    double mtest_; // current mass hypothesis
  };
  // for markov chain integration
  void map_xMarkovChain_LFV(const double*, bool, bool, double*);
  class MCObjectiveFunctionAdapterLFV : public ROOT::Math::Functor
  {
   public:
    void SetIsLep(bool isLep) { isLep_ = isLep; }
    void SetShiftVisMassAndPt(bool shiftVisMassAndPt) { shiftVisMassAndPt_ = shiftVisMassAndPt; }
    void SetNDim(int nDim) { nDim_ = nDim; }
    unsigned int NDim() const { return nDim_; }
   private:
    virtual double DoEval(const double* x) const
    {
      map_xMarkovChain_LFV(x, isLep_, shiftVisMassAndPt_, x_mapped_);
      double prob = SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV->prob(x_mapped_);
      if ( TMath::IsNaN(prob) ) prob = 0.;
      return prob;
    } 
    mutable double x_mapped_[5];
    int nDim_;
    bool isLep_;
    bool shiftVisMassAndPt_;
  };
  class MCPtEtaPhiMassAdapterLFV : public MCPtEtaPhiMassAdapter
  {
   public:
    MCPtEtaPhiMassAdapterLFV() 
      : MCPtEtaPhiMassAdapter()
    {}
    ~MCPtEtaPhiMassAdapterLFV() {}
    void SetIsLep(bool isLep) { isLep_ = isLep; }
   private:    
    virtual double DoEval(const double* x) const
    {
      map_xMarkovChain_LFV(x, isLep_, shiftVisMassAndPt_, x_mappedLFV_);
      SVfitStandaloneLikelihoodLFV::gSVfitStandaloneLikelihoodLFV->results(fittedTauLeptons_, x_mappedLFV_);
      fittedDiTauSystem_ = fittedTauLeptons_[0] + fittedTauLeptons_[1];
      //std::cout << "<MCPtEtaPhiMassAdapterLFV::DoEval>" << std::endl;
      //std::cout << " Pt = " << fittedDiTauSystem_.pt() << "," 
      //	  << " eta = " << fittedDiTauSystem_.eta() << "," 
      //	  << " phi = " << fittedDiTauSystem_.phi() << ","
      //	  << " mass = " << fittedDiTauSystem_.mass() << std::endl;
      histogramPt_->Fill(fittedDiTauSystem_.pt());
      histogramEta_->Fill(fittedDiTauSystem_.eta());
      histogramPhi_->Fill(fittedDiTauSystem_.phi());
      histogramMass_->Fill(fittedDiTauSystem_.mass());
      return 0.;
    } 
    mutable double x_mappedLFV_[5];
    bool isLep_;
  };
}

/**
   \class   SVfitStandaloneAlgorithmLFV SVfitStandaloneAlgorithmLFV.h "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
     
   \brief   Standalone version of the SVfitAlgorithm.
            Special version for lepton flavor violating Higgs decays H -> mu tau -> mu tau_had, H -> mu tau -> mu tau_e
            and H -> e tau -> e tau_had, H -> e tau -> e tau_mu.

   This class is a standalone version of the SVfitAlgorithm to perform the full reconstruction of a di-tau resonance system. The 
   implementation is supposed to deal with any combination of leptonic or hadronic tau decays. It exploits likelihood functions 
   as defined in interface/LikelihoodFunctions.h of this package, which are combined into a single likelihood function as defined 
   interface/SVfitStandaloneLikelihood.h in this package. The combined likelihood function depends on the following variables: 

   \var nunuMass   : the invariant mass of the neutrino system for each decay branch (two parameters)
   \var decayAngle : the decay angle in the restframe of each decay branch (two parameters)
   \var visMass    : the mass of the visible component of the di-tau system (two parameters)

   The actual fit parameters are:

   \var nunuMass   : the invariant mass of the neutrino system for each decay branch (two parameters)
   \var xFrac      : the fraction of the visible energy on the energy of the tau lepton in the labframe (two parameters)
   \var phi        : the azimuthal angle of the tau lepton (two parameters)

   In the fit mode. The azimuthal angle of each tau lepton is not constraint by measurement. It is limited to the physical values 
   from -Math::Phi to Math::Phi in the likelihood function of the combined likelihood class. The parameter nunuMass is constraint 
   to the tau lepton mass minus the mass of the visible part of the decay (which is itself constraint to values below the tau 
   lepton mass) in the setup function of this class. The parameter xFrac is constraint to values between 0. and 1. in the setup 
   function of this class. The invariant mass of the neutrino system is fixed to be zero for hadronic tau lepton decays as only 
   one (tau-) neutrino is involved in the decay. The original number of free parameters of 3 is therefore reduced by one for each 
   hadronic tau decay within the resonance. All information about the negative log likelihood is stored in the SVfitStandaloneLikelihood 
   class as defined in the same package. This class interfaces the combined likelihood to the ROOT::Math::Minuit minimization program. 
   It does setup/initialize the fit parameters as defined in interface/SVfitStandaloneLikelihood.h in this package, initializes the 
   minimization procedure, executes the fit algorithm and returns the fit result. The fit result consists of the fully reconstructed 
   di-tau system, from which also the invariant mass can be derived.

   In the integration mode xFrac for the second leptons is determiend from xFrac of the first lepton for given di-tau mass, thus reducing 
   the number of parameters to be integrated out wrt. to the fit version by one. The di-tau mass is scanned for the highest likelihood 
   starting from the visible mass of the two leptons. The return value is just the di-tau mass. 

   Common usage is: 
   
   // construct the class object from the minimal necessary information
   SVfitStandaloneAlgorithmLFV algo(measuredTauLeptons, measuredMET, covMET);
   // apply customized configurations if wanted (examples are given below)
   //algo.maxObjFunctionCalls(10000); // only applies for fit mode
   //algo.addLogM(false);             // applies for fit and integration mode
   //algo.metPower(0.5);              // only applies for fit mode
   // run the fit in fit mode
   algo.fit();
   // retrieve the results upon success
   if ( algo.isValidSolution() ) {
     std::cout << algo.mass();
   }
   // run the integration in integration mode
   algo.integrate();
   std::cout << algo.mass();

   The following optional parameters can be applied after initialization but before running the fit in fit mode: 

   \var metPower : indicating an additional power to enhance the MET likelihood (default is 1.)
   \var addLogM : specifying whether to use the LogM penalty term or not (default is true)     
   \var maxObjFunctionCalls : the maximum of function calls before the minimization procedure is terminated (default is 5000)
*/

class SVfitStandaloneAlgorithmLFV : public SVfitStandaloneAlgorithm
{
 public:
  /// constructor from a minimal set of configurables
  SVfitStandaloneAlgorithmLFV(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, unsigned int verbose = 0);
  /// destructor
  ~SVfitStandaloneAlgorithmLFV();

  /// fit to be called from outside
  void fit();
  /// integration by VEGAS (kept for legacy)
  void integrate() { return integrateVEGAS(""); }
  /// integration by VEGAS to be called from outside
  void integrateVEGAS(const std::string& likelihoodFileName = "");
  /// integration by Markov Chain MC to be called from outside
  void integrateMarkovChain();

 private:
  /// setup the starting values for the minimization (default values for the fit parameters are taken from src/SVFitParameters.cc in the same package)
  void setup();

 private:
  /// standalone combined likelihood
  svFitStandalone::SVfitStandaloneLikelihoodLFV* nllLFV_;
  /// needed to make the fit function callable from within minuit
  svFitStandalone::ObjectiveFunctionAdapterMINUIT_LFV standaloneObjectiveFunctionAdapterMINUIT_LFV_;

  /// needed for VEGAS integration
  svFitStandalone::ObjectiveFunctionAdapterVEGAS_LFV* standaloneObjectiveFunctionAdapterVEGAS_LFV_;   

  gsl_monte_function* integrand_;
  gsl_monte_vegas_state* workspace_;
  mutable gsl_rng* rnd_;

  /// needed for markov chain integration
  svFitStandalone::MCObjectiveFunctionAdapterLFV* mcObjectiveFunctionAdapterLFV_;
  svFitStandalone::MCPtEtaPhiMassAdapterLFV* mcPtEtaPhiMassAdapterLFV_;

  bool isLep_;
};

#endif
