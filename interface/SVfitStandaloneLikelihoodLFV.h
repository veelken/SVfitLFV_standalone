#ifndef TauAnalysis_SVfitStandalone_SVfitStandaloneLikelihoodLFV_h
#define TauAnalysis_SVfitStandalone_SVfitStandaloneLikelihoodLFV_h

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"
#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "Math/Minimizer.h"

#include <vector>

/**
   \class   SVfitStandaloneLikelihoodLFV SVfitStandaloneLikelihoodLFV.h "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"
     
   \brief   Negative log likelihood for a resonance X to decay via the lepton flavor violating decay modes X -> mu tau -> mu tau_had, X -> mu tau -> mu tau_e
            and X -> e tau -> e tau_had, X -> e tau -> e tau_mu.
     
   Negative log likelihood for a resonance decay into two tau leptons that may themselves decay hadronically or leptonically 
   Depending on the configuration during object creation it will be a combination of MET, TauToHad, TauToLep and additional
   penalty terms, e.g. to suppress tails in m(tau,tau) (logM). Configurables during creation time are:
     
   \var measuredTauLeptons : the vector of the reconstructed tau decay products and of the prompt lepton
   \var measuredMET        : the spacial vector of the measured MET
   \var covMET             : the covariance matrix of the MET (as determined from the MEt significance for instance)
   \verbose                : indicating the verbosity level 

   In fit mode additional optional values may be set before the fit is performed. During construction the class is initialized with 
   default values as indicated in braces (below):

   \var metPower : indicating an additional power to enhance the MET likelihood (default is 1.)
   \var addLogM  : specifying whether to use the LogM penalty term or not (default is true)     

   A typical way to obtain the covariance matrix of the MET is to follow the MET significance algorithm as provided by RecoMET.
   The SVfitStandaloneLikelihood class is for internal use only. The general use calse is to access it from the class 
   SVfitStandaloneAlgorithm as defined in interface/SVfitStandaloneAlgorithm.h in the same package. The SVfitLikelihood class 
   keeps all necessary information to calculate the combined likelihood but does not perform any fit nor integration. It is 
   interfaced to the ROOT minuit minimization package or to the VEGAS integration packages via the global function pointer 
   gSVfitStandaloneLikelihood as defined in src/SVfitStandaloneLikelihood.cc in the same package. 
*/
 
namespace svFitStandalone
{
  class SVfitStandaloneLikelihoodLFV : public SVfitStandaloneLikelihood
  {
   public:
    /// constructor with a minimla set of configurables 
    SVfitStandaloneLikelihoodLFV(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, bool verbose);
    /// default destructor
    ~SVfitStandaloneLikelihoodLFV() {};
    /// static pointer to this (needed for the minuit function calls)
    static const SVfitStandaloneLikelihoodLFV* gSVfitStandaloneLikelihoodLFV;
       
    /// marginalize unknown mass of hadronic tau decay products (ATLAS case)
    void marginalizeVisMass(bool value, const TH1* l1lutVisMass);  

    /// take resolution on energy and mass of hadronic tau decays into account
    void shiftVisMass(bool value, const TH1* lutVisMassRes);    
    void shiftVisPt(bool value, const TH1* lutVisPtRes);    

    /// fit function to be called from outside. Has to be const to be usable by minuit. This function will call the actual 
    /// functions transform and prob internally 
    double prob(const double* x, bool fixToMtest = false, double mtest = -1.) const;

    /// return vector of fitted tau leptons, which will be the actual fit result. This function is a subset of transform.
    /// It needs to be factored out though as transform has to be const to be usable by minuit and therefore is not allowed 
    /// change the class members.  
    void results(std::vector<LorentzVector>& fittedTauLeptons, const double* x) const;

    /// return measured tau lepton
    const MeasuredTauLepton& measuredTauLepton() const { return measuredTauLepton_; }
    /// return prompt lepton
    const MeasuredTauLepton& promptLepton() const { return promptLepton_; }

   private:
    /// transformation from x to xPrime, x are the actual fit parameters, xPrime are the transformed parameters that go into 
    /// the prob function. Has to be const to be usable by minuit.
    const double* transform(double* xPrime, const double* x, bool fixToMtest, double mtest) const;
    /// combined likelihood function. The same function os called for fit and integratino mode. Has to be const to be usable 
    /// by minuit/VEGAS/MarkovChain. The additional boolean phiPenalty is added to prevent singularities at the +/-pi boundaries 
    /// of kPhi within the fit parameters (kFitParams). It is only used in fit mode. In integration mode the passed on value 
    /// is always 0. 
    double prob(const double* xPrime, double phiPenalty) const;
    
   private:
    /// prompt lepton
    MeasuredTauLepton promptLepton_;
    /// measured tau lepton
    MeasuredTauLepton measuredTauLepton_;

    /// resolution on energy and mass of hadronic taus
    const TH1* lutVisMass_;  
    const TH1* lutVisMassRes_;
    const TH1* lutVisPtRes_;
  };
}

#endif
