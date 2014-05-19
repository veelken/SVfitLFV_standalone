
/**
   \class testSVfitStandaloneLFV testSVfitStandaloneLFV.cc "TauAnalysis/SVfitStandalone/bin/testSVfitStandaloneLFV.cc"
   \brief Basic example of the use of the standalone version of SVfit
          for lepton flavor violating Higgs decays Higgs -> e tau -> e tau_had, Higgs -> e tau -> e tau_mu
	  and Higgs -> mu tau -> mu tau_had, Higgs -> mu tau -> mu tau_e

   This is an example executable to show the use of the standalone version of SVfit 
   from a flat n-tuple or single event.
*/

#include "TauAnalysis/SVfitStandaloneLFV/interface/SVfitStandaloneAlgorithmLFV.h"

#include "TTree.h"
#include "TFile.h"

void singleEvent()
{
  /* 
     This is a single event for testing in the integration mode.
  */
  // define MET
  double metPx =  17.6851;
  double metPy =  23.5161;
  Vector MET(metPx, metPy, 0.);
  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] =  284;
  covMET[1][0] =   13.4;
  covMET[0][1] = covMET[1][0];
  covMET[1][1] =  255.6;
  // define lepton four vectors
  double l1Px   =  -47.5987;
  double l1Py   =  -13.6761;
  double l1Pz   =  -61.7664;
  double l1Mass =    0.105658;
  double l1En   = TMath::Sqrt(l1Px*l1Px + l1Py*l1Py + l1Pz*l1Pz + l1Mass*l1Mass);
  svFitStandalone::LorentzVector l1(l1Px, l1Py, l1Pz, l1En); // "prompt" muon (muon originating from LFV Higgs decay directly)
  double l2Px   =   35.0206;
  double l2Py   =    9.57334;
  double l2Pz   =    9.49413;
  double l2Mass =    1.00231;
  double l2En   = TMath::Sqrt(l2Px*l2Px + l2Py*l2Py + l2Pz*l2Pz + l2Mass*l2Mass);
  svFitStandalone::LorentzVector l2(l2Px, l2Py, l2Pz, l2En); // tau -> had decay
  int l2decayMode = 10;
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kPrompt, l1));
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, l2, l2decayMode));
  // define algorithm (set the debug level to 3 for testing)
  SVfitStandaloneAlgorithmLFV algo(measuredTauLeptons, MET, covMET, 2);
  algo.addLogM(false);
  TFile* inputFile = new TFile("../data/svFitVisMassAndPtResolutionPDF.root");  
  algo.shiftVisMassAndPt(true, inputFile);
  /* 
     the following lines show how to use the different methods on a single event
  */
  // minuit fit method
  //algo.fit();
  // integration by VEGAS (default)
  algo.integrateVEGAS();
  // integration by markov chain MC
  //algo.integrateMarkovChain();

  double mass = algo.getMass(); // return value is in units of GeV
  if ( algo.isValidSolution() ) {
    std::cout << "found mass = " << mass << " (expected value = 115.68)" << std::endl;
  } else {
    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  }
  return;

  delete inputFile;
}

void eventsFromTree(int argc, char* argv[]) 
{
  // parse arguments
  if ( argc < 3 ) {
    std::cout << "Usage : " << argv[0] << " [inputfile.root] [tree_name]" << std::endl;
    return;
  }
  // get intput directory up to one before mass points
  TFile* file = new TFile(argv[1]); 
  // access tree in file
  TTree* tree = (TTree*) file->Get(argv[2]);
  // input variables
  float met, metPhi;
  float covMet11, covMet12; 
  float covMet21, covMet22;
  float l1M, l1Px, l1Py, l1Pz;
  float l2M, l2Px, l2Py, l2Pz;
  float mTrue;
  // branch adresses
   tree->SetBranchAddress("met", &met);
  tree->SetBranchAddress("mphi", &metPhi);
  tree->SetBranchAddress("mcov_11", &covMet11);
  tree->SetBranchAddress("mcov_12", &covMet12);
  tree->SetBranchAddress("mcov_21", &covMet21);
  tree->SetBranchAddress("mcov_22", &covMet22);
  tree->SetBranchAddress("l1_M", &l1M);
  tree->SetBranchAddress("l1_Px", &l1Px);
  tree->SetBranchAddress("l1_Py", &l1Py);
  tree->SetBranchAddress("l1_Pz", &l1Pz);
  tree->SetBranchAddress("l2_M", &l2M);
  tree->SetBranchAddress("l2_Px", &l2Px);
  tree->SetBranchAddress("l2_Py", &l2Py);
  tree->SetBranchAddress("l2_Pz", &l2Pz);
  tree->SetBranchAddress("m_true", &mTrue);
  int nevent = tree->GetEntries();
  for ( int i = 0; i < nevent; ++i ) {
    tree->GetEvent(i);
    std::cout << "event " << (i + 1) << std::endl;
    // setup MET input vector
    svFitStandalone::Vector measuredMET(met *TMath::Sin(metPhi), met *TMath::Cos(metPhi), 0); 
    // setup the MET significance
    TMatrixD covMET(2,2);
    covMET[0][0] = covMet11;
    covMET[0][1] = covMet12;
    covMET[1][0] = covMet21;
    covMET[1][1] = covMet22;
    // setup measure tau lepton vectors 
    svFitStandalone::LorentzVector l1(l1Px, l1Py, l1Pz, TMath::Sqrt(l1M*l1M + l1Px*l1Px + l1Py*l1Py + l1Pz*l1Pz));
    svFitStandalone::LorentzVector l2(l2Px, l2Py, l2Pz, TMath::Sqrt(l2M*l2M + l2Px*l2Px + l2Py*l2Py + l2Pz*l2Pz));    
    svFitStandalone::kDecayType l1Type, l2Type;
    if ( std::string(argv[2]) == "EMu" ) {
      l1Type = svFitStandalone::kPrompt;
      l2Type = svFitStandalone::kTauToMuDecay;
    } else if ( std::string(argv[2]) == "MuE" ) {
      l1Type = svFitStandalone::kPrompt;
      l2Type = svFitStandalone::kTauToElecDecay;
    } else if ( std::string(argv[2]) == "MuTau" || std::string(argv[2]) == "ETau" ) {
      l1Type = svFitStandalone::kPrompt;
      l2Type = svFitStandalone::kTauToHadDecay;
    } else {
      std::cerr << "Error: Invalid channel = " << std::string(argv[2]) << " !!" << std::endl;
      std::cerr << "(some customization of this code will be needed for your analysis)" << std::endl;
      assert(0);
    }
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2));
    // construct the class object from the minimal necesarry information
    SVfitStandaloneAlgorithmLFV algo(measuredTauLeptons, measuredMET, covMET, 1);
    // apply customized configurations if wanted (examples are given below)
    algo.maxObjFunctionCalls(5000);
    //algo.addLogM(false);
    //algo.metPower(0.5)
    // minuit fit method
    //algo.fit();
    // integration by VEGAS (default)
    algo.integrateVEGAS();
    // integration by markov chain MC
    //algo.integrateMarkovChain();
    // retrieve the results upon success
    std::cout << "... m truth : " << mTrue << std::endl;
    if ( algo.isValidSolution() ) {
      std::cout << "... m svfit : " << algo.mass() << "+/-" << algo.massUncert() << std::endl; // return value is in units of GeV
    } else {
      std::cout << "... m svfit : ---" << std::endl;
    }
  }
  return;
}

int main(int argc, char* argv[]) 
{
  //eventsFromTree(argc, argv);
  singleEvent();
  return 0;
}
