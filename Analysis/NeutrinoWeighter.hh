#ifndef NeutrinoWeighter_h
#define NeutrinoWeighter_h

#include <iostream>
#include <random>
#include "Rivet/Projections/FinalState.hh"

//#include "TObject.h"
//#include "TRandom3.h"
//#include "FourMomentum.h"
//#include "TH1F.h"

namespace Rivet {

class NeutrinoWeighter: public FinalState {
private:

  std::vector< double > m_nu_eta;
  std::vector< double > m_nu_sinh;
  std::vector< double > m_nu_cosh;
  std::vector< double > m_nubar_eta;
  std::vector< double > m_nubar_sinh;
  std::vector< double > m_nubar_cosh;
  std::vector< double > m_top_smears;
  std::vector< double > m_W_smears;
  FourMomentum m_top;
  FourMomentum m_tbar;
  FourMomentum m_ttbar;
  FourMomentum m_b;
  FourMomentum m_bbar;
  FourMomentum m_nu;
  FourMomentum m_nubar;
  double m_weight_max;
//  TRandom3 m_random;

  bool m_do_both_pairings;
  bool m_do_chi2;
  bool m_do_crystalball_smearing;
  bool m_include_x;
  bool m_include_y;
  bool m_include_phi;
  bool m_include_mWp;
  bool m_include_mWn;
  bool m_include_mtop;
  bool m_include_mtbar;

//  TH1F pt_25_50;
//  TH1F pt_50_100;
//  TH1F pt_100_200;
//  TH1F pt_200_inf;

public:

  std::vector< double > GetNuEta(){ return m_nu_eta; };
  std::vector< double > GetNubarEta(){ return m_nubar_eta; };
  void SetNuEta(std::vector< double > etas)     { m_nu_eta = etas; m_nubar_eta = etas;};
  void SetNuSinhEta(std::vector< double > etas) { m_nu_sinh = etas; m_nubar_sinh = etas;};
  void SetNuCoshEta(std::vector< double > etas) { m_nu_cosh = etas; m_nubar_cosh = etas;};
  void SetTopMass(std::vector< double > top_mass){ m_top_smears = top_mass;};
  void SetWMass(std::vector< double > W_mass){ m_W_smears = W_mass;};

  NeutrinoWeighter(int setting, int event_number);  
  virtual ~NeutrinoWeighter(){};  
  double Reconstruct(FourMomentum lepton_pos, FourMomentum lepton_neg, FourMomentum jet_1, FourMomentum jet_2, double met_ex, double met_ey, double met_phi);
  void calculateWeight(FourMomentum lepton_pos, FourMomentum lepton_neg, FourMomentum b1, FourMomentum b2, double met_ex, double met_ey, double met_phi, double mtop, double mtbar, double mWp, double mWn);
  std::vector<FourMomentum> solveForNeutrinoEta(FourMomentum* lepton, FourMomentum* bJet, int index, int index_type, double mtop, double mW);
  //double neutrino_weight(FourMomentum neutrino1, FourMomentum neutrino2, double met_ex, double met_ey, double met_phi);
  double neutrino_weight(FourMomentum neutrino1, FourMomentum neutrino2, double met_ex, double met_ey, double met_phi, FourMomentum lep_pos, FourMomentum lep_neg, FourMomentum b1, FourMomentum b2, double mtop, double mtbar, double mWp, double mWn);

//  TString randomString(size_t l, std::string charIndex);

  FourMomentum GetTop(){   return m_top;   };
  FourMomentum GetTbar(){  return m_tbar;  };
  FourMomentum GetTtbar(){ return m_ttbar; };
  FourMomentum GetB(){     return m_b;     };
  FourMomentum GetBbar(){  return m_bbar;  };
  FourMomentum GetNu(){    return m_nu;    };
  FourMomentum GetNubar(){ return m_nubar; };
  double GetWeight(){ return m_weight_max; };
  void Reset();
  void RecalculateEtas(double pos_lep_eta, double neg_lep_eta);
  void DoBothPairings(bool do_pairing){ m_do_both_pairings = do_pairing; };
  //void SetTrueParameters(double true_nu_eta, double true_nubar_eta, double true_top_m, double true_tbar_m, double true_wp_m, double true_wp_n);
  float GetCrystalBallWeight(float jet_pt);
};

}

#endif
 
 
