// -*- C++ -*-
#include <cassert>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"

#include "WFinder221.hh"
#include "PseudoTop.hh"
#include "NeutrinoWeighter.hh"

// #define _SHOW_MASS_CUT_MSGS_ _SHOW_MASS_CUT_MSGS_
// #define _NW_DEBUG_MSGS_ _NW_DEBUG_MSGS_

#define _EcmsTag_        13000.
#define _TopQMass_       173.1
#define _NconsidJetsTag_     4


namespace Rivet {


  /// CMS 8 TeV dilepton channel ttbar spin correlations and polarisation analysis
  class MC13TeV_SPINCORR_emu_TOP_tt : public Analysis {

#include "NLOHisto1D.cc"

  private:

    NLOHisto1DPtr _h_dphil, _h_mll, _h_njets, _h_nbjets, _h_pTb1, _h_pTb2, _h_pTj1, _h_pTj2, _h_DRj, _h_dphij, _h_mjj, _h_ht, _h_ht2, _h_meff, _h_MT, _h_MT0j, _h_MT1j, _h_MT2j, _h_MT3j, _h_MT4j, _h_MTooj, _h_met, _h_eta_el, _h_eta_mu, _h_topPt, _h_topPtLead, _h_topPtSubLead, _h_tMass, _h_ttbarMass, _h_v_topPt, _h_v_topPtLead, _h_v_topPtSubLead, _h_v_tMass, _h_v_ttbarMass, _h_v_ttbarMassATLAS, _h_v_ttbarMassATLAS2, _h_selected_events, _h_xs, _h_xscut, _h_eff;
    Histo1DPtr _h_topreco;

    std::vector<double> v_bins_dphi;
    std::vector<double> v_bins_mtt;
    std::vector<double> v_bins_mtt2;
    //const vector<double> _bins_dphi = {0., 5.*M_PI/60., 10.*M_PI/60., 15.*M_PI/60., 20.*M_PI/60., 25.*M_PI/60., 30.*M_PI/60., 35.*M_PI/60., 40.*M_PI/60., 45.*M_PI/60., 50.*M_PI/60., 55.*M_PI/60., M_PI};

    std::vector<std::pair<double,double> > v_cutmT, v_cutmtt;
    std::vector<std::pair<double,double> > v_cutpTlb, v_cutpTt;
    std::vector<string> v_cutmT_spec, v_cutmtt_spec, v_cutpT_spec;

    NLOHisto1DPtr              h_cutmT_xsecs,     h_cutmtt_xsecs;
    std::vector<NLOHisto1DPtr> v_cutmT_histdphil, v_cutmtt_histdphil;
    std::vector<NLOHisto1DPtr> v_cutmT_histmT,    v_cutmtt_histmT;
    std::vector<NLOHisto1DPtr> v_cutmT_histmtt,   v_cutmtt_histmtt;
    std::vector<NLOHisto1DPtr> v_cutmT_histpTlb,  v_cutmtt_histpTlb;
    std::vector<NLOHisto1DPtr> v_cutmT_histpTt,   v_cutmtt_histpTt;

    // Inclusive definitions of MT
    NLOHisto1DPtr              h_cutmT0j_xsecs, h_cutmT1j_xsecs, h_cutmT2j_xsecs, h_cutmT3j_xsecs, h_cutmT4j_xsecs, h_cutmTooj_xsecs;
    std::vector<NLOHisto1DPtr> v_cutmT0j_histdphil, v_cutmT1j_histdphil, v_cutmT2j_histdphil, v_cutmT3j_histdphil, v_cutmT4j_histdphil, v_cutmTooj_histdphil;
    std::vector<NLOHisto1DPtr> v_cutmT0j_histmT0j, v_cutmT1j_histmT1j, v_cutmT2j_histmT2j, v_cutmT3j_histmT3j, v_cutmT4j_histmT4j, v_cutmTooj_histmTooj;

    NLOHisto1DPtr              h_cutpTlb_xsecs,     h_cutpTt_xsecs;
    std::vector<NLOHisto1DPtr> v_cutpTlb_histdphil, v_cutpTt_histdphil;
    std::vector<NLOHisto1DPtr> v_cutpTlb_histmT,    v_cutpTt_histmT;
    std::vector<NLOHisto1DPtr> v_cutpTlb_histmtt,   v_cutpTt_histmtt;
    std::vector<NLOHisto1DPtr> v_cutpTlb_histpTlb,  v_cutpTt_histpTlb;
    std::vector<NLOHisto1DPtr> v_cutpTlb_histpTt,   v_cutpTt_histpTt;

    // Top reco alternatives: pseudo-top, neutrino-weighted
    NLOHisto1DPtr              h_cutmtt_ps_xsecs,     h_cutmtt_nw_xsecs;
    std::vector<NLOHisto1DPtr> v_cutmtt_ps_histdphil, v_cutmtt_nw_histdphil;
    std::vector<NLOHisto1DPtr> v_cutmtt_ps_histmT,    v_cutmtt_nw_histmT;
    std::vector<NLOHisto1DPtr> v_cutmtt_ps_histmtt,   v_cutmtt_nw_histmtt;
    std::vector<NLOHisto1DPtr> v_cutmtt_ps_histpTlb,  v_cutmtt_nw_histpTlb;
    std::vector<NLOHisto1DPtr> v_cutmtt_ps_histpTt,   v_cutmtt_nw_histpTt;

    NLOHisto1DPtr              h_cutpTt_ps_xsecs,     h_cutpTt_nw_xsecs;
    std::vector<NLOHisto1DPtr> v_cutpTt_ps_histdphil, v_cutpTt_nw_histdphil;
    std::vector<NLOHisto1DPtr> v_cutpTt_ps_histmT,    v_cutpTt_nw_histmT;
    std::vector<NLOHisto1DPtr> v_cutpTt_ps_histmtt,   v_cutpTt_nw_histmtt;
    std::vector<NLOHisto1DPtr> v_cutpTt_ps_histpTlb,  v_cutpTt_nw_histpTlb;
    std::vector<NLOHisto1DPtr> v_cutpTt_ps_histpTt,   v_cutpTt_nw_histpTt;

    std::map<std::string, NLOHisto1DPtr> m_histfracs;


    void initMassCutHistos() {
      v_cutpT_spec ={"0_40", "0_60", "0_80", "0_100", "100_300", "300_oo", "0_200"};
      v_cutmtt_spec={"0_450", "450_550", "550_800", "800_oo", "0_400", "0_500"};
      v_cutmT_spec ={"0_400", "400_500", "500_750", "750_oo", "0_350", "0_450",
	"0_oo"};
      v_cutpTt={std::make_pair(0.,40.), std::make_pair(0.,60.), std::make_pair(0.,80.),
	std::make_pair(0.,100.), std::make_pair(100.,300.),
	std::make_pair(300.,_EcmsTag_),	std::make_pair(0.,200.)};
      v_cutpTlb={std::make_pair(0.,40.), std::make_pair(0.,60.), std::make_pair(0.,80.),
	std::make_pair(0.,100.), std::make_pair(100.,300.),
	std::make_pair(300.,_EcmsTag_),	std::make_pair(0.,200.)};
      v_cutmtt={std::make_pair(0.,450.), std::make_pair(450.,550.),
	std::make_pair(550.,800.), std::make_pair(800.,_EcmsTag_),
	std::make_pair(0.,400.), std::make_pair(0.,500.)};
      v_cutmT={std::make_pair(0.,400.), std::make_pair(400.,500.),
	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
	std::make_pair(0.,350.), std::make_pair(0.,450.),
	std::make_pair(0.,_EcmsTag_)};
//      v_cutmT0j={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
//      v_cutmT1j={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
//      v_cutmT2j={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
//      v_cutmT3j={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
//      v_cutmT4j={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
//      v_cutmTooj={std::make_pair(0.,400.), std::make_pair(400.,500.),
//	std::make_pair(500.,750.), std::make_pair(750.,_EcmsTag_),
//	std::make_pair(0.,350.), std::make_pair(0.,450.),
//	std::make_pair(0.,_EcmsTag_)};
      size_t num(v_cutmtt.size());
      h_cutmtt_xsecs=bookNLOHisto1D("Axsecs_Cmtt", num, 0., num);
      h_cutmtt_ps_xsecs=bookNLOHisto1D("Axsecs_Cmtt_ps", num, 0., num);
      h_cutmtt_nw_xsecs=bookNLOHisto1D("Axsecs_Cmtt_nw", num, 0., num);
      for(size_t i=0; i<num; ++i) {
	v_cutmtt_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_Cmtt_"+v_cutmtt_spec[i], v_bins_dphi));
	v_cutmtt_histmtt.push_back
	  (bookNLOHisto1D("Amtt_Cmtt_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_histmT.push_back
	  (bookNLOHisto1D("AmT_Cmtt_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_histpTt.push_back
	  (bookNLOHisto1D("ApTt_Cmtt_"+v_cutmtt_spec[i], 100, 0., 600.));
	v_cutmtt_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_Cmtt_"+v_cutmtt_spec[i], 100, 0., 600.));
	// Pseudo
	v_cutmtt_ps_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_Cmtt_ps_"+v_cutmtt_spec[i], v_bins_dphi));
	v_cutmtt_ps_histmtt.push_back
	  (bookNLOHisto1D("Amtt_Cmtt_ps_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_ps_histmT.push_back
	  (bookNLOHisto1D("AmT_Cmtt_ps_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_ps_histpTt.push_back
	  (bookNLOHisto1D("ApTt_Cmtt_ps_"+v_cutmtt_spec[i], 100, 0., 600.));
	v_cutmtt_ps_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_Cmtt_ps_"+v_cutmtt_spec[i], 100, 0., 600.));
	// Neutrino-weighted
	v_cutmtt_nw_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_Cmtt_nw_"+v_cutmtt_spec[i], v_bins_dphi));
	v_cutmtt_nw_histmtt.push_back
	  (bookNLOHisto1D("Amtt_Cmtt_nw_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_nw_histmT.push_back
	  (bookNLOHisto1D("AmT_Cmtt_nw_"+v_cutmtt_spec[i], 100, 0., 1000.));
	v_cutmtt_nw_histpTt.push_back
	  (bookNLOHisto1D("ApTt_Cmtt_nw_"+v_cutmtt_spec[i], 100, 0., 600.));
	v_cutmtt_nw_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_Cmtt_nw_"+v_cutmtt_spec[i], 100, 0., 600.));
      }
      num=v_cutmT.size();
      h_cutmT_xsecs=bookNLOHisto1D("Axsecs_CmT", num, 0., num);
      for(size_t i=0; i<num; ++i) {
	v_cutmT_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmT_histmtt.push_back
	  (bookNLOHisto1D("Amtt_CmT_"+v_cutmT_spec[i], 100, 0., 1000.));
	v_cutmT_histmT.push_back
	  (bookNLOHisto1D("AmT_CmT_"+v_cutmT_spec[i], 100, 0., 1000.));
	v_cutmT_histpTt.push_back
	  (bookNLOHisto1D("ApTt_CmT_"+v_cutmT_spec[i], 100, 0., 600.));
	v_cutmT_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_CmT_"+v_cutmT_spec[i], 100, 0., 600.));
      }
      num=v_cutmT.size();
      h_cutmT0j_xsecs=bookNLOHisto1D("Axsecs_CmT0j", num, 0., num);
      h_cutmT1j_xsecs=bookNLOHisto1D("Axsecs_CmT1j", num, 0., num);
      h_cutmT2j_xsecs=bookNLOHisto1D("Axsecs_CmT2j", num, 0., num);
      h_cutmT3j_xsecs=bookNLOHisto1D("Axsecs_CmT3j", num, 0., num);
      h_cutmT4j_xsecs=bookNLOHisto1D("Axsecs_CmT4j", num, 0., num);
      h_cutmTooj_xsecs=bookNLOHisto1D("Axsecs_CmTooj", num, 0., num);
      for(size_t i=0; i<num; ++i) {
	v_cutmT0j_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT0j_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmT1j_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT1j_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmT2j_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT2j_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmT3j_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT3j_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmT4j_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmT4j_"+v_cutmT_spec[i], v_bins_dphi));
	v_cutmTooj_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CmTooj_"+v_cutmT_spec[i], v_bins_dphi));
        v_cutmT0j_histmT0j.push_back
          (bookNLOHisto1D("AmT0j_CmT0j_"+v_cutmT_spec[i], 100, 0., 1000.));
        v_cutmT1j_histmT1j.push_back
          (bookNLOHisto1D("AmT1j_CmT1j_"+v_cutmT_spec[i], 100, 0., 1000.));
        v_cutmT2j_histmT2j.push_back
          (bookNLOHisto1D("AmT2j_CmT2j_"+v_cutmT_spec[i], 100, 0., 1000.));
        v_cutmT3j_histmT3j.push_back
          (bookNLOHisto1D("AmT3j_CmT3j_"+v_cutmT_spec[i], 100, 0., 1000.));
        v_cutmT4j_histmT4j.push_back
          (bookNLOHisto1D("AmT4j_CmT4j_"+v_cutmT_spec[i], 100, 0., 1000.));
        v_cutmTooj_histmTooj.push_back
          (bookNLOHisto1D("AmTooj_CmTooj_"+v_cutmT_spec[i], 100, 0., 1000.));
      }
      num=v_cutpTt.size();
      h_cutpTt_xsecs=bookNLOHisto1D("Axsecs_CpTt", num, 0., num);
      h_cutpTt_ps_xsecs=bookNLOHisto1D("Axsecs_CpTt_ps", num, 0., num);
      h_cutpTt_nw_xsecs=bookNLOHisto1D("Axsecs_CpTt_nw", num, 0., num);
      for(size_t i=0; i<num; ++i) {
	v_cutpTt_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CpTt_"+v_cutpT_spec[i], v_bins_dphi));
	v_cutpTt_histmtt.push_back
	  (bookNLOHisto1D("Amtt_CpTt_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_histmT.push_back
	  (bookNLOHisto1D("AmT_CpTt_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_histpTt.push_back
	  (bookNLOHisto1D("ApTt_CpTt_"+v_cutpT_spec[i], 100, 0., 600.));
	v_cutpTt_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_CpTt_"+v_cutpT_spec[i], 100, 0., 600.));
	// Pseudo
	v_cutpTt_ps_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CpTt_ps_"+v_cutpT_spec[i], v_bins_dphi));
	v_cutpTt_ps_histmtt.push_back
	  (bookNLOHisto1D("Amtt_CpTt_ps_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_ps_histmT.push_back
	  (bookNLOHisto1D("AmT_CpTt_ps_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_ps_histpTt.push_back
	  (bookNLOHisto1D("ApTt_CpTt_ps_"+v_cutpT_spec[i], 100, 0., 600.));
	v_cutpTt_ps_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_CpTt_ps_"+v_cutpT_spec[i], 100, 0., 600.));
	// Neutrino-weighted	
	v_cutpTt_nw_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CpTt_nw_"+v_cutpT_spec[i], v_bins_dphi));
	v_cutpTt_nw_histmtt.push_back
	  (bookNLOHisto1D("Amtt_CpTt_nw_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_nw_histmT.push_back
	  (bookNLOHisto1D("AmT_CpTt_nw_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTt_nw_histpTt.push_back
	  (bookNLOHisto1D("ApTt_CpTt_nw_"+v_cutpT_spec[i], 100, 0., 600.));
	v_cutpTt_nw_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_CpTt_nw_"+v_cutpT_spec[i], 100, 0., 600.));
      }
      num=v_cutpTlb.size();
      h_cutpTlb_xsecs=bookNLOHisto1D("Axsecs_CpTlb", num, 0., num);
      for(size_t i=0; i<num; ++i) {
	v_cutpTlb_histdphil.push_back
	  (bookNLOHisto1D("Adphi_l_CpTlb_"+v_cutpT_spec[i], v_bins_dphi));
	v_cutpTlb_histmtt.push_back
	  (bookNLOHisto1D("Amtt_CpTlb_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTlb_histmT.push_back
	  (bookNLOHisto1D("AmT_CpTlb_"+v_cutpT_spec[i], 100, 0., 1000.));
	v_cutpTlb_histpTt.push_back
	  (bookNLOHisto1D("ApTt_CpTlb_"+v_cutpT_spec[i], 100, 0., 600.));
	v_cutpTlb_histpTlb.push_back
	  (bookNLOHisto1D("ApTlb_CpTlb_"+v_cutpT_spec[i], 100, 0., 600.));
      }

      cout<<" CutData info:  "<<v_cutmtt.size()<<" "<<v_cutmT.size()<<" : "
	  <<v_cutmtt_spec.size()<<" "<<v_cutmT_spec.size()<<" : "
	  <<v_cutmtt_histmT.size()<<" "<<v_cutmT_histmT.size()<<endl<<"    ";
      for(size_t i=0; i<v_cutmtt_histmT.size(); ++i) {
	cout<<v_cutmtt_histmT[i]->path()<<", ";} cout<<endl;
      for(size_t i=0; i<v_cutmT_histdphil.size(); ++i) {
	cout<<v_cutmT_histdphil[i]->path()<<", ";} cout<<endl;

      m_histfracs["DRatio_mT:mtt"]        = bookNLOHisto1D("DRatio_mT:mtt",200,0.,2.);
      m_histfracs["DRatio_mttT:mtt"]      = bookNLOHisto1D("DRatio_mttT:mtt",200,0.,2.);
      m_histfracs["DRatio_mTnn:mttT"]     =
	bookNLOHisto1D("DRatio_mTnn:mttT",200,0.,2.);
      m_histfracs["DRatio_mTnn:mT"]       = bookNLOHisto1D("DRatio_mTnn:mT",200,0.,2.);
      m_histfracs["DRatio_MT:mT"]         = bookNLOHisto1D("DRatio_MT:mT",200,0.,2.);
      m_histfracs["DRatio_MTp:mT"]        = bookNLOHisto1D("DRatio_MTp:mT",200,0.,2.);
      m_histfracs["DRatio_pTnn:pTmiss"]   =
	bookNLOHisto1D("DRatio_pTnn:pTmiss",200,0.,2.);
      m_histfracs["DRatio_pTllbb:pTmiss"] =
	bookNLOHisto1D("DRatio_pTllbb:pTmiss",200,0.,2.);
      m_histfracs["DRatio_mtt:mtt_ps"]    = bookNLOHisto1D("DRatio_mtt:mtt_ps",200,0.,2.);
      m_histfracs["DRatio_mtt:mtt_nw"]    = bookNLOHisto1D("DRatio_mtt:mtt_nw",200,0.,2.);
      m_histfracs["DRatio_mtt_nw:mtt_truth"]=bookNLOHisto1D("DRatio_mtt_nw:mtt_truth",200,0.,2.);
      m_histfracs["DRatio_pTt:pTt_ps"]    = bookNLOHisto1D("DRatio_pTt:pTt_ps",200,0.,2.);
      m_histfracs["DRatio_pTt:pTt_nw"]    = bookNLOHisto1D("DRatio_pTt:pTt_nw",200,0.,2.);

      m_histfracs["DRatio_mT:mT0j"] = bookNLOHisto1D("DRatio_mT:mT0j", 200, 0., 2.);
      m_histfracs["DRatio_mT:mT1j"] = bookNLOHisto1D("DRatio_mT:mT1j", 200, 0., 2.);
      m_histfracs["DRatio_mT:mT2j"] = bookNLOHisto1D("DRatio_mT:mT2j", 200, 0., 2.);
      m_histfracs["DRatio_mT:mT3j"] = bookNLOHisto1D("DRatio_mT:mT3j", 200, 0., 2.);
      m_histfracs["DRatio_mT:mT4j"] = bookNLOHisto1D("DRatio_mT:mT4j", 200, 0., 2.);
      m_histfracs["DRatio_mT:mTooj"] = bookNLOHisto1D("DRatio_mT:mTooj", 200, 0., 2.);

    }


    void fillMassCutHistos(const Event& evt, const double& dphi,
			   const double& mtt, const double& mT,
                           const double& mT0j, const double& mT1j,
                           const double& mT2j, const double& mT3j,
                           const double& mT4j, const double& mTooj,
			   const double& pTt, const double& pTlb,
			   const double& mtt_ps, const double& mtt_nw,
			   const double& pTt_ps, const double& pTt_nw) {
      for(size_t i=0; i<v_cutmtt.size(); ++i) {
	if(v_cutmtt[i].first<=mtt && mtt<v_cutmtt[i].second) {
	  h_cutmtt_xsecs->fill(i, evt);
	  v_cutmtt_histdphil[i]->fill(dphi, evt);
	  v_cutmtt_histmtt[i]->fill(mtt, evt);
	  v_cutmtt_histmT[i]->fill(mT, evt);
	  v_cutmtt_histpTt[i]->fill(pTt, evt);
	  v_cutmtt_histpTlb[i]->fill(pTlb, evt);
	}
	if(v_cutmtt[i].first<=mtt_ps && mtt_ps<v_cutmtt[i].second && mtt_ps > 0.) {
	  h_cutmtt_ps_xsecs->fill(i, evt);
	  v_cutmtt_ps_histdphil[i]->fill(dphi, evt);
	  v_cutmtt_ps_histmtt[i]->fill(mtt_ps, evt);
	  v_cutmtt_ps_histmT[i]->fill(mT, evt);
	  v_cutmtt_ps_histpTt[i]->fill(pTt_ps, evt);
	  v_cutmtt_ps_histpTlb[i]->fill(pTlb, evt);
	}
	if(v_cutmtt[i].first<=mtt_nw && mtt_nw<v_cutmtt[i].second && mtt_nw > 0.) {
	  h_cutmtt_nw_xsecs->fill(i, evt);
	  v_cutmtt_nw_histdphil[i]->fill(dphi, evt);
	  v_cutmtt_nw_histmtt[i]->fill(mtt_nw, evt);
	  v_cutmtt_nw_histmT[i]->fill(mT, evt);
	  v_cutmtt_nw_histpTt[i]->fill(pTt_nw, evt);
	  v_cutmtt_nw_histpTlb[i]->fill(pTlb, evt);
	}
	}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT && mT<v_cutmT[i].second) {
	  h_cutmT_xsecs->fill(i, evt);
	  v_cutmT_histdphil[i]->fill(dphi, evt);
	  v_cutmT_histmtt[i]->fill(mtt, evt);
	  v_cutmT_histmT[i]->fill(mT, evt);
	  v_cutmT_histpTt[i]->fill(pTt, evt);
	  v_cutmT_histpTlb[i]->fill(pTlb, evt);
	}}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT0j && mT0j<v_cutmT[i].second) {
	  h_cutmT0j_xsecs->fill(i, evt);
	  v_cutmT0j_histdphil[i]->fill(dphi, evt);
	  v_cutmT0j_histmT0j[i]->fill(mT0j, evt);
        }}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT1j && mT1j<v_cutmT[i].second) {
	  h_cutmT1j_xsecs->fill(i, evt);
	  v_cutmT1j_histdphil[i]->fill(dphi, evt);
	  v_cutmT1j_histmT1j[i]->fill(mT1j, evt);
        }}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT2j && mT2j<v_cutmT[i].second) {
	  h_cutmT2j_xsecs->fill(i, evt);
	  v_cutmT2j_histdphil[i]->fill(dphi, evt);
	  v_cutmT2j_histmT2j[i]->fill(mT2j, evt);
        }}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT3j && mT3j<v_cutmT[i].second) {
	  h_cutmT3j_xsecs->fill(i, evt);
	  v_cutmT3j_histdphil[i]->fill(dphi, evt);
	  v_cutmT3j_histmT3j[i]->fill(mT3j, evt);
        }}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mT4j && mT4j<v_cutmT[i].second) {
	  h_cutmT4j_xsecs->fill(i, evt);
	  v_cutmT4j_histdphil[i]->fill(dphi, evt);
	  v_cutmT4j_histmT4j[i]->fill(mT4j, evt);
        }}
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	if(v_cutmT[i].first<=mTooj && mTooj<v_cutmT[i].second) {
	  h_cutmTooj_xsecs->fill(i, evt);
	  v_cutmTooj_histdphil[i]->fill(dphi, evt);
	  v_cutmTooj_histmTooj[i]->fill(mTooj, evt);
        }}
      for(size_t i=0; i<v_cutpTt.size(); ++i) {
	if(v_cutpTt[i].first<=pTt && pTt<v_cutpTt[i].second) {
	  h_cutpTt_xsecs->fill(i, evt);
	  v_cutpTt_histdphil[i]->fill(dphi, evt);
	  v_cutpTt_histmtt[i]->fill(mtt, evt);
	  v_cutpTt_histmT[i]->fill(mT, evt);
	  v_cutpTt_histpTt[i]->fill(pTt, evt);
	  v_cutpTt_histpTlb[i]->fill(pTlb, evt);
	}
	if(v_cutpTt[i].first<=pTt_ps && pTt_ps<v_cutpTt[i].second && pTt_ps > 0.) {
	  h_cutpTt_ps_xsecs->fill(i, evt);
	  v_cutpTt_ps_histdphil[i]->fill(dphi, evt);
	  v_cutpTt_ps_histmtt[i]->fill(mtt_ps, evt);
	  v_cutpTt_ps_histmT[i]->fill(mT, evt);
	  v_cutpTt_ps_histpTt[i]->fill(pTt_ps, evt);
	  v_cutpTt_ps_histpTlb[i]->fill(pTlb, evt);
	}
	if(v_cutpTt[i].first<=pTt_nw && pTt_nw<v_cutpTt[i].second && pTt_nw > 0.) {
	  h_cutpTt_nw_xsecs->fill(i, evt);
	  v_cutpTt_nw_histdphil[i]->fill(dphi, evt);
	  v_cutpTt_nw_histmtt[i]->fill(mtt_nw, evt);
	  v_cutpTt_nw_histmT[i]->fill(mT, evt);
	  v_cutpTt_nw_histpTt[i]->fill(pTt_nw, evt);
	  v_cutpTt_nw_histpTlb[i]->fill(pTlb, evt);
	}
	}
      for(size_t i=0; i<v_cutpTlb.size(); ++i) {
	if(v_cutpTlb[i].first<=pTlb && pTlb<v_cutpTlb[i].second) {
	  h_cutpTlb_xsecs->fill(i, evt);
	  v_cutpTlb_histdphil[i]->fill(dphi, evt);
	  v_cutpTlb_histmtt[i]->fill(mtt, evt);
	  v_cutpTlb_histmT[i]->fill(mT, evt);
	  v_cutpTlb_histpTt[i]->fill(pTt, evt);
	  v_cutpTlb_histpTlb[i]->fill(pTlb, evt);
	}}
    }




  public:

    /// Constructor
    /// Constructor
    MC13TeV_SPINCORR_emu_TOP_tt()
      : Analysis("MC13TeV_SPINCORR_emu_TOP_tt")
    {    }

    /// Book histograms and initialise projections
    void init() {

      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for dressed electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons,
				       0.1, Cuts::open(), true, false);
      addProjection(dressed_electrons, "DressedElectrons");

      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons,
				   0.1, Cuts::open(), true, false);
      addProjection(dressed_muons, "DressedMuons");

      // W finders for identifying the ttbar system
      WFinder221 wfel(fs, Cuts::open(), PID::ELECTRON, 0*GeV, _EcmsTag_*GeV,
		      0*GeV, 0.1, WFinder221::CLUSTERNODECAY,
		      WFinder221::NOTRACK);
      addProjection(wfel, "WFinder_el");
//      WFinder221 wfek(wfel.remainingFinalState(),
//		      Cuts::open(), PID::ELECTRON, PID::MUON, 0*GeV, _EcmsTag_*GeV, 0*GeV,
//		      0.1, WFinder221::CLUSTERNODECAY, WFinder221::NOTRACK);
//      addProjection(wfek, "WFinder_ek");
      WFinder221 wfmu(fs, Cuts::open(), PID::MUON, 0*GeV, _EcmsTag_*GeV, 0*GeV,
		      0.1, WFinder221::CLUSTERNODECAY, WFinder221::NOTRACK);
      addProjection(wfmu, "WFinder_mu");

      // Pseudo-top finder
      declare(PseudoTop(0.1, 25, 2.5, 0.4, 25, 2.5), "ttbar");

      // Projection for anti-kt jets (use fs without dressed leptons / neutrinos as input)
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(dressed_electrons);
      vfs.addVetoOnThisFinalState(dressed_muons);
      vfs.addVetoPairId(PID::NU_E);
      vfs.addVetoPairId(PID::NU_MU);
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.4), "Jets");

      addProjection(MissingMomentum(fs), "MissingET");

      // Visible final state to compute HT, MT, m_eff, ...
      addProjection(VisibleFinalState(-5, 5),"visfs");

      // Booking of histograms
//      v_bins_dphi = {0.,5.*PI/60.,10.*PI/60.,15.*PI/60.,20.*PI/60.,
//	25.*PI/60.,30.*PI/60.,35.*PI/60.,40.*PI/60.,45.*PI/60.,
//	50.*PI/60.,55.*PI/60.,PI};
      v_bins_dphi = {0.,1.*PI/10.,2.*PI/10.,3.*PI/10.,4.*PI/10.,
	5.*PI/10.,6.*PI/10.,7.*PI/10.,8.*PI/10.,9.*PI/10.,
	PI};
      //const vector<double> _bins_dphi = {0., 5.*M_PI/60., 10.*M_PI/60., 15.*M_PI/60., 20.*M_PI/60., 25.*M_PI/60., 30.*M_PI/60., 35.*M_PI/60., 40.*M_PI/60., 45.*M_PI/60., 50.*M_PI/60., 55.*M_PI/60., M_PI};
      v_bins_mtt = {350.,450.,700.,1000.,1500.,6000.};
      v_bins_mtt2 = {350.,450.,550.,800.,2000.};
      
      // This histogram is independent of the parton-level information, and is an addition to the original analysis.
      // It is compared to the same data as the parton-level delta_phi histogram d02-x01-y01.
      _h_dphil = bookNLOHisto1D("dphi_l", v_bins_dphi);
      _h_mll = bookNLOHisto1D("mll", 70., 0., 700.);
      _h_njets = bookNLOHisto1D("njets", 8., -0.5, 7.5);
      _h_nbjets = bookNLOHisto1D("nbjets", 6., -0.5, 5.5);
      _h_pTb1 = bookNLOHisto1D("pTb1", 40, 0., 800.);
      _h_pTb2 = bookNLOHisto1D("pTb2", 40, 0., 800.);
      _h_pTj1 = bookNLOHisto1D("pTj1", 40, 0., 800.);
      _h_pTj2 = bookNLOHisto1D("pTj2", 50, 20., 270.);
      _h_DRj = bookNLOHisto1D("DRj", 20., 0., 7.0);
      _h_dphij = bookNLOHisto1D("dphij", 20., 0., M_PI);
      _h_mjj = bookNLOHisto1D("mjj", 70., 0., 700.);
      _h_ht = bookNLOHisto1D("ht", 100., 0., 1000.);
      _h_ht2 = bookNLOHisto1D("ht2", 100., 0., 1000.);
      _h_meff = bookNLOHisto1D("meff", 100., 0., 1000.);
      _h_MT = bookNLOHisto1D("MT", 100., 0., 1000.);
      _h_MT0j = bookNLOHisto1D("mT0j", 100., 0., 1000.);
      _h_MT1j = bookNLOHisto1D("mT1j", 100., 0., 1000.);
      _h_MT2j = bookNLOHisto1D("mT2j", 100., 0., 1000.);
      _h_MT3j = bookNLOHisto1D("mT3j", 100., 0., 1000.);
      _h_MT4j = bookNLOHisto1D("mT4j", 100., 0., 1000.);
      _h_MTooj = bookNLOHisto1D("mTooj", 100., 0., 1000.);

      // Pseudo-top hists
      _h_topPt = bookNLOHisto1D("pseudo_pt_t", 100, 0., 600.);
      _h_topPtLead = bookNLOHisto1D("pseudo_pt_t_lead", 100, 0., 600.);
      _h_topPtSubLead = bookNLOHisto1D("pseudo_pt_t_sublead", 100, 0., 600.);
      _h_tMass = bookNLOHisto1D("pseudo_mass_t", 100, 0., 400.);
      _h_ttbarMass = bookNLOHisto1D("pseudo_ttbar_mass", 100., 0., 1000.);
      // Neutrino-weighted tops
      _h_v_topPt = bookNLOHisto1D("vw_pt_t", 100, 0., 600.);
      _h_v_topPtLead = bookNLOHisto1D("vw_pt_t_lead", 100, 0., 600.);
      _h_v_topPtSubLead = bookNLOHisto1D("vw_pt_t_sublead", 100, 0., 600.);
      _h_v_tMass = bookNLOHisto1D("vw_mass_t", 100, 0., 400.);
      _h_v_ttbarMass = bookNLOHisto1D("vw_ttbar_mass", 100., 0., 1000.);
      _h_v_ttbarMassATLAS = bookNLOHisto1D("vw_ttbar_massATLAS", v_bins_mtt);
      _h_v_ttbarMassATLAS2 = bookNLOHisto1D("vw_ttbar_massATLAS2", v_bins_mtt2);

      _h_met = bookNLOHisto1D("ETmiss", 40., 0, 250.);
      _h_eta_el = bookNLOHisto1D("eta_el", 20, -2.5, 2.5);
      _h_eta_mu = bookNLOHisto1D("eta_mu", 20, -2.5, 2.5);
      _h_selected_events = bookNLOHisto1D("h_selected_events", 7, 0.0, 7.0);

      _h_xs = bookNLOHisto1D("xs", 1, 0., 1.);
      _h_xscut = bookNLOHisto1D("xscut", 1, 0., 1.);
      _h_eff = bookNLOHisto1D("efficiency", 1, 0., 1.);
      _h_topreco = bookNLOHisto1D("reco_yield", 1, 0., 1.);

      // Initialize the cut info and associated histograms:
      initMassCutHistos();

    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Event& weight = event;
      _h_xs->fill( 0.5, weight );
      // Use particle-level leptons for the first histogram

      const DressedLeptons& _dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& _dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");

      const vector<DressedLepton> dressedels = _dressed_electrons.dressedLeptons();
      const vector<DressedLepton> dressedmus = _dressed_muons.dressedLeptons();

      _h_selected_events->fill(0.5, weight);

      // Jet requirements
      Jets goodJets;
      Jets bJets;

      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt();

      foreach (const Jet& jet, alljets) {
        if (jet.momentum().pT() > 25*GeV && fabs(jet.momentum().eta()) < 2.5) {
	  goodJets.push_back(jet);
	  if ( jet.containsBottom() ) bJets.push_back(jet);
        }
      }

      // Require two good jets, and at least one b-jet
      if (goodJets.size() < 2 || bJets.size() < 1) vetoEvent;
      _h_selected_events->fill(1.5, weight);

      // Remove electrons and muons that overlap with jets
      vector<DressedLepton> elecs, muons;

      foreach (const DressedLepton& elec, dressedels) {
	bool overlaps = false;
	foreach (const Jet& jet, goodJets) {
	  if (deltaR(elec, jet) < 0.4) overlaps = true;
	}
	if (!overlaps) elecs.push_back(elec);
      }
      foreach (const DressedLepton& muon, dressedmus) {
	bool overlaps = false;
        foreach (const Jet& jet, goodJets) {
          if (deltaR(muon, jet) < 0.4) overlaps = true;
        }
	if (!overlaps) muons.push_back(muon);
      }

      const size_t ndressedel = elecs.size();
      const size_t ndressedmu = muons.size();

      // For the particle-level histogram, require exactly one electron and exactly one muon, to select
      // the ttbar->emu channel. Note this means ttbar->emu events with additional PromptFinalState
      // dilepton pairs from the shower are vetoed - for PYTHIA8, this affects ~0.5% of events, so the
      // effect is well below the level of sensitivity of the measured distribution.
      if ( ndressedel == 1 && ndressedmu == 1 ) {

	_h_selected_events->fill(2.5, weight);

        // Opposite-charge leptons only
        if ( sameSign(elecs[0],muons[0]) ) {
          MSG_INFO("Error, e and mu have same charge, skipping event");
	  vetoEvent;
        }
	_h_selected_events->fill(3.5, weight);
	// pT and rapidity requirements for leptons
	if (elecs[0].pT() < 25*GeV || fabs(elecs[0].eta()) > 2.47 || ( fabs(elecs[0].eta()) > 1.37 && fabs(elecs[0].eta()) < 1.52 )) {
	  vetoEvent;
	}
	_h_selected_events->fill(4.5, weight);
	if (muons[0].pT() < 25*GeV || fabs(muons[0].eta()) > 2.5) vetoEvent;
	_h_selected_events->fill(5.5, weight);
	if (elecs[0].pT() < 27*GeV && muons[0].pT() < 27*GeV) vetoEvent;
        else {

	  const WFinder221& wfel=
	    applyProjection<WFinder221>(event, "WFinder_el");
	  const WFinder221& wfmu=
	    applyProjection<WFinder221>(event, "WFinder_mu");
#ifdef _SHOW_MASS_CUT_MSGS_
	  foreach(const Particle& p, wfel.particles()) {
	    const GenParticle* part=p.genParticle(); cout<<*part<<endl;}
	  cout<<"-------"<<endl;
	  foreach(const Particle& p, wfmu.particles()) {
	    const GenParticle* part=p.genParticle(); cout<<*part<<endl;}
#endif
	  if(wfel.empty() || wfmu.empty()) {
	    MSG_DEBUG("Less than two W bosons found. Discard event.");
	    vetoEvent;
	  }
	  if(wfel.size()==2 && wfmu.size()==2); else{
	    MSG_DEBUG("W boson reconstruction has failed. Discard event.");
	    vetoEvent;
	  }
	  if(elecs[0].charge()==wfel.constituentLepton().charge()) {
	    assert(elecs[0].momentum()==wfel.constituentLepton().momentum());
	    assert(muons[0].momentum()==wfmu.constituentLepton().momentum());
          }
	  else {
	    MSG_DEBUG("Wrong charge. Problem in the reco. Discard event.");
	    vetoEvent;
//	    assert(elecs[0].momentum()==wfek.constituentLepton().momentum());
//	    assert(elecs[1].momentum()==wfel.constituentLepton().momentum());
	  }
#ifdef _SHOW_MASS_CUT_MSGS_
	  cout<<"nW="<<wfel.size()<<"  Bs="<<wfel.bosons().size()
	      <<"  Ls="<<wfel.constituentLeptons().size()
	      <<"  Ns="<<wfel.constituentNeutrinos().size()<<endl;
#endif
	  const Particle& w1=wfel.boson(); const Particle w2=wfmu.boson();
	  const Particle& n1=wfel.constituentNeutrinos()[0];
	  const Particle& n2=wfmu.constituentNeutrinos()[0];
	  const MissingMomentum& missP = applyProjection<MissingMomentum>(event, "MissingET");
#ifdef _SHOW_MASS_CUT_MSGS_
	  cout<<"found: "<<w1.pid()<<" ("<<n1.pid()<<" "
	      <<wfel.constituentLepton().pid()<<")  and  "<<w2.pid()<<" ("
	      <<n2.pid()<<" "<<wfmu.constituentLepton().pid()<<")."<<endl;
	  cout<<"   w1:"<<w1.momentum()<<"\t m="<<w1.momentum().mass()<<endl;
	  cout<<"   n1:"<<n1.momentum()<<"\t m="<<n1.momentum().mass()<<endl;
	  cout<<"   e1:"<<wfel.constituentLepton().momentum()<<endl;
	  cout<<"-------- --------"<<endl;
	  cout<<"   w2:"<<w2.momentum()<<"\t m="<<w2.momentum().mass()<<endl;
	  cout<<"   n2:"<<n2.momentum()<<"\t m="<<n2.momentum().mass()<<endl;
	  cout<<"   e2:"<<wfmu.constituentLepton().momentum()<<endl;
	  cout<<"-------- --------"<<endl;
	  cout<<"     n1+n2="<<n1.momentum()+n2.momentum()
	      <<"  m="<<(n1.momentum()+n2.momentum()).mass()
	      <<"  pT="<<(n1.momentum()+n2.momentum()).pT()
	      <<"    versus    visP="<<missP.visibleMomentum()
	      <<"  m="<<missP.visibleMomentum().mass()
	      <<"  pT="<<missP.visibleMomentum().pT()<<endl;
#endif

	  // Fill the histograms
	  _h_xscut->fill( 0.5, weight );
          _h_eff->fill( 0.5, weight );
	  _h_selected_events->fill(6.5, weight);

	  double ptjet = 0.;
	  foreach ( const Jet& jet, goodJets) {
    //        std::cout << jet.momentum() << ", m2 = " << jet.mass2() << std::endl;
	    ptjet += jet.pT();
	  }
	  double ht = ptjet + elecs[0].pT() + muons[0].pT();
	  _h_ht->fill( ht, weight );

	  _h_njets->fill( goodJets.size(), weight);
	  _h_nbjets->fill( bJets.size(), weight);

	  double pTb_lead = bJets[0].pT();
	  double pTj_lead = goodJets[0].pT();
	  double pTj_sub = goodJets[1].pT();
	  double DRjets = deltaR(goodJets[0], goodJets[1]);
	  double dphijets = deltaPhi(goodJets[0], goodJets[1]);
	  double mjj = (goodJets[0].momentum() + goodJets[1].momentum()).mass();
	  _h_pTb1->fill( pTb_lead, weight);
	  if ( bJets.size() > 1 ) _h_pTb2->fill( bJets[1].pT(), weight);
	  _h_pTj1->fill( pTj_lead, weight);
	  _h_pTj2->fill( pTj_sub, weight);
	  _h_DRj->fill( DRjets, weight);
	  _h_dphij->fill( dphijets, weight);
	  _h_mjj->fill( mjj, weight);

	  double eta_el = elecs[0].eta();
	  double eta_mu = muons[0].eta();
	  _h_eta_el->fill( eta_el, weight);
	  _h_eta_mu->fill( eta_mu, weight);

	  // Missing momentum observables
	  double met = missP.vectorEt().mod();
	  _h_met->fill( met, weight);

	  ParticleVector visfs_particles = applyProjection<VisibleFinalState>(event, "visfs").particles();
	  double pT_vis = 0.;
	  foreach ( const Particle& p, visfs_particles ) {
	    pT_vis += p.pT();
	  }
	  _h_ht2->fill( pT_vis, weight);
	  double m_eff = pT_vis + (- missP.visibleMomentum()).pT();
	  _h_meff->fill( m_eff, weight);

	  FourMomentum p_llbb = elecs[0].momentum() + muons[0].momentum() + bJets[0].momentum();
	  if (bJets.size() > 1) p_llbb += bJets[1].momentum();
	  double MT = sqrt( pow(p_llbb.mass(),2) + pow(p_llbb.pT(), 2) ) + (-missP.visibleMomentum()).pT();
	  _h_MT->fill( MT, weight);

          // Get 4-momenta of the positively- and negatively-charged leptons
          const FourMomentum& lepPlus = elecs[0].charge()>0 ? elecs[0]:muons[0];
          const FourMomentum& lepMinus = elecs[0].charge()>0 ? muons[0]:elecs[0];

	  double mll = (lepPlus + lepMinus).mass();
	  _h_mll->fill( mll, weight);

          // Now calculate the variable
          double dphi_temp = deltaPhi(lepPlus,lepMinus);
          _h_dphil->fill( dphi_temp, weight );

	  // Evaluate mT and mtt from the 2 sel leps and max 4 candidate jets
#ifdef _SHOW_MASS_CUT_MSGS_
	  cout<<"j="<<goodJets.size()<<" b="<<bJets.size()<<endl;
#endif
	  size_t fnb(bJets.size()); std::vector<const Jet*> vajets;
	  if(fnb>1) for(size_t i=0; i<fnb && i<_NconsidJetsTag_; ++i)
		      vajets.push_back(&bJets[i]);
	  else {
	    vajets.push_back(&bJets[0]);
	    for(size_t i=0;
		i<goodJets.size() && vajets.size()<_NconsidJetsTag_; ++i) {
	      if(goodJets[i].momentum()==bJets[0].momentum()) continue;
	      vajets.push_back(&goodJets[i]);
	    }
	  }
	  assert(vajets.size()>=2 && vajets.size()<=_NconsidJetsTag_);
#ifdef _SHOW_MASS_CUT_MSGS_
	  for(size_t i=0; i<vajets.size(); ++i)
	    cout<<"  "<<i<<" : "<<*vajets[i]
		<<"\t pT="<<vajets[i]->momentum().pT()
		<<"\t isBj="<<vajets[i]->containsBottom()<<endl;
#endif

	  const FourMomentum& wPlus = w1.charge()>0 ? w1:w2;
          const FourMomentum& wMinus = w1.charge()>0 ? w2:w1;
	  std::pair<size_t,size_t> mlbjsel(0,1); double mlbmin(_EcmsTag_);
	  std::pair<size_t,size_t> mttjsel(0,1); double mttmin(_EcmsTag_);
	  //In all cases pairing decision needed to distinguish top from antitop.
	  for(size_t i=0; i<vajets.size(); ++i) {
	    if(fnb==1 && i>0) break;    //Want this one bjet if >2 jets in total.
	    for(size_t j=i+1; j<vajets.size(); ++j) {
	      double mlbA=(lepPlus+vajets[i]->momentum()).mass() +
		(lepMinus+vajets[j]->momentum()).mass();
	      double mlbB=(lepPlus+vajets[j]->momentum()).mass() +
		(lepMinus+vajets[i]->momentum()).mass();
	      if(mlbA<mlbB) { if(mlbA<mlbmin) {
		  mlbmin=mlbA; mlbjsel.first=i; mlbjsel.second=j;}}
	      else { if(mlbB<mlbmin) {
		  mlbmin=mlbB; mlbjsel.first=j; mlbjsel.second=i;}}
	      double mA=(wPlus+vajets[i]->momentum()).mass();
	      mA-=_TopQMass_; if(mA>0.0) mA*=1.1;       //Small penalty for
	      double mtmp=(wMinus+vajets[j]->momentum()).mass();
	      mtmp-=_TopQMass_; if(mtmp>0.0) mtmp*=1.1; //being above
	      mA=fabs(mA)+fabs(mtmp);
	      double mB=(wPlus+vajets[j]->momentum()).mass();
	      mB-=_TopQMass_; if(mB>0.0) mB*=1.1;       //mass shell.
	      mtmp=(wMinus+vajets[i]->momentum()).mass();
	      mtmp-=_TopQMass_; if(mtmp>0.0) mtmp*=1.1;
	      mB=fabs(mB)+fabs(mtmp);
#ifdef _SHOW_MASS_CUT_MSGS_
	      cout<<"      "<<i<<j<<": "
		  <<vajets[i]->momentum().pT()<<"/"<<vajets[j]->momentum().pT()
		  <<"\t "<<mlbA<<" , "<<mlbB<<"\t ::\t "<<mA<<" , "<<mB<<endl;
#endif
	      if(mA<mB) { if(mA<mttmin) {
		  mttmin=mA; mttjsel.first=i; mttjsel.second=j;}}
	      else { if(mB<mttmin) {
		  mttmin=mB; mttjsel.first=j; mttjsel.second=i;}}
	    }
	  }
#ifdef _SHOW_MASS_CUT_MSGS_
	  cout<<"    win(mT)="<<mlbjsel.first<<mlbjsel.second<<" has "<<mlbmin
	      <<"    || "
	      <<"   win(mtt)="<<mttjsel.first<<mttjsel.second<<" has "<<mttmin;
	  cout<<endl;
	  cout<<"     e+:"<<lepPlus<<endl;
	  cout<<"     e-:"<<lepMinus<<endl;
	  cout<<"   visP:"<<missP.visibleMomentum()
	      <<"  Et="<<missP.vectorEt()<<endl;
#endif
	  //Vec4D vis(vectors[0]); for(int i(1); i<n-2; ++i) vis+=vectors[i];
	  //Always pair `Plus' with 1st sel jet and `Minus' with 2nd one.
	  FourMomentum vis(lepPlus); vis+=vajets[mlbjsel.first]->momentum();
	  FourMomentum vit(lepMinus); vit+=vajets[mlbjsel.second]->momentum();
	  double pTlb=vis.pT()+vit.pT(); pTlb/=2.0;
	  vis+=vit;
	  //Vec4D ivis(vectors[n-1]); ivis+=vectors[n-2];
	  FourMomentum ivis(-missP.visibleMomentum());
	  mlbmin=vis.mass2();
	  mlbmin+=2*(sqrt(mlbmin+vis[1]*vis[1]+vis[2]*vis[2]) *
		     sqrt(ivis[1]*ivis[1]+ivis[2]*ivis[2]) -
		     vis[1]*ivis[1]-vis[2]*ivis[2]);
	  mlbmin=sqrt(mlbmin);
	  //double mass(vis.Abs2());
	  //mass+=2*(sqrt(mass+vis[1]*vis[1]+vis[2]*vis[2]) *
	  //sqrt(ivis[1]*ivis[1]+ivis[2]*ivis[2]) -
	  //vis[1]*ivis[1]-vis[2]*ivis[2]);
	  FourMomentum nunu(n1.momentum()+n2.momentum());
	  mttmin=vis.mass2();    //all transverse masses based on mlb j sel
	  mttmin+=2*(sqrt(mttmin+vis[1]*vis[1]+vis[2]*vis[2]) *
		     sqrt(nunu[1]*nunu[1]+nunu[2]*nunu[2]) -
		     vis[1]*nunu[1]-vis[2]*nunu[2]);
	  mttmin=sqrt(mttmin);
	  mll=sqrt(vis.mass2()+vis.pT()*vis.pT());    //mll -> for MT evaluation
	  double MTp=mll+vis.pT(); mll+=ivis.pT();    //reevaluate MT with j sel
	  vit=wPlus; vit+=vajets[mttjsel.first]->momentum();
	  double pTt=vit.pT();
	  vit=wMinus; vit+=vajets[mttjsel.second]->momentum();
	  pTt+=vit.pT(); pTt/=2.0;
	  double mtt=(wPlus+vajets[mttjsel.first]->momentum()+vit).mass();
//	  std::cout << "top mom. = " << wPlus+vajets[mttjsel.first]->momentum() << std::endl;
//	  std::cout << "antitop mom. = " << vit << std::endl;
	  double mttT=(vis+nunu).mass();
#ifdef _SHOW_MASS_CUT_MSGS_
	  cout<<"   mttT="<<mttT<<"\t mtt="<<mtt
	      <<"\t mTnn="<<mttmin<<"\t mT="<<mlbmin
	      <<"\t\t MT="<<mll<<"\t MTp="<<MTp<<"\t [MTLudo="<<MT<<"]"<<endl
	      <<"       ="<<(wPlus+wMinus+vajets[mlbjsel.first]->momentum()+
			     vajets[mlbjsel.second]->momentum()).mass()
	      <<"\t\t\t pTnn="<<nunu.pT()<<"\t pTmiss="<<ivis.pT()
	      <<"\t\t pTllbb="<<vis.pT()<<"\t\t\t [pTllbbLudo="<<p_llbb.pT()<<"]"
	      <<endl<<endl;
#endif


	  m_histfracs["DRatio_mT:mtt"]->fill(mlbmin/mtt, weight);
	  m_histfracs["DRatio_mttT:mtt"]->fill(mttT/mtt, weight);
	  m_histfracs["DRatio_mTnn:mttT"]->fill(mttmin/mttT, weight);
	  m_histfracs["DRatio_mTnn:mT"]->fill(mttmin/mlbmin, weight);
	  m_histfracs["DRatio_MT:mT"]->fill(mll/mlbmin, weight);
	  m_histfracs["DRatio_MTp:mT"]->fill(MTp/mlbmin, weight);
	  m_histfracs["DRatio_pTnn:pTmiss"]->fill(nunu.pT()/ivis.pT(), weight);
	  m_histfracs["DRatio_pTllbb:pTmiss"]->fill(vis.pT()/ivis.pT(), weight);

          // Inclusive definitions of MT
          FourMomentum allmom = -ivis + vis;
          double mT0j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
          std::vector<std::pair<int,string> > i0, i1, i2, i3, i4, ioo;
          unsigned int index = 1;
          double mT1j, mT2j, mT3j, mT4j, mTooj;
          // Not efficient, but let's keep it like that for now
          foreach (const Jet& jet, goodJets) {
              allmom += jet.momentum();
              // std::cout << "allmom " << index << " = " << allmom << std::endl;
              if (index == 1) { mT1j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                mT2j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                mT3j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                mT4j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));}
              else if (index == 2) { mT2j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                     mT3j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                     mT4j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));}
              else if (index == 3) { mT3j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
                                     mT4j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));}
              else if (index == 4) mT4j = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
              ++index;
          }
          mTooj = sqrt(pow(allmom.E(),2) - pow(allmom.pz(),2));
          //std::cout << "MT0j = " << mT0j
          //          << "MT1j = " << mT1j
          //          << "MT2j = " << mT2j
          //          << "MT3j = " << mT3j
          //          << "MT4j = " << mT4j
          //          << "MTooj = " << mTooj << std::endl;
          if (mT0j != 0) m_histfracs["DRatio_mT:mT0j"]->fill(mlbmin/mT0j,weight);
          if (mT1j != 0) m_histfracs["DRatio_mT:mT1j"]->fill(mlbmin/mT1j,weight);
          if (mT2j != 0) m_histfracs["DRatio_mT:mT2j"]->fill(mlbmin/mT2j,weight);
          if (mT3j != 0) m_histfracs["DRatio_mT:mT3j"]->fill(mlbmin/mT3j,weight);
          if (mT4j != 0) m_histfracs["DRatio_mT:mT4j"]->fill(mlbmin/mT4j,weight);
          if (mTooj != 0) m_histfracs["DRatio_mT:mTooj"]->fill(mlbmin/mTooj,weight);
          _h_MT0j->fill(mT0j,weight);
          _h_MT1j->fill(mT1j,weight);
          _h_MT2j->fill(mT2j,weight);
          _h_MT3j->fill(mT3j,weight);
          _h_MT4j->fill(mT4j,weight);
          _h_MTooj->fill(mTooj,weight);


	  // PSEUDO-TOP histograms

//	  const PseudoTop& ttbar = apply<PseudoTop>(event, "ttbar");
//
//          const FourMomentum& t1P4 = ttbar.t1().momentum();
//          const FourMomentum& t2P4 = ttbar.t2().momentum();
//          const double pt1 = std::max(t1P4.pT(), t2P4.pT());
//          const double pt2 = std::min(t1P4.pT(), t2P4.pT());
//          const FourMomentum ttP4 = t1P4 + t2P4;
//
//	  double mtt_ps = ttP4.mass();
//	  double pTt_ps = (t1P4.pT() + t2P4.pT())/2.0;
          double mtt_ps = 0.;
          double pTt_ps = 0.;
//
//	  _h_topPt->fill(pTt_ps, weight);
//	  _h_topPtLead->fill(pt1, weight);
//	  _h_topPtSubLead->fill(pt2, weight);
//	  _h_tMass->fill((t1P4.mass() + t2P4.mass())/2.0, weight);
//	  _h_ttbarMass->fill(mtt_ps, weight);
//
//	  m_histfracs["DRatio_mtt:mtt_ps"]->fill(mtt/mtt_ps, weight);
//	  m_histfracs["DRatio_pTt:pTt_ps"]->fill(pTt/pTt_ps, weight);
//
	  // Neutrino-weighted top reconstruction

	  double m_weight_max = -99.;
	  if ( bJets.size() > 1 ) {

	  NeutrinoWeighter nuW = NeutrinoWeighter(1, lepPlus.pt() + lepPlus.phi());
	  NeutrinoWeighter nuW_alt = NeutrinoWeighter(1, lepPlus.pt() + lepPlus.eta() + bJets.size());//2

	  Vector3 metV = missP.vectorEt();
	  m_weight_max  = nuW.Reconstruct(lepPlus, lepMinus, bJets[0], bJets[1], -metV.x(), -metV.y(), -metV.phi());
          double m_weight_max_alt = nuW_alt.Reconstruct(lepPlus, lepMinus, bJets[1], bJets[0], -metV.x(), -metV.y(), -metV.phi());

	  //std::cout << "weight = " << m_weight_max << " ... " << m_weight_max_alt << std::endl;
          if( m_weight_max_alt > m_weight_max){
            nuW = nuW_alt;
            m_weight_max = m_weight_max_alt;
          }

	  if(m_weight_max > 0.){

            FourMomentum top   = nuW.GetTop();
            FourMomentum tbar  = nuW.GetTbar();
            FourMomentum ttbar = nuW.GetTtbar();
            FourMomentum b     = nuW.GetB();
            FourMomentum bbar  = nuW.GetBbar();
            FourMomentum nu    = nuW.GetNu();
            FourMomentum nubar = nuW.GetNubar();

	    _h_topreco->fill(0.5, 1.);
//	    std::cout << "vw_top mom. = " << top << std::endl;
//	    std::cout << "vw_antitop mom. = " << tbar << std::endl;

	    const double ptlead = std::max(top.pT(), tbar.pT());
	    const double ptsub  = std::min(top.pT(), tbar.pT());
	    double mtt_nw = ttbar.mass();
	    double pTt_nw = (top.pT() + tbar.pT())/2.0;

	    _h_v_topPt->fill(pTt_nw, weight);
	    _h_v_topPtLead->fill(ptlead, weight);
	    _h_v_topPtSubLead->fill(ptsub, weight);
	    _h_v_tMass->fill((top.mass() + tbar.mass())/2.0, weight);
	    _h_v_ttbarMass->fill(mtt_nw, weight);
	    _h_v_ttbarMassATLAS->fill(mtt_nw, weight);
	    _h_v_ttbarMassATLAS2->fill(mtt_nw, weight);

	    m_histfracs["DRatio_mtt:mtt_nw"]->fill(mtt/mtt_nw, weight);
	    m_histfracs["DRatio_pTt:pTt_nw"]->fill(pTt/pTt_nw, weight);

            double mtt1 = 1.0;
            if (bJets.size() >= 2) {

            // int index = -1; int indexb = -1;
             //if bJets[0].containsParticleId(5) { index = 0; indexb = 1; }
           //  else if bJets[0].containsParticleId(-5) { index = 1; indexb = 0; }
            // else std::cout << "First jet does not contain any b-quark." << std::endl;

               mtt1 = (bJets[0].momentum() + wPlus + bJets[1].momentum() + wMinus).mass();
            }

#ifdef _SHOW_MASS_CUT_MSGS_
            std::cout << "mtt nw = " << mtt_nw << ", mtt_truth = " << mtt1 << std::endl;
#endif

            m_histfracs["DRatio_mtt_nw:mtt_truth"]->fill(mtt_nw/mtt1, weight);


            FourMomentum Wp, Wm;
            Wp = lepPlus + nu;
            Wm = lepMinus + nubar;

	    fillMassCutHistos(event, dphi_temp, mtt, mlbmin,
                              mT0j, mT1j, mT2j, mT3j, mT4j, mTooj,
                              pTt, pTlb, mtt_ps, mtt_nw, pTt_ps, pTt_nw);

#ifdef _NW_DEBUG_MSGS_
	    std::cout << " ----------------- 4-momenta -------------------\n"
		      << " **************** neutrino-wgt. ****************\n"
		      << " *    nu: " << nu << "\n"
		      << " * nubar: " << nubar << "\n"
		      << " *  lep+: " << lepPlus << "\n"
		      << " *  lep-: " << lepMinus << "\n"
		      << " ***********************************************\n"
	              << " *    W+: " << Wp << "\n"
		      << " *    W-: " << Wm << "\n"
		      << " *     b: " << b << "\n"
		      << " *  bbar: " << bbar << "\n"
		      << " ***********************************************\n"
		      << " *     t: " << top << "\n"
		      << " *  tbar: " << tbar << std::endl;
#endif

	  }
	  else {
            std::cout << "No mtt_ps matching found!" << std::endl;
	    fillMassCutHistos(event, dphi_temp, mtt, mlbmin,
                              mT0j, mT1j, mT2j, mT3j, mT4j, mTooj,
                              pTt, pTlb, mtt_ps, -99., pTt_ps, -99.);
	  }
	  }
	  else {
	    fillMassCutHistos(event, dphi_temp, mtt, mlbmin,
                              mT0j, mT1j, mT2j, mT3j, mT4j, mTooj,
			      pTt, pTlb, mtt_ps, -99., pTt_ps, -99.);
	  }

	}

      }


    }


    /// Normalise histograms to unit area
    void finalize() {

      const double fctr = crossSection()/sumOfWeights();

      h_cutmtt_xsecs->finalize(); scale(h_cutmtt_xsecs, fctr);
      h_cutmtt_ps_xsecs->finalize(); scale(h_cutmtt_ps_xsecs, fctr);
      h_cutmtt_nw_xsecs->finalize(); scale(h_cutmtt_nw_xsecs, fctr);
      h_cutmT_xsecs->finalize();  scale(h_cutmT_xsecs, fctr);
      h_cutmT0j_xsecs->finalize();  scale(h_cutmT0j_xsecs, fctr);
      h_cutmT1j_xsecs->finalize();  scale(h_cutmT1j_xsecs, fctr);
      h_cutmT2j_xsecs->finalize();  scale(h_cutmT2j_xsecs, fctr);
      h_cutmT3j_xsecs->finalize();  scale(h_cutmT3j_xsecs, fctr);
      h_cutmT4j_xsecs->finalize();  scale(h_cutmT4j_xsecs, fctr);
      h_cutmTooj_xsecs->finalize();  scale(h_cutmTooj_xsecs, fctr);
      h_cutpTt_xsecs->finalize();  scale(h_cutpTt_xsecs, fctr);
      h_cutpTt_ps_xsecs->finalize();  scale(h_cutpTt_ps_xsecs, fctr);
      h_cutpTt_nw_xsecs->finalize();  scale(h_cutpTt_nw_xsecs, fctr);
      h_cutpTlb_xsecs->finalize(); scale(h_cutpTlb_xsecs, fctr);
      for(size_t i=0; i<v_cutmtt.size(); ++i) {
	v_cutmtt_histdphil[i]->finalize(); scale(v_cutmtt_histdphil[i], fctr);
	v_cutmtt_histmtt[i]->finalize();   scale(v_cutmtt_histmtt[i], fctr);
	v_cutmtt_histmT[i]->finalize();    scale(v_cutmtt_histmT[i], fctr);
	v_cutmtt_histpTt[i]->finalize();   scale(v_cutmtt_histpTt[i], fctr);
	v_cutmtt_histpTlb[i]->finalize();  scale(v_cutmtt_histpTlb[i], fctr);
	v_cutmtt_ps_histdphil[i]->finalize(); scale(v_cutmtt_ps_histdphil[i], fctr);
	v_cutmtt_ps_histmtt[i]->finalize();   scale(v_cutmtt_ps_histmtt[i], fctr);
	v_cutmtt_ps_histmT[i]->finalize();    scale(v_cutmtt_ps_histmT[i], fctr);
	v_cutmtt_ps_histpTt[i]->finalize();   scale(v_cutmtt_ps_histpTt[i], fctr);
	v_cutmtt_ps_histpTlb[i]->finalize();  scale(v_cutmtt_ps_histpTlb[i], fctr);
	v_cutmtt_nw_histdphil[i]->finalize(); scale(v_cutmtt_nw_histdphil[i], fctr);
	v_cutmtt_nw_histmtt[i]->finalize();   scale(v_cutmtt_nw_histmtt[i], fctr);
	v_cutmtt_nw_histmT[i]->finalize();    scale(v_cutmtt_nw_histmT[i], fctr);
	v_cutmtt_nw_histpTt[i]->finalize();   scale(v_cutmtt_nw_histpTt[i], fctr);
	v_cutmtt_nw_histpTlb[i]->finalize();  scale(v_cutmtt_nw_histpTlb[i], fctr);
      }
      for(size_t i=0; i<v_cutmT.size(); ++i) {
	v_cutmT_histdphil[i]->finalize(); scale(v_cutmT_histdphil[i], fctr);
	v_cutmT_histmtt[i]->finalize();   scale(v_cutmT_histmtt[i], fctr);
	v_cutmT_histmT[i]->finalize();    scale(v_cutmT_histmT[i], fctr);
	v_cutmT_histpTt[i]->finalize();   scale(v_cutmT_histpTt[i], fctr);
	v_cutmT_histpTlb[i]->finalize();  scale(v_cutmT_histpTlb[i], fctr);
	v_cutmT0j_histdphil[i]->finalize(); scale(v_cutmT0j_histdphil[i], fctr);
	v_cutmT0j_histmT0j[i]->finalize();    scale(v_cutmT0j_histmT0j[i], fctr);
	v_cutmT1j_histdphil[i]->finalize(); scale(v_cutmT1j_histdphil[i], fctr);
	v_cutmT1j_histmT1j[i]->finalize();    scale(v_cutmT1j_histmT1j[i], fctr);
	v_cutmT2j_histdphil[i]->finalize(); scale(v_cutmT2j_histdphil[i], fctr);
	v_cutmT2j_histmT2j[i]->finalize();    scale(v_cutmT2j_histmT2j[i], fctr);
	v_cutmT3j_histdphil[i]->finalize(); scale(v_cutmT3j_histdphil[i], fctr);
	v_cutmT3j_histmT3j[i]->finalize();    scale(v_cutmT3j_histmT3j[i], fctr);
	v_cutmT4j_histdphil[i]->finalize(); scale(v_cutmT4j_histdphil[i], fctr);
	v_cutmT4j_histmT4j[i]->finalize();    scale(v_cutmT4j_histmT4j[i], fctr);
	v_cutmTooj_histdphil[i]->finalize(); scale(v_cutmTooj_histdphil[i], fctr);
	v_cutmTooj_histmTooj[i]->finalize();    scale(v_cutmTooj_histmTooj[i], fctr);
      }
      for(size_t i=0; i<v_cutpTt.size(); ++i) {
	v_cutpTt_histdphil[i]->finalize(); scale(v_cutpTt_histdphil[i], fctr);
	v_cutpTt_histmtt[i]->finalize();   scale(v_cutpTt_histmtt[i], fctr);
	v_cutpTt_histmT[i]->finalize();    scale(v_cutpTt_histmT[i], fctr);
	v_cutpTt_histpTt[i]->finalize();   scale(v_cutpTt_histpTt[i], fctr);
	v_cutpTt_histpTlb[i]->finalize();  scale(v_cutpTt_histpTlb[i], fctr);
	v_cutpTt_ps_histdphil[i]->finalize(); scale(v_cutpTt_ps_histdphil[i], fctr);
	v_cutpTt_ps_histmtt[i]->finalize();   scale(v_cutpTt_ps_histmtt[i], fctr);
	v_cutpTt_ps_histmT[i]->finalize();    scale(v_cutpTt_ps_histmT[i], fctr);
	v_cutpTt_ps_histpTt[i]->finalize();   scale(v_cutpTt_ps_histpTt[i], fctr);
	v_cutpTt_ps_histpTlb[i]->finalize();  scale(v_cutpTt_ps_histpTlb[i], fctr);
	v_cutpTt_nw_histdphil[i]->finalize(); scale(v_cutpTt_nw_histdphil[i], fctr);
	v_cutpTt_nw_histmtt[i]->finalize();   scale(v_cutpTt_nw_histmtt[i], fctr);
	v_cutpTt_nw_histmT[i]->finalize();    scale(v_cutpTt_nw_histmT[i], fctr);
	v_cutpTt_nw_histpTt[i]->finalize();   scale(v_cutpTt_nw_histpTt[i], fctr);
	v_cutpTt_nw_histpTlb[i]->finalize();  scale(v_cutpTt_nw_histpTlb[i], fctr);
      }
      for(size_t i=0; i<v_cutpTlb.size(); ++i) {
	v_cutpTlb_histdphil[i]->finalize(); scale(v_cutpTlb_histdphil[i], fctr);
	v_cutpTlb_histmtt[i]->finalize();   scale(v_cutpTlb_histmtt[i], fctr);
	v_cutpTlb_histmT[i]->finalize();    scale(v_cutpTlb_histmT[i], fctr);
	v_cutpTlb_histpTt[i]->finalize();   scale(v_cutpTlb_histpTt[i], fctr);
	v_cutpTlb_histpTlb[i]->finalize();  scale(v_cutpTlb_histpTlb[i], fctr);
      }
      for(std::map<std::string,NLOHisto1DPtr>::iterator hit=m_histfracs.begin();
          hit!=m_histfracs.end(); ++hit) {
        hit->second->finalize(); scale(hit->second, fctr);
      }

      _h_dphil->finalize();           scale(_h_dphil, fctr);
      _h_mll->finalize();             scale(_h_mll, fctr);
      _h_njets->finalize();           scale(_h_njets, fctr);
      _h_nbjets->finalize();          scale(_h_nbjets, fctr);
      _h_pTb1->finalize();            scale(_h_pTb1, fctr);
      _h_pTb2->finalize();            scale(_h_pTb2, fctr);
      _h_pTj1->finalize();            scale(_h_pTj1, fctr);
      _h_pTj2->finalize();            scale(_h_pTj2, fctr);
      _h_DRj->finalize();             scale(_h_DRj, fctr);
      _h_dphij->finalize();           scale(_h_dphij, fctr);
      _h_mjj->finalize();             scale(_h_mjj, fctr);
      _h_ht->finalize();              scale(_h_ht, fctr);
      _h_ht2->finalize();             scale(_h_ht2, fctr);
      _h_meff->finalize();            scale(_h_meff, fctr);
      _h_MT->finalize();              scale(_h_MT, fctr);
      _h_MT0j->finalize();              scale(_h_MT0j, fctr);
      _h_MT1j->finalize();              scale(_h_MT1j, fctr);
      _h_MT2j->finalize();              scale(_h_MT2j, fctr);
      _h_MT3j->finalize();              scale(_h_MT3j, fctr);
      _h_MT4j->finalize();              scale(_h_MT4j, fctr);
      _h_MTooj->finalize();              scale(_h_MTooj, fctr);
      _h_met->finalize();             scale(_h_met, fctr);
      _h_eta_el->finalize();          scale(_h_eta_el, fctr);
      _h_eta_mu->finalize();          scale(_h_eta_mu, fctr);
      _h_topPt->finalize();           scale(_h_topPt, fctr);
      _h_topPtLead->finalize();       scale(_h_topPtLead, fctr);
      _h_topPtSubLead->finalize();    scale(_h_topPtSubLead, fctr);
      _h_tMass->finalize();           scale(_h_tMass, fctr);
      _h_ttbarMass->finalize();       scale(_h_ttbarMass, fctr);
      _h_v_topPt->finalize();         scale(_h_v_topPt, fctr);
      _h_v_topPtLead->finalize();     scale(_h_v_topPtLead, fctr);
      _h_v_topPtSubLead->finalize();  scale(_h_v_topPtSubLead, fctr);
      _h_v_tMass->finalize();         scale(_h_v_tMass, fctr);
      _h_v_ttbarMass->finalize();     scale(_h_v_ttbarMass, fctr);
      _h_v_ttbarMassATLAS->finalize();scale(_h_v_ttbarMassATLAS, fctr);
      _h_v_ttbarMassATLAS2->finalize();scale(_h_v_ttbarMassATLAS2, fctr);
      _h_selected_events->finalize(); scale(_h_selected_events, fctr);
      _h_xs->finalize();              scale(_h_xs, fctr);
      _h_xscut->finalize();           scale(_h_xscut, fctr);
      scale(_h_topreco, 1./_h_xscut->numEntries()); // numEvent()
      _h_eff->finalize();             scale(_h_eff, fctr);
      scale(_h_eff, 1./_h_xs->integral());

    }

  
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC13TeV_SPINCORR_emu_TOP_tt);


}



namespace Rivet {
#include "WFinder221.cc"
#include "PseudoTop.cc"
#include "NeutrinoWeighter.cc"
}
