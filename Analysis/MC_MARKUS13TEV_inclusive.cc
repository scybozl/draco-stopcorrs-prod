// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
//#include "Rivet/Tools/RivetMT2.hh"

namespace Rivet {


  class MC_MARKUS13TEV_inclusive : public Analysis {

#include "NLOHisto1D.cc"

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_MARKUS13TEV_inclusive()
      : Analysis("MC_MARKUS13TEV_inclusive")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

       VetoedFinalState fs;
       fs.addVetoPairId(PID::ELECTRON);
       fs.addVetoPairId(PID::MUON);
       fs.addVetoPairId(PID::NU_E);
       fs.addVetoPairId(PID::NU_MU);
       addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

       IdentifiedFinalState elecs;
       elecs.acceptIdPair(PID::ELECTRON);
       addProjection(elecs, "elecs");

       IdentifiedFinalState muons;
       muons.acceptIdPair(PID::MUON);
       addProjection(muons, "muons");

       IdentifiedFinalState nue;
       nue.acceptIdPair(PID::NU_E);
       addProjection(nue, "nue");

       IdentifiedFinalState numu;
       numu.acceptIdPair(PID::NU_MU);
       addProjection(numu, "numu");

       addProjection(VisibleFinalState(-5, 5),"vfs");


       _h_xs = bookNLOHisto1D("xs", 1, 0.5, 1.5);
       _h_xscut = bookNLOHisto1D("xscut", 1, 0.5, 1.5);
       _h_pTbmin = bookNLOHisto1D("pTbmin", 50, 20., 270.);
       _h_pTbmax = bookNLOHisto1D("pTbmax", 40, 0., 800.);
       _h_pTmu = bookNLOHisto1D("pTmu", 40, 0., 1000.);
       _h_etamu = bookNLOHisto1D("etamu", 50, -5, 5);
       _h_mll = bookNLOHisto1D("mll", 500, 0., 1000.);
       _h_mlb = bookNLOHisto1D("mlb", 100, 0., 200.);
       _h_mlb_diff = bookNLOHisto1D("mlbdiff", 100, -100., 100.);
       _h_mwb = bookNLOHisto1D("mwb", 100, 150., 200.);
       _h_etdr = bookNLOHisto1D("etdr", 200, 0., 400.);
       _h_malt = bookNLOHisto1D("malt", 400, 0., 800.); //not lorentz invariant
       _h_pTlb = bookNLOHisto1D("pTlb", 400, 0., 800.);
       _h_pTl = bookNLOHisto1D("pTl", 400, 0., 800.); //two entries per event
       _h_pTb = bookNLOHisto1D("pTb", 400, 0., 800.); //two entries per event
       _h_selected_events = bookNLOHisto1D("h_selected_events", 5, 0.0, 5.0);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
//      const double weight = event.weight();

        _h_xs->fill(1,event);
        _h_selected_events->fill(0.5,event);
//beginning of cuts ---------------------------------------
	// TODO: Change for hadron level running. More than one electron.
       ParticleVector elec = applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt();
	if (elec.size() !=1){
		vetoEvent;
	}
	if (elec[0].momentum().pT()<28 || fabs(elec[0].momentum().eta())>2.47 || (fabs(elec[0].momentum().eta())>1.37 && fabs(elec[0].momentum().eta())<1.52)){
//	if (elec[0].momentum().pT()<28 || fabs(elec[0].momentum().eta())>2.47 ){
		vetoEvent;
	}
       _h_selected_events->fill(1.5,event);
       ParticleVector muon = applyProjection<IdentifiedFinalState>(event, "muons").particlesByPt();
	if (muon.size() !=1){
		vetoEvent;
	}
	if (muon[0].momentum().pT()<28 || fabs(muon[0].momentum().eta())>2.5){
		vetoEvent;
	}
       _h_selected_events->fill(2.5,event);

	Jets goodjets;
	Jets bjets;
	const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
	const Jets alljets = jetpro.jetsByPt();
//cuts for ttbarjet-----------------------------------------
//if (jetpro.jetsByPt(15.*GeV).size()<3){
//	vetoEvent;
//}
//cuts for ttbarjet-----------------------------------------

	double jetpt = 0;
	foreach (const Jet& jet, alljets) {
		if ( fabs( jet.momentum().eta() ) < 2.5 && jet.momentum().pT() > 25  ) {
			//if  (deltaR(jet.momentum(), elec[0].momentum())>0.4 && deltaR(jet.momentum(), muon[0].momentum())>0.4) {
			jetpt += jet.momentum().pT(); //for H_T cut
			goodjets.push_back(jet);
			if (jet.containsBottom() ) {
				bjets.push_back(jet);
			}
			//}
		}
	}

	foreach(const Jet& jet, goodjets){
		if(deltaR(jet.momentum(), elec[0].momentum())<0.4 || deltaR(jet.momentum(), muon[0].momentum())<0.4){
			vetoEvent;
		}
	}

	if (bjets.size() < 2) {
		vetoEvent;
	} 
        _h_selected_events->fill(3.5,event);
       

//plb cut
	FourMomentum p00 = elec[0].momentum() + bjets[0].momentum();
	FourMomentum p10 = muon[0].momentum() + bjets[0].momentum();
	FourMomentum p01 = elec[0].momentum() + bjets[1].momentum();
	FourMomentum p11 = muon[0].momentum() + bjets[1].momentum();

//	if((p00.mass()+p11.mass())<(p01.mass()+p10.mass())){
//		if((p00.pT()+p11.pT())/2<120){
//			vetoEvent;
//		}
//	} else {
//		if((p01.pT()+p10.pT())/2<120){
//			vetoEvent;
//		}
//	}
        _h_selected_events->fill(4.5,event);

//end of cuts ----------------------------------------------

//end of cuts ----------------------------------------------
//if(bjets.size()>2){
//cout << bjets.size() << "\n";
//cout << event.genEvent().event_number() << "\n";
//vetoEvent;
//}


	ParticleVector vfs_particles = applyProjection<VisibleFinalState>(event, "vfs").particles();
	FourMomentum pTmiss;
	foreach ( const Particle & p, vfs_particles ) {
		pTmiss -= p.momentum();
	}

	if( bjets[0].momentum().pT() > bjets[1].momentum().pT() ){
		_h_pTbmax->fill(bjets[0].momentum().pT(), event);
		_h_pTbmin->fill(bjets[1].momentum().pT(), event);
	}else{
		_h_pTbmax->fill(bjets[1].momentum().pT(), event);
		_h_pTbmin->fill(bjets[0].momentum().pT(), event);
	}
	_h_pTmu->fill(muon[0].momentum().pT(), event);
	_h_etamu->fill(muon[0].momentum().eta(), event);
	_h_mll->fill((muon[0].momentum()+elec[0].momentum()).mass(), event);

	_h_pTl->fill(elec[0].momentum().pT(), event);
	_h_pTl->fill(muon[0].momentum().pT(), event);
	_h_pTb->fill(bjets[0].momentum().pT(), event);
	_h_pTb->fill(bjets[1].momentum().pT(), event);

        ParticleVector nue = applyProjection<IdentifiedFinalState>(event, "nue").particlesByPt();
        ParticleVector numu = applyProjection<IdentifiedFinalState>(event, "numu").particlesByPt();

	double mlb;
	FourMomentum pwb1;
	FourMomentum pwb2;
	double mwb;
	double etdr;
	double m_T2;
	if((p00.mass()+p11.mass())<(p01.mass()+p10.mass())){
		mlb =  (p00.mass() + p11.mass())/2; 
		pwb1 = p00+nue[0].momentum();
		pwb2 = p11+numu[0].momentum();
		mwb =  (pwb1.mass() + pwb2.mass())/2;
		double pa[3]    = {p00.mass() , p00.x() , p00.y() };
		double pb[3]    = {p11.mass() , p11.x() , p11.y() };
		double pmiss[3] = {-999.999 , pTmiss.x() , pTmiss.y() };
		etdr=(elec[0].momentum().Et()*deltaR(elec[0].momentum(),bjets[0].momentum()) + muon[0].momentum().Et()*deltaR(muon[0].momentum(),bjets[1].momentum()))/2;
		_h_pTlb->fill((p00.pT()+p11.pT())/2, event);
	} else {
		mlb =  (p10.mass() + p01.mass())/2;
		pwb1 = p10+numu[0].momentum();
		pwb2 = p01+nue[0].momentum();
		mwb =  (pwb1.mass() + pwb2.mass())/2;
		double pa[3]    = {p10.mass() , p10.x() , p10.y() };
		double pb[3]    = {p01.mass() , p01.x() , p01.y() };
		double pmiss[3] = {-999.999 , pTmiss.x() , pTmiss.y() };
		etdr=(elec[0].momentum().Et()*deltaR(elec[0].momentum(),bjets[1].momentum()) + muon[0].momentum().Et()*deltaR(muon[0].momentum(),bjets[0].momentum()))/2;
		_h_pTlb->fill((p01.pT()+p10.pT())/2, event);
	}
	double mlb_diff = (p00.mass()+p11.mass())-(p01.mass()+p10.mass());
	_h_mlb_diff->fill(mlb_diff, event);
	_h_mlb->fill(mlb , event);
	_h_mwb->fill(mwb , event);
	_h_etdr->fill(etdr , event);

	double malt;
	if(elec[0].momentum().vector3().unit().dot(bjets[0].momentum().vector3().unit())+muon[0].momentum().vector3().unit().dot(bjets[1].momentum().vector3().unit()) > elec[0].momentum().vector3().unit().dot(bjets[1].momentum().vector3().unit())+muon[0].momentum().vector3().unit().dot(bjets[0].momentum().vector3().unit())){
		malt=sqrt((elec[0].momentum().E()-elec[0].momentum().vector3().dot(bjets[0].momentum().vector3().unit()))*(muon[0].momentum().E()-muon[0].momentum().vector3().dot(bjets[1].momentum().vector3().unit())));
	} else {
		malt=sqrt((elec[0].momentum().E()-elec[0].momentum().vector3().dot(bjets[1].momentum().vector3().unit()))*(muon[0].momentum().E()-muon[0].momentum().vector3().dot(bjets[0].momentum().vector3().unit())));
	}
	_h_malt->fill(malt,event);
        _h_xscut->fill(1,event);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

double factor = crossSection()/sumOfWeights();

//test
//cout << "factor: " << factor <<"\n";
//cout << "xs: " << crossSection() <<"\n";
//cout << "sumofweights: " << sumOfWeights() <<"\n";
//test
 _h_xs->finalize();
scale(_h_xs, factor);
 _h_xscut->finalize();
scale(_h_xscut, factor);
 _h_pTbmin->finalize();
scale(_h_pTbmin, factor);
 _h_pTbmax->finalize();
scale(_h_pTbmax, factor);
 _h_pTmu->finalize();
scale(_h_pTmu, factor);
 _h_etamu->finalize();
scale(_h_etamu, factor);
 _h_mll->finalize();
scale(_h_mll, factor);
 _h_mlb->finalize();
scale(_h_mlb, factor);
 _h_mlb_diff->finalize();
scale(_h_mlb_diff, factor);
 _h_mwb->finalize();
scale(_h_mwb, factor);
 _h_etdr->finalize();
scale(_h_etdr, factor);
 _h_malt->finalize();
scale(_h_malt, factor);
 _h_pTlb->finalize();
scale(_h_pTlb, factor);
 _h_pTl->finalize();
scale(_h_pTl, factor);
 _h_pTb->finalize();
scale(_h_pTb, factor);
 _h_selected_events->finalize();
scale(_h_selected_events, factor);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
  NLOHisto1DPtr  _h_xs;
  NLOHisto1DPtr  _h_xscut;
  NLOHisto1DPtr  _h_pTbmin;
  NLOHisto1DPtr  _h_pTbmax;
  NLOHisto1DPtr  _h_pTmu;
  NLOHisto1DPtr  _h_etamu;
  NLOHisto1DPtr  _h_mll;
  NLOHisto1DPtr  _h_mlb;
  NLOHisto1DPtr  _h_mlb_diff;
  NLOHisto1DPtr  _h_mwb;
  NLOHisto1DPtr  _h_etdr;
  NLOHisto1DPtr  _h_malt;
  NLOHisto1DPtr  _h_pTlb;
  NLOHisto1DPtr  _h_pTl;
  NLOHisto1DPtr  _h_pTb;
  NLOHisto1DPtr  _h_selected_events;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_MARKUS13TEV_inclusive);

}
