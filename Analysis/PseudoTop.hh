// -*- C++ -*-
#ifndef RIVET_PseudoTop_HH
#define RIVET_PseudoTop_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

    // Trace event record to see if particle came from a hadron (or a tau from a hadron decay)
    // Based on fromDecay() function
    bool fromHadronDecay( const Particle& p ) {
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == NULL) return false;
      foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
        const PdgId pid = ancestor->pdg_id();
        if (ancestor->status() == 2 && PID::isHadron(pid)) return true;
        if (ancestor->status() == 2 && (abs(pid) == PID::TAU && fromHadronDecay(ancestor))) return true;
      }
      return false;
    }

    /// Clone on the heap.
  /// @brief Convenience finder of leptonically decaying Ws
  ///
  /// Chain together different projections as convenience for finding W's
  /// from two leptons in the final state, including photon clustering.
  ///
  /// @todo Inherit directly from ParticleFinder, not FinalState
  class PseudoTop : public FinalState {
  public:

    /// @name Constructors
    //@{
    PseudoTop(double lepR = 0.1, double lepMinPt = 20, double lepMaxEta = 2.4,
              double jetR = 0.4, double jetMinPt = 30, double jetMaxEta = 4.7)
      : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
        _lepR(lepR), _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta),
        _jetR(jetR), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta)
    {
      setName("PseudoTop");
    }

    enum TTbarMode {CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON};
    enum DecayMode {CH_HADRON = 0, CH_MUON, CH_ELECTRON};

    TTbarMode mode() const {
      if (!_isValid) return CH_NONE;
      if (_mode1 == CH_HADRON && _mode2 == CH_HADRON) return CH_FULLHADRON;
      else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON) return CH_FULLLEPTON;
      else return CH_SEMILEPTON;
    }
    DecayMode mode1() const {return _mode1;}
    DecayMode mode2() const {return _mode2;}

    virtual unique_ptr<Projection> clone() const {
      return unique_ptr<Projection>(new PseudoTop(*this));
    }


    public:
      Particle t1() const {return _t1;}
      Particle t2() const {return _t2;}
      Particle b1() const {return _b1;}
      Particle b2() const {return _b2;}
      ParticleVector wDecays1() const {return _wDecays1;}
      ParticleVector wDecays2() const {return _wDecays2;}
      Jets jets() const {return _jets;}
      Jets bjets() const {return _bjets;}
      Jets ljets() const {return _ljets;}

    protected:
      // Apply the projection to the event
      void project(const Event& e); // override; ///< @todo Re-enable when C++11 allowed
      void cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup=false) const;

    private:
      const double _lepR, _lepMinPt, _lepMaxEta;
      const double _jetR, _jetMinPt, _jetMaxEta;

      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _tMass; // = 172.5*GeV; ///< @todo Re-enable when C++11 allowed
      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _wMass; // = 80.4*GeV; ///< @todo Re-enable when C++11 allowed

    private:
      bool _isValid;
      DecayMode _mode1, _mode2;

      Particle _t1, _t2;
      Particle _b1, _b2;
      ParticleVector _wDecays1, _wDecays2;
      Jets _jets, _bjets, _ljets;


  };


}


#endif
