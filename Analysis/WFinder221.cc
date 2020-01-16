// -*- C++ -*-
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"


WFinder221::WFinder221(const FinalState& inputfs,
		       const Cut& fsCut,
		       PdgId pid,
		       double minmass, double maxmass,
		       double missingET,
		       double dRmax,
		       ClusterPhotons clusterPhotons,
		       PhotonTracking trackPhotons,
		       MassWindow masstype,
		       double masstarget) {
  setName("WFinder221");

  _minmass = minmass;
  _maxmass = maxmass;
  _masstarget = masstarget;
  _pid = pid;
  _trackPhotons = trackPhotons;
  _useTransverseMass = (masstype == TRANSMASS);

  // Check that the arguments are legal
  assert(abs(_pid) == PID::ELECTRON || abs(_pid) == PID::MUON);
  _nu_pid = abs(_pid) + 1;
  assert(abs(_nu_pid) == PID::NU_E || abs(_nu_pid) == PID::NU_MU);

  // Don't make pT or eta cuts on the neutrino
  IdentifiedFinalState neutrinos(inputfs);
  neutrinos.acceptNeutrinos();
  addProjection(neutrinos, "Neutrinos");

  // Lepton clusters
  IdentifiedFinalState bareleptons(inputfs);
  bareleptons.acceptIdPair(pid);
  const bool doClustering = (clusterPhotons != NOCLUSTER);
  const bool useDecayPhotons = (clusterPhotons == CLUSTERALL);
  DressedLeptons leptons(inputfs, bareleptons, dRmax, fsCut, doClustering, useDecayPhotons);
  addProjection(leptons, "DressedLeptons");

  // Add MissingMomentum proj to calc MET
  MissingMomentum vismom(inputfs);
  addProjection(vismom, "MissingET");
  // Set ETmiss
  _etMiss = missingET;

  VetoedFinalState remainingFS;
  remainingFS.addVetoOnThisFinalState(*this);
  addProjection(remainingFS, "RFS");
}


/////////////////////////////////////////////////////


const FinalState& WFinder221::remainingFinalState() const {
  return getProjection<FinalState>("RFS");
}


int WFinder221::compare(const Projection& p) const {
  PCmp dlcmp = mkNamedPCmp(p, "DressedLeptons");
  if (dlcmp != EQUIVALENT) return dlcmp;

  const WFinder221& other = dynamic_cast<const WFinder221&>(p);
  return (cmp(_minmass, other._minmass) || cmp(_maxmass, other._maxmass) ||
	  cmp(_useTransverseMass, other._useTransverseMass) ||
	  cmp(_etMiss, other._etMiss) ||
	  cmp(_pid, other._pid) || cmp(_trackPhotons, other._trackPhotons));
}


void WFinder221::project(const Event& e) {
  clear();

  const DressedLeptons& leptons = applyProjection<DressedLeptons>(e, "DressedLeptons");
  const FinalState& neutrinos = applyProjection<FinalState>(e, "Neutrinos");

  // Make and register an invariant mass final state for the W decay leptons
  vector<pair<PdgId, PdgId> > l_nu_ids;
  l_nu_ids += make_pair(abs(_pid), -abs(_nu_pid));
  l_nu_ids += make_pair(-abs(_pid), abs(_nu_pid));
  InvMassFinalState imfs(l_nu_ids, _minmass, _maxmass, _masstarget);
  imfs.useTransverseMass(_useTransverseMass);
  Particles tmp;
  tmp.insert(tmp.end(), leptons.dressedLeptons().begin(), leptons.dressedLeptons().end());
  tmp.insert(tmp.end(), neutrinos.particles().begin(), neutrinos.particles().end());
  imfs.calc(tmp);

  if (imfs.particlePairs().size() < 1) return;

  ParticlePair Wconstituents(imfs.particlePairs()[0]);
  Particle p1(Wconstituents.first), p2(Wconstituents.second);

  if (threeCharge(p1) == 0) {
    _constituentLeptons += p2;
    _constituentNeutrinos += p1;
  } else {
    _constituentLeptons += p1;
    _constituentNeutrinos += p2;
  }

  FourMomentum pW = p1.momentum() + p2.momentum();
  const int w3charge = threeCharge(p1) + threeCharge(p2);
  assert(abs(w3charge) == 3);
  const int wcharge = w3charge/3;

  stringstream msg;
  string wsign = (wcharge == 1) ? "+" : "-";
  string wstr = "W" + wsign;
  msg << wstr << " reconstructed from: " << "\n"
      << "   " << p1.momentum() << " " << p1.pid() << "\n"
      << " + " << p2.momentum() << " " << p2.pid();
  MSG_DEBUG(msg.str());

  // Check missing ET
  const MissingMomentum& vismom = applyProjection<MissingMomentum>(e, "MissingET");
  const double met = vismom.vectorEt().mod();
  MSG_TRACE("MET = " << met/GeV << " GeV vs. required = " << _etMiss/GeV << " GeV");
  if (met < _etMiss) {
    MSG_DEBUG("Not enough missing ET: " << met/GeV << " GeV vs. " << _etMiss/GeV << " GeV");
    return;
  }

  // Make W Particle and insert into particles list
  const PdgId wpid = (wcharge == 1) ? PID::WPLUSBOSON : PID::WMINUSBOSON;
  _bosons.push_back(Particle(wpid, pW));

  // Find the DressedLeptons and neutrinos which survived the IMFS cut such that we can
  // extract their original particles
  foreach (const Particle& p, _constituentNeutrinos) {
    _theParticles.push_back(p);
  }
  foreach (const Particle& p, _constituentLeptons) {
    foreach (const DressedLepton& l, leptons.dressedLeptons()) {
      if (p.pid() == l.pid() && p.momentum() == l.momentum()) {
	_theParticles.push_back(l.constituentLepton());
	if (_trackPhotons) {
	  _theParticles.insert(_theParticles.end(), l.constituentPhotons().begin(), l.constituentPhotons().end());
	}
      }
    }
  }
}
