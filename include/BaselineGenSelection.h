#pragma once

#include "AnalysisTools.h"
#include "GenTools.h"
#include "TextIO.h"
#include "HHCore.h"


void AssignHadronicVCand(VCand& cand, int cand_idx, const RVecI& GenPart_pdgId, const RVecVecI& GenPart_daughters,
                         const RVecF& GenPart_pt, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenPart_mass)
{
  cand.leg_kind[0] = Vleg::Jet;
  cand.leg_kind[1] = Vleg::Jet;

  RVecI daughters = GenPart_daughters.at(cand_idx);
  int first = GenPart_pt[daughters[0]] > GenPart_pt[daughters[1]] ? 0 : 1;
  int second = GenPart_pt[daughters[0]] > GenPart_pt[daughters[1]] ? 1 : 0;
  cand.leg_index[0] = daughters[first];
  cand.leg_index[1] = daughters[second];

  cand.cand_p4 = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, cand_idx);

  cand.leg_p4[0] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, cand.leg_index[0]); // use get mass function
  cand.leg_p4[1] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, cand.leg_index[1]); // fix later

  cand.leg_vis_p4[0] = LorentzVectorM{}; // fix later
  cand.leg_vis_p4[1] = LorentzVectorM{};
}

void AssignLeptonicVCand(VCand& cand, int cand_idx, std::vector<reco_tau::gen_truth::GenLepton> const& gen_leptons,
                         const RVecI& GenPart_pdgId, const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags,
                         const RVecF& GenPart_pt, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenPart_mass)
{
  cand.cand_p4 = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, cand_idx);

  RVecI daughters = GenPart_daughters.at(cand_idx);
  std::sort(daughters.begin(), daughters.end(), [&GenPart_pdgId](int x, int y){ return std::abs(GenPart_pdgId[x]) < std::abs(GenPart_pdgId[y]); });

  bool invisible = PdG::isNeutrino(GenPart_pdgId[daughters[0]]) && PdG::isNeutrino(GenPart_pdgId[daughters[1]]);
  bool lep_pair = PdG::isLepton(GenPart_pdgId[daughters[0]]) && PdG::isLepton(GenPart_pdgId[daughters[1]]);

  if (invisible)
  {
    cand.leg_index[0] = daughters[0];
    cand.leg_kind[0] = Vleg::Nu;
    cand.leg_p4[0] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, daughters[0]);
    cand.leg_vis_p4[0] = LorentzVectorM(0, 0, 0, 0);

    cand.leg_index[1] = daughters[1];
    cand.leg_kind[1] = Vleg::Nu;
    cand.leg_p4[1] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass,  daughters[1]);
    cand.leg_vis_p4[1] = LorentzVectorM(0, 0, 0, 0);

    return;
  }

  if (lep_pair)
  {
    reco_tau::gen_truth::GenLepton const* l1_ptr = findLeptonByIndex(gen_leptons, daughters[0]);
    if (!l1_ptr)
    {
      throw analysis::exception("Could not find first gen lepton by index");
    }

    cand.leg_index[0] = daughters[0];
    cand.leg_kind[0] = static_cast<Vleg>(l1_ptr->kind());
    cand.leg_p4[0] = l1_ptr->lastCopy().p4;
    cand.leg_vis_p4[0] = LorentzVectorM(l1_ptr->visibleP4());

    reco_tau::gen_truth::GenLepton const* l2_ptr = findLeptonByIndex(gen_leptons, daughters[1]);
    if (!l2_ptr)
    {
      throw analysis::exception("Could not find second gen lepton by index");
    }

    cand.leg_index[1] = daughters[1];
    cand.leg_kind[1] = static_cast<Vleg>(l2_ptr->kind());
    cand.leg_p4[1] = l2_ptr->lastCopy().p4;
    cand.leg_vis_p4[1] = LorentzVectorM(l2_ptr->visibleP4());

    return;
  }

  int nu_idx = daughters.at(1);
  int lep_idx = daughters.at(0);

  reco_tau::gen_truth::GenLepton const* lep_ptr = findLeptonByIndex(gen_leptons, lep_idx);
  if (!lep_ptr)
  {
    throw analysis::exception("Could not find gen lepton by index");
  }

  cand.leg_index[1] = nu_idx;
  cand.leg_kind[1] = Vleg::Nu;
  cand.leg_p4[1] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, nu_idx);
  cand.leg_vis_p4[1] = LorentzVectorM(0, 0, 0, 0);

  cand.leg_index[0] = lep_idx;
  cand.leg_kind[0] = static_cast<Vleg>(lep_ptr->kind());
  cand.leg_p4[0] = lep_ptr->lastCopy().p4;
  cand.leg_vis_p4[0] = LorentzVectorM(lep_ptr->visibleP4());
}


VCand GetGenVCand(int evt, int v_idx, std::vector<reco_tau::gen_truth::GenLepton> const& gen_leptons,
                  const RVecI& GenPart_pdgId, const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags,
                  const RVecF& GenPart_pt, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenPart_mass)
{
  VCand res;
  RVecI daughters = GenPart_daughters.at(v_idx);
  if (daughters.size() != 2)
  {
    throw analysis::exception("Vector boson candidate has wrong number of daughters");
  }

  bool decays_to_hadrons = std::all_of(daughters.begin(), daughters.end(), [&](int idx){ return std::abs(GenPart_pdgId[idx]) <= PdG::b(); });
  if (decays_to_hadrons)
  {
    AssignHadronicVCand(res, v_idx, GenPart_pdgId, GenPart_daughters,
                        GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);
  }

  bool decays_to_leptons = std::all_of(daughters.begin(), daughters.end(), [&](int idx){ return PdG::isNeutrino(GenPart_pdgId[idx]) || PdG::isLepton(GenPart_pdgId[idx]); });
  if (decays_to_leptons)
  {
    AssignLeptonicVCand(res, v_idx, gen_leptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags,
                        GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);
  }

  if (!decays_to_hadrons && !decays_to_leptons)
  {
    throw analysis::exception("Encountered neither leptonic nor hadronic vector boson decay");
  }

  return res;
}

HVVCand GetGenHVVCandidate(int evt, std::vector<reco_tau::gen_truth::GenLepton> const& gen_leptons,
                           const RVecI& GenPart_pdgId, const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags,
                           const RVecF& GenPart_pt, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenPart_mass,
                           bool throw_error_if_not_found)
{
  try
  {
    int sz = GenPart_pdgId.size();
    std::set<int> HVV_indices;
    for (int i = 0; i < sz; ++i)
    {
      const GenStatusFlags status(GenPart_statusFlags.at(i));
      bool is_higgs = GenPart_pdgId[i] == PdG::Higgs();
      if (!(is_higgs && status.isLastCopy()))
      {
        continue;
      }
      auto const& daughters = GenPart_daughters.at(i);
      auto IsV = [&](int idx){ return GenPart_pdgId[idx] == PdG::Wplus() || GenPart_pdgId[idx] == PdG::Wminus() || GenPart_pdgId[idx] == PdG::Z(); };
      int n_V_from_H = std::count_if(daughters.begin(), daughters.end(), IsV);
      if (n_V_from_H == 2)
      {
        HVV_indices.insert(i);
      }
    }

    if(HVV_indices.empty())
    {
      throw analysis::exception("H->VV not found.");
    }
    if(HVV_indices.size() != 1)
    {
      throw analysis::exception("Multiple H->VV candidates.");
    }

    int const HVV_index = *HVV_indices.begin();
    HVVCand HVV_cand;

    HVV_cand.cand_p4 = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, HVV_index);

    RVecI V_from_H = GenPart_daughters.at(HVV_index);
    int v1_idx = V_from_H.at(0);
    int v2_idx = V_from_H.at(1);

    v1_idx = GetLastCopy(v1_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_daughters);
    v2_idx = GetLastCopy(v2_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_daughters);

    RVecI daughters_v1 = GenPart_daughters.at(v1_idx);
    bool v1_to_hadrons = std::all_of(daughters_v1.begin(), daughters_v1.end(),
                                      [&](int idx){ return std::abs(GenPart_pdgId[idx]) <= PdG::b(); });
    bool v1_to_leptons = std::all_of(daughters_v1.begin(), daughters_v1.end(),
                                      [&](int idx){ return PdG::isNeutrino(GenPart_pdgId[idx]) || PdG::isLepton(GenPart_pdgId[idx]); });

    RVecI daughters_v2 = GenPart_daughters.at(v2_idx);
    bool v2_to_hadrons = std::all_of(daughters_v2.begin(), daughters_v2.end(),
                                      [&](int idx){ return std::abs(GenPart_pdgId[idx]) <= PdG::b(); });
    bool v2_to_leptons = std::all_of(daughters_v2.begin(), daughters_v2.end(),
                                      [&](int idx){ return PdG::isNeutrino(GenPart_pdgId[idx]) || PdG::isLepton(GenPart_pdgId[idx]); });

    bool single_lepton = (v1_to_hadrons && v2_to_leptons) || (v1_to_leptons && v2_to_hadrons);
    bool double_lepton = v2_to_leptons && v1_to_leptons;

    if (!single_lepton && !double_lepton)
    {
      throw analysis::exception("Encountered decay neither in single lepton nor double lepton channel");
    }

    HVV_cand.legs.at(0) = GetGenVCand(evt, v1_idx, gen_leptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags,
                                      GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);
    HVV_cand.legs.at(1) = GetGenVCand(evt, v2_idx, gen_leptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags,
                                      GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);

    if (single_lepton)
    {
      // let leg 0 be always W->qq and leg 1 leg be W->lv
      if (HVV_cand.legs.at(0).leg_kind.at(0) != Vleg::Jet)
      {
        std::swap(HVV_cand.legs.at(0), HVV_cand.legs.at(1));
      }
    }

    return HVV_cand;
  }
  catch (analysis::exception& e)
  {
    if(throw_error_if_not_found)
    {
      throw analysis::exception("GetGenHVVCandidate (event=%1%): %2%") % evt % e.message();
    }
    return {};
  }
}

HBBCand GetGenHBBCandidate(int evt, const RVecI& GenPart_pdgId, const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags,
                           const RVecF& GenPart_pt, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenPart_mass,
                           bool throw_error_if_not_found)
{
  try
  {
    int sz = GenPart_pdgId.size();
    std::set<int> Hbb_indices;
    for (int i = 0; i < sz; ++i)
    {
      const GenStatusFlags status(GenPart_statusFlags.at(i));
      bool is_higgs = GenPart_pdgId[i] == PdG::Higgs();
      if (!(is_higgs && status.isLastCopy()))
      {
        continue;
      }
      auto const& daughters = GenPart_daughters.at(i);
      auto IsBquark = [&](int idx){ return std::abs(GenPart_pdgId[idx]) == PdG::b(); };
      int n_b_from_H = std::count_if(daughters.begin(), daughters.end(), IsBquark);
      if (n_b_from_H == 2)
      {
        Hbb_indices.insert(i);
      }
    }

    if(Hbb_indices.empty())
    {
      throw analysis::exception("H->bb not found.");
    }
    if(Hbb_indices.size() != 1)
    {
      throw analysis::exception("Multiple H->bb candidates.");
    }

    int const Hbb_index = *Hbb_indices.begin();
    HBBCand Hbb_cand;

    Hbb_cand.cand_p4 = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, Hbb_index);

    RVecI b_from_H = GenPart_daughters.at(Hbb_index);
    int first = GenPart_pt[b_from_H[0]] > GenPart_pt[b_from_H[1]] ? 0 : 1;
    int second = GenPart_pt[b_from_H[0]] > GenPart_pt[b_from_H[1]] ? 1 : 0;
    Hbb_cand.leg_index[0] = b_from_H.at(first);
    Hbb_cand.leg_index[1] = b_from_H.at(second);

    Hbb_cand.leg_p4[0] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, Hbb_cand.leg_index[0]);
    Hbb_cand.leg_p4[1] = GetP4(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, Hbb_cand.leg_index[1]);

    Hbb_cand.leg_vis_p4[0] = LorentzVectorM{};
    Hbb_cand.leg_vis_p4[1] = LorentzVectorM{};

    return Hbb_cand;
  }
  catch (analysis::exception& e)
  {
    if(throw_error_if_not_found)
    {
      throw analysis::exception("GetGenHBBCandidate (event=%1%): %2%") % evt % e.message();
    }
    return {};
  }
}

std::shared_ptr<HTTCand<2>> GetGenHTTCandidate(int evt, const RVecI& GenPart_pdgId,
                                               const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags,
                                               const RVecF& GenPart_pt, const RVecF& GenPart_eta,
                                               const RVecF& GenPart_phi, const RVecF& GenPart_mass,
                                               bool throw_error_if_not_found)
{
  try {
    std::set<int> htt_indices;
    for(int n = 0; n < GenPart_pdgId.size(); ++n) {
        const GenStatusFlags status(GenPart_statusFlags.at(n));
        if(!(GenPart_pdgId[n] == PdG::Higgs() && status.isLastCopy())) continue;
        const auto& daughters = GenPart_daughters.at(n);
        int n_tau_daughters = std::count_if(daughters.begin(), daughters.end(), [&](int idx) {
        return std::abs(GenPart_pdgId.at(idx)) == PdG::tau();
        });
        if(n_tau_daughters == 0) continue;
        if(n_tau_daughters != 2)
        throw analysis::exception("Invalid H->tautau decay. n_tau_daughters = %1%, higgs_idx = %2%")
            % n_tau_daughters % n;
        htt_indices.insert(n);
    }
    if(htt_indices.empty())
        throw analysis::exception("H->tautau not found.");
    if(htt_indices.size() != 1)
        throw analysis::exception("Multiple H->tautau candidates.");
    const int htt_index = *htt_indices.begin();
    HTTCand<2> htt_cand;
    int leg_idx = 0;
    for(int tau_idx : GenPart_daughters.at(htt_index)) {
      if(std::abs(GenPart_pdgId.at(tau_idx)) != PdG::tau()) continue;
      tau_idx = GetLastCopy(tau_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_daughters);

      int lepton_idx = tau_idx;
      size_t n_neutrinos = 0;
      for(int tau_daughter : GenPart_daughters.at(tau_idx)) {
          if(PdG::isNeutrino(GenPart_pdgId[tau_daughter]))
              ++n_neutrinos;
          const GenStatusFlags status(GenPart_statusFlags.at(tau_daughter));
          const int tau_daughter_pdg = std::abs(GenPart_pdgId[tau_daughter]);
          if(!((tau_daughter_pdg == PdG::e() || tau_daughter_pdg == PdG::mu())
              && status.isDirectPromptTauDecayProduct())) continue;
          if(lepton_idx != tau_idx)
              throw analysis::exception("Invalid tau decay. tau_idx = %1%") % tau_idx;
          lepton_idx = tau_daughter;
      }

      if(!((lepton_idx == tau_idx && n_neutrinos == 1) || (lepton_idx != tau_idx && n_neutrinos == 2)))
          throw analysis::exception("Invalid number of neutrinos = %1% in tau decay. tau_idx = %2%")
              % n_neutrinos % tau_idx;

      lepton_idx = GetLastCopy(lepton_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_daughters);
      htt_cand.leg_index.at(leg_idx) = lepton_idx;
      ++leg_idx;
    }
    int leg0_pdg = std::abs(GenPart_pdgId.at(htt_cand.leg_index.at(0)));
    int leg1_pdg = std::abs(GenPart_pdgId.at(htt_cand.leg_index.at(1)));
    if(leg0_pdg > leg1_pdg || (leg0_pdg == leg1_pdg
            && GenPart_pt.at(htt_cand.leg_index.at(0)) < GenPart_pt.at(htt_cand.leg_index.at(1))))
        std::swap(htt_cand.leg_index.at(0), htt_cand.leg_index.at(1));

    for(leg_idx = 0; leg_idx < htt_cand.leg_index.size(); ++leg_idx) {
        const int genPart_index = htt_cand.leg_index.at(leg_idx);
        const int genPart_pdg = GenPart_pdgId.at(genPart_index);
        const auto& genPart_info = ParticleDB::GetParticleInfo(genPart_pdg);
        htt_cand.leg_type[leg_idx] = PdGToLeg(genPart_pdg);
        htt_cand.leg_charge[leg_idx] = genPart_info.charge;
        htt_cand.leg_p4[leg_idx] = GetVisibleP4(genPart_index, GenPart_pdgId, GenPart_daughters,
                                                GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);
    }

    return std::make_shared<HTTCand<2>>(htt_cand);
  } catch(analysis::exception& e) {
    if(throw_error_if_not_found)
      throw analysis::exception("GetGenHTTCandidate (event=%1%): %2%") % evt % e.message();
    return std::make_shared<HTTCand<2>>();
  }
}

int GetGenHBBIndex(int evt, const RVecI& GenPart_pdgId,
                   const RVecVecI& GenPart_daughters, const RVecI& GenPart_statusFlags)
{
  try {
    std::set<int> hbb_indices;
    for(int n = 0; n < GenPart_pdgId.size(); ++n) {
        const GenStatusFlags status(GenPart_statusFlags.at(n));
        if(!(GenPart_pdgId[n] == PdG::Higgs() && status.isLastCopy())) continue;
        const auto& daughters = GenPart_daughters.at(n);
        int n_b_daughters = std::count_if(daughters.begin(), daughters.end(), [&](int idx) {
        return std::abs(GenPart_pdgId.at(idx)) == PdG::b();
        });
        if(n_b_daughters == 0) continue;
        if(n_b_daughters != 2)
        throw analysis::exception("Invalid H->bb decay. n_b_daughters = %1%, higgs_idx = %2%")
            % n_b_daughters % n;
        hbb_indices.insert(n);
    }
    if(hbb_indices.empty())
        throw analysis::exception("H->bb not found.");
    if(hbb_indices.size() != 1)
        throw analysis::exception("Multiple H->bb candidates.");
    const int hbb_index = *hbb_indices.begin();
    return hbb_index;
  }
  catch(analysis::exception& e) {
      throw analysis::exception("GetGenHBBCandidate (event=%1%): %2%") % evt % e.message();
    }
}

bool PassGenAcceptance(const HTTCand<2>& HTT_Cand){
    for(size_t i = 0; i < HTT_Cand.leg_p4.size(); ++i){
        if(!(HTT_Cand.leg_p4.at(i).pt()>20 && std::abs(HTT_Cand.leg_p4.at(i).eta())<2.3 )){
            return false;
        }
    }
    return true;
}

RVecB FindTwoJetsClosestToMPV(float mpv, const RVecLV& GenJet_p4, const RVecB& pre_sel)
{
  int i_min=-1, j_min = -1;
  float delta_min = std::numeric_limits<float>::infinity();
  for(int i = 0; i < GenJet_p4.size(); i++) {
    if(!pre_sel[i]) continue;
    for(int j = 0; j < i; j++) {
      if(!pre_sel[j]) continue;
      const float inv_mass = (GenJet_p4[i]+GenJet_p4[j]).M();
      const float delta_mass = std::abs(inv_mass - mpv);
      if(delta_mass < delta_min) {
        i_min = i;
        j_min = j;
        delta_min = delta_mass;
      }
    }
  }
  RVecB result(pre_sel.size(), false);
  if(i_min >= 0 && j_min>=0) {
    result[i_min] = true;
    result[j_min] = true;
  }
  return result;
}
RVecB FindGenJetAK8(const RVecF& GenJetAK8_mass, const RVecB& pre_sel){

  int i_max = -1;
  float max_mass = -1.;
  for(int i = 0; i < GenJetAK8_mass.size(); i++) {
    if(!pre_sel[i]) continue;
    if(GenJetAK8_mass[i]>max_mass){
      i_max=i;
      max_mass = GenJetAK8_mass[i];
    }
  }
  RVecB result(pre_sel.size(), false);
  if(i_max >= 0 ) {
    result[i_max] = true;
  }
  return result;
}

