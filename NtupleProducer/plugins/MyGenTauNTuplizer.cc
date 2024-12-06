// -*- C++ -*-
//
// Package:    Giovanni/NTuplizer
// Class:      NTuplizer
// 
/**\class NTuplizer NTuplizer.cc Giovanni/NTuplizer/plugins/NTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Petrucciani
//         Created:  Thu, 01 Sep 2016 11:30:38 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"

#include <cstdint>
#include <TTree.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


namespace gentautools{

bool useoffsets = true;
const float& catchInfs(const float& in,const float& replace_value){
        if(in==in){
            if(std::isinf(in))
                return replace_value;
            else if(in < -1e32 || in > 1e32)
                return replace_value;
            return in;
        }
        return replace_value;
}


float catchInfsAndBound(const float& in,const float& replace_value,
        const float& lowerbound, const float& upperbound,const float offset=0){
    float withoutinfs=catchInfs(in,replace_value);
    if(withoutinfs+offset<lowerbound) return lowerbound;
    if(withoutinfs+offset>upperbound) return upperbound;
    if(useoffsets)
        withoutinfs+=offset;
    return withoutinfs;
}


double etaRel(const math::XYZVector& dir, const math::XYZVector& track) {
    double momPar = dir.Dot(track);
    // double energy = std::sqrt(track.Mag2() + ROOT::Math::Square(reco::ParticleMasses::piPlus));
    double energy = std::sqrt(track.Mag2() + ROOT::Math::Square(0.13957));

    return 0.5 * std::log((energy + momPar) / (energy - momPar));
}

  // Sorters to order object collections in decreasing order of pT
  template<typename T> 
  class PatPtSorter {
  public:
    bool operator()(const T& i, const T& j) const {
      return (i.pt() > j.pt());
    }
  };

  PatPtSorter<l1t::PFCandidate>  l1PFCandidateSorter;


  template<typename T> 
  class PatRefPtSorter {
  public:
    bool operator()(const T& i, const T& j) const {
      return (i->pt() > j->pt());
    }
  };
  PatRefPtSorter<l1t::PFJetRef>    jetRefSorter;
  PatRefPtSorter<l1t::PFTauRef>    tauRefSorter;
  PatRefPtSorter<reco::GenJetRef>  genJetRefSorter;

}

using namespace gentautools;

class MyGenTauNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
    public:
        explicit MyGenTauNTuplizer(const edm::ParameterSet&);
        ~MyGenTauNTuplizer();

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {
            // bZ_ = 3.8112; // avoid loading the event setup
        }
        virtual void endRun(edm::Run const&, edm::EventSetup const& iSetup) override { } // framework wants this to be implemented

        void fill_genParticles(const edm::Event& iEvent);

        template<unsigned int bits=10>
        static float zip(float f) {
            return MiniFloatConverter::reduceMantissaToNbitsRounding<bits>(f);
        }

        edm::EDGetTokenT<std::vector<reco::GenJet>> genjets_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        edm::EDGetTokenT<std::vector<l1t::PFJet>> scjets_;
        edm::EDGetTokenT<std::vector<l1t::PFTau>> nntaus_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavour_;
        edm::EDGetTokenT<std::vector<l1t::VertexWord>> const fVtxEmu_;
        TTree *tree_;
        uint32_t run_, lumi_; uint64_t event_;

        float dRJetGenMatch_ = 0.4;
        bool  isMC_;
        float jetPtMin_ = 10;
        float jetEtaMin_ = -2.4;
        float jetEtaMax_ = 2.4;
        float jetPFCandidatePtMin_ = 0.;

        static constexpr size_t max_pfcand_ = 16;
        bool applyJetIDFlag = false;
        /// thresholds for matching
        float dRCone        = 0.2;
        float ptGenLeptonMin = 8;
        float ptGenTauVisibleMin = 15;



    // --------------------
    std::vector<reco::GenParticle> gToBB_;
    std::vector<reco::GenParticle> gToCC_;
    std::vector<reco::GenParticle> neutrinosLepB_;
    std::vector<reco::GenParticle> neutrinosLepB_C_;
    std::vector<reco::GenParticle> alltaus_;
    std::vector<reco::GenParticle> Bhadron_;
    std::vector<reco::GenParticle> Bhadron_daughter_;

    // --------------------
    // from PNET
    // Generator-level information (GEN particles)
    std::vector<float>  gen_particle_pt;
    std::vector<float>  gen_particle_eta;
    std::vector<float>  gen_particle_phi;
    std::vector<float>  gen_particle_mass;
    std::vector<int>    gen_particle_id;
    std::vector<unsigned int>  gen_particle_status;
    std::vector<int>    gen_particle_daughters_id;
    std::vector<unsigned int> gen_particle_daughters_igen;
    std::vector<unsigned int> gen_particle_daughters_status;
    std::vector<float>  gen_particle_daughters_pt;
    std::vector<float>  gen_particle_daughters_eta;
    std::vector<float>  gen_particle_daughters_phi;
    std::vector<float>  gen_particle_daughters_mass;
    std::vector<int>    gen_particle_daughters_charge;

    // Gen leptons from resonance decay 
    std::vector<TLorentzVector> genLepFromResonance4V_;
    std::vector<TLorentzVector> genMuonsFromResonance4V_;
    std::vector<TLorentzVector> genElectronsFromResonance4V_;
    std::vector<TLorentzVector> tau_gen_visible_;
    std::vector<TLorentzVector> tau_gen_;
    std::vector<int> tau_gen_charge_;
    std::vector<unsigned int> tau_gen_nch_;
    std::vector<unsigned int> tau_gen_np0_;
    std::vector<unsigned int> tau_gen_nnh_;

    int gentau_tauflav_, gentau_muflav_, gentau_elflav_, gentau_taudecaymode_, tau_genlepflav_;
    int gentau_taucharge_;
    float gentau_genmatch_lep_pt_;
    float gentau_genmatch_lep_vis_pt_;
    // --------------------
    bool gentau_reject_;

    float gentau_eta_;
    float gentau_phi_;
    float gentau_pt_;
    float gentau_mass_;
    float gentau_energy_;

    float recojet_pt_;
    float recojet_eta_;
    float recojet_phi_;
    float recojet_mass_;
    float recojet_match_dR_;

    float recotau_tauscore_;
    float recotau_pt_;
    float recotau_eta_;
    float recotau_phi_;
    float recotau_mass_;
    float recotau_match_dR_;


    // --------------------

};

MyGenTauNTuplizer::MyGenTauNTuplizer(const edm::ParameterSet& iConfig) :
    genjets_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    scjets_(consumes<std::vector<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("scPuppiJets"))), // l1tSCPFL1PuppiEmulator
    nntaus_(consumes<std::vector<l1t::PFTau>>(iConfig.getParameter<edm::InputTag>("nnTaus"))), 
    genJetsFlavour_   (consumes<reco::JetFlavourInfoMatchingCollection >    (iConfig.getParameter<edm::InputTag>("genJetsFlavour"))),
    fVtxEmu_(consumes<std::vector<l1t::VertexWord>>(iConfig.getParameter<edm::InputTag>("vtx")))
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("GenTaus","GenTaus");

    tree_->Branch("gentau_eta", &gentau_eta_);
    tree_->Branch("gentau_phi", &gentau_phi_);
    tree_->Branch("gentau_pt", &gentau_pt_);
    tree_->Branch("gentau_mass", &gentau_mass_);
    tree_->Branch("gentau_energy", &gentau_energy_);

    tree_->Branch("recojet_pt", &recojet_pt_);
    tree_->Branch("recojet_eta", &recojet_eta_);
    tree_->Branch("recojet_phi", &recojet_phi_);
    tree_->Branch("recojet_mass", &recojet_mass_);
    tree_->Branch("recojet_match_dR", &recojet_match_dR_);

    tree_->Branch("recotau_tauscore", &recotau_tauscore_);
    tree_->Branch("recotau_pt", &recotau_pt_);
    tree_->Branch("recotau_eta", &recotau_eta_);
    tree_->Branch("recotau_phi", &recotau_phi_);
    tree_->Branch("recotau_mass", &recotau_mass_);
    tree_->Branch("recotau_match_dR", &recotau_match_dR_);
    
    tree_->Branch("gentau_taudecaymode", &gentau_taudecaymode_);
    tree_->Branch("gentau_taucharge", &gentau_taucharge_);

    tree_->Branch("gentau_lep_pt", &gentau_genmatch_lep_pt_);
    tree_->Branch("gentau_lep_vis_pt", &gentau_genmatch_lep_vis_pt_);

    // -------------------------------------
    // settings for output TFile and TTree optimized for nanoAOD
    fs->file().SetCompressionAlgorithm(ROOT::ECompressionAlgorithm::kLZ4);
    fs->file().SetCompressionLevel(4);
    for (int idx = 0; idx < tree_->GetListOfBranches()->GetEntries(); ++idx) {
        TBranch* br = dynamic_cast<TBranch*>(tree_->GetListOfBranches()->At(idx));
        if (br) {
            br->SetBasketSize(1024 * 1024);
        }
    }
    if (tree_->GetListOfBranches()->GetEntries() > 0) {
        tree_->SetAutoFlush(-1024 * 1024 * tree_->GetListOfBranches()->GetEntries());
    }
    // -------------------------------------
}

MyGenTauNTuplizer::~MyGenTauNTuplizer() { }

// ------------ method called once each job just before starting event loop  ------------
void 
MyGenTauNTuplizer::beginJob()
{

}


// ------------ method called for each event  ------------
void
MyGenTauNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<std::vector<reco::GenJet>> genjets;
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    edm::Handle<std::vector<l1t::PFJet>> scjets;
    edm::Handle<std::vector<l1t::PFTau>> nntaus;

    edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetsFlavour;
    
    if (iEvent.isRealData()){
        isMC_ = false;
    }else{
        isMC_ = true;
    }
    if (isMC_){
        iEvent.getByToken(genjets_, genjets);
        iEvent.getByToken(genparticles_, genparticles);
        // fill genparticle categories once per event
        fill_genParticles(iEvent);

        iEvent.getByToken(genJetsFlavour_, genJetsFlavour);
        
    }
    iEvent.getByToken(scjets_, scjets);
    iEvent.getByToken(nntaus_, nntaus);


    std::vector<reco::GenJetRef> jetv_gen;  
    if(genjets.isValid()){
        for (auto jets_iter = genjets->begin(); jets_iter != genjets->end(); ++jets_iter) {                                                                                                   
        reco::GenJetRef jref (genjets, jets_iter - genjets->begin());                                                                                                                      
        jetv_gen.push_back(jref);                                                                                                                                                              
        }
        sort(jetv_gen.begin(), jetv_gen.end(), genJetRefSorter);
    }

    std::vector<l1t::PFJetRef> jetv_l1;
    for (auto jets_iter = scjets->begin(); jets_iter != scjets->end(); ++jets_iter) {                                                                                                   
        l1t::PFJetRef jref(scjets, jets_iter - scjets->begin());                                                                                                                             
        // if (jref->pt() < 10.) continue;             
        if (jref->pt() < 1.) continue;             
        jetv_l1.push_back(jref);                                                                                                                                                              
    }
    sort(jetv_l1.begin(), jetv_l1.end(), jetRefSorter);

    std::vector<l1t::PFTauRef> tauv_l1;
    for (auto tau_iter = nntaus->begin(); tau_iter != nntaus->end(); ++tau_iter) {                                                                                                   
        l1t::PFTauRef tref(nntaus, tau_iter - nntaus->begin());
        // if (tau_iter->pt() < 10.) continue;                                                                                                                                  
        if (tau_iter->pt() < 1.) continue;                                                                                                                                  
        // if (tref->pt() < jetPtMin_) continue;
        // if (fabs(tref->eta()) > jetEtaMax_) continue;                 
        // if (fabs(tref->eta()) < jetEtaMin_) continue;                  
        tauv_l1.push_back(tref);                                                                                                                                                              
    }
    sort(tauv_l1.begin(), tauv_l1.end(), tauRefSorter);

    // collection of gen visible taus: tau_gen_visible_

    // size_t itau = 0; itau < tau_gen_visible_.size(); itau++
    for (size_t i = 0; i < tau_gen_visible_.size(); i++) {
        
        if (tau_gen_visible_[i].Pt() < 3.) continue;
        if (fabs(tau_gen_visible_[i].Eta()) > 3.) continue;                 
        // if (fabs(tau_gen_visible_[i].Eta()) < 3.) continue;    

        gentau_pt_ = tau_gen_visible_[i].Pt();
        gentau_eta_ = tau_gen_visible_[i].Eta();
        gentau_phi_ = tau_gen_visible_[i].Phi();
        gentau_mass_ = tau_gen_visible_[i].M();
        gentau_energy_ = tau_gen_visible_[i].Energy();

        // gentau_tauscore_ = tauv_l1[i]->chargedIso();

        // match to GEN
        // int   pos_matched = -1;
        // float minDR = dRJetGenMatch_;
        // for(size_t igen = 0; igen < jetv_gen.size(); igen++){
        //     if(tauv_l1[i]->pt() <= 0.1 * jetv_gen[igen]->pt()) continue;
        //     if(reco::deltaR(jetv_gen[igen]->p4(),tauv_l1[i]->p4()) < minDR){
        //         pos_matched = igen;
        //         minDR = reco::deltaR(jetv_gen[igen]->p4(),tauv_l1[i]->p4());
        //     }
        // }
        
        // if(pos_matched >= 0){
        //     gentau_genmatch_pt_ = jetv_gen[pos_matched]->pt();
        //     gentau_genmatch_eta_ = jetv_gen[pos_matched]->eta();
        //     gentau_genmatch_phi_ = jetv_gen[pos_matched]->phi();
        //     gentau_genmatch_mass_ = jetv_gen[pos_matched]->mass();
        //     gentau_genmatch_dR_ = minDR;
        //     gentau_genmatch_hflav_ = (*genJetsFlavour)[edm::RefToBase<reco::Jet>(jetv_gen[pos_matched])].getHadronFlavour();
        //     gentau_genmatch_pflav_ = (*genJetsFlavour)[edm::RefToBase<reco::Jet>(jetv_gen[pos_matched])].getPartonFlavour();        
        // }
        // else{
        //     gentau_genmatch_pt_ = 0;
        //     gentau_genmatch_eta_ = 0;
        //     gentau_genmatch_phi_ = 0;
        //     gentau_genmatch_mass_ = 0;
        //     gentau_genmatch_dR_ = 0;
        //     gentau_genmatch_hflav_ = 0;
        //     gentau_genmatch_pflav_ = 0;
        // }


        // matching with gen-leptons (muons/electrons/hadronic taus)
        // minDR = 1000;
        // int nlep_in_cone  = 0;
        // int pos_matched_genmu = -1;
        // int pos_matched_genele = -1;
        // int pos_matched_tauh = -1;
        // int gentau_decaymode = -1; 	  
        // TLorentzVector genLepton4V;
        // TLorentzVector genLeptonVis4V;
        // TLorentzVector jet4V;
        // jet4V.SetPtEtaPhiM(tauv_l1[i]->pt(), tauv_l1[i]->eta(), tauv_l1[i]->phi(), tauv_l1[i]->mass());
        // for(size_t igen = 0; igen < genMuonsFromResonance4V_.size(); igen++){
        //     float dR = jet4V.DeltaR(genMuonsFromResonance4V_.at(igen));	      
        //     if(dR < dRCone) nlep_in_cone++;
        //     if(dR < dRCone and dR < minDR){
        //     pos_matched_genmu = igen;
        //     minDR = dR;
        //     genLepton4V = genMuonsFromResonance4V_.at(igen);
        //     genLeptonVis4V = genMuonsFromResonance4V_.at(igen);
        //     }
        // }

        // for(size_t igen = 0; igen < genElectronsFromResonance4V_.size(); igen++){
        //     float dR = jet4V.DeltaR(genElectronsFromResonance4V_.at(igen));	      
        //     if(dR < dRCone) nlep_in_cone++;
        //     if(dR < dRCone and dR < minDR){
        //     pos_matched_genmu  = -1;
        //     pos_matched_genele = igen;
        //     minDR = dR;
        //     genLepton4V = genElectronsFromResonance4V_.at(igen);
        //     genLeptonVis4V = genElectronsFromResonance4V_.at(igen);
        //     }
        // }

        // for(size_t itau = 0; itau < tau_gen_visible_.size(); itau++){
        //     float dR = tau_gen_visible_.at(itau).DeltaR(jet4V); 
        //     if(dR < dRCone) nlep_in_cone++;
        //     if(dR < dRCone and dR < minDR){
        //     pos_matched_genmu  = -1;
        //     pos_matched_genele = -1;
        //     pos_matched_tauh = itau;
        //     minDR = dR;
            // gentau_decaymode = 5*(tau_gen_nch_.at(itau)-1)+tau_gen_np0_.at(itau);
        //     genLepton4V = tau_gen_.at(itau);
        //     genLeptonVis4V = tau_gen_visible_.at(itau);
        //     }
        // }
        

        // gentau_reject_  = false;
        // // exclude, when a jet is matched with a lepton, those for which the matched lepton is below the chosen pt threshold
        // // Jet id applied only to jets not overlapping with gen-leptons
        // if(applyJetIDFlag && pos_matched_genmu  == -1 && pos_matched_genele == -1 && pos_matched_tauh == -1){
        //     gentau_reject_ = true;
        // }
        // if(pos_matched_genmu != -1 and genLeptonVis4V.Pt() < ptGenLeptonMin){
        //     gentau_reject_ = true;
        // }
        // if(pos_matched_genele != -1 and genLeptonVis4V.Pt() < ptGenLeptonMin){
        //     gentau_reject_ = true;
        // }
        // if(pos_matched_tauh != -1 and genLeptonVis4V.Pt() < ptGenTauVisibleMin){
        //     gentau_reject_ = true;
        // } 


        // if(pos_matched_genmu >= 0){
        //     gentau_muflav_  = 1;
        // }
        // else{
        //     gentau_muflav_  = 0;
        // }

        // if(pos_matched_genele >= 0){
        //     gentau_elflav_  = 1;
        // }
        // else{
        //     gentau_elflav_  = 0;
        // }

        // if(pos_matched_tauh >= 0){
            gentau_tauflav_ = 1;
            gentau_taudecaymode_ = 5*(tau_gen_nch_.at(i)-1)+tau_gen_np0_.at(i);;
            gentau_taucharge_ = tau_gen_charge_.at(i);
        // }
        // else{
        //     gentau_tauflav_ = 0;
        //     gentau_taudecaymode_ = -1;
        //     gentau_taucharge_ = 0;
        // }

        // gentau_lepflav_ = nlep_in_cone;

        gentau_genmatch_lep_pt_ = tau_gen_[i].Pt();
        gentau_genmatch_lep_vis_pt_ = tau_gen_visible_[i].Pt();

        // std::cout<<gentau_pt_<<" "<<gentau_genmatch_lep_vis_pt_<<" "<<tau_gen_visible_[i].Pt()<<std::endl;


        TLorentzVector temp4V;

        // match to reco JETS
        int   pos_matched_jet = -1;
        float minDR_jet = 0.4;
        for(size_t ijet = 0; ijet < jetv_l1.size(); ijet++){
            temp4V.SetPtEtaPhiM(jetv_l1[ijet]->pt(), jetv_l1[ijet]->eta(), jetv_l1[ijet]->phi(), jetv_l1[ijet]->mass());
            // if(reco::deltaR(jetv_l1[ijet]->p4(), tau_gen_visible_[i]) < minDR_jet){
            if(tau_gen_visible_[i].DeltaR(temp4V) < minDR_jet){
                pos_matched_jet = ijet;
                minDR_jet = tau_gen_visible_[i].DeltaR(temp4V);
            }
        }
        if (pos_matched_jet > -1){
            recojet_pt_ = jetv_l1[pos_matched_jet]->pt();
            recojet_eta_ = jetv_l1[pos_matched_jet]->eta();
            recojet_phi_ = jetv_l1[pos_matched_jet]->phi();
            recojet_mass_ = jetv_l1[pos_matched_jet]->mass();
            recojet_match_dR_ = minDR_jet;
        }else{
            recojet_pt_ = 0.;
            recojet_eta_ = 0.;
            recojet_phi_ = 0.;
            recojet_mass_ = 0.;
            recojet_match_dR_ = 999.;
        }

        // match to reco TAUS
        int   pos_matched_tau = -1;
        float minDR_tau = 0.4;
        for(size_t itau = 0; itau < tauv_l1.size(); itau++){
            temp4V.SetPtEtaPhiM(tauv_l1[itau]->pt(),tauv_l1[itau]->eta(),tauv_l1[itau]->phi(),tauv_l1[itau]->mass());
            // if(reco::deltaR(tauv_l1[itau]->p4(), tau_gen_visible_[i]) < minDR_tau){
            if(tau_gen_visible_[i].DeltaR(temp4V) < minDR_tau){
                pos_matched_tau = itau;
                minDR_tau = tau_gen_visible_[i].DeltaR(temp4V);
            }
        }
        if (pos_matched_tau > -1){
            recotau_pt_ =  tauv_l1[pos_matched_tau]->pt();
            recotau_eta_ =  tauv_l1[pos_matched_tau]->eta();
            recotau_phi_ =  tauv_l1[pos_matched_tau]->phi();
            recotau_mass_ =  tauv_l1[pos_matched_tau]->mass();
            recotau_match_dR_ = minDR_tau;
            recotau_tauscore_ =tauv_l1[pos_matched_tau]->chargedIso();
        }else{
            recotau_pt_ = 0.;
            recotau_eta_ = 0.;
            recotau_phi_ = 0.;
            recotau_mass_ = 0.;
            recotau_match_dR_ = 999.;
            recotau_tauscore_ = -1.;
        }


        tree_->Fill();
    }


}





void MyGenTauNTuplizer::fill_genParticles(const edm::Event& iEvent)
{
    gToBB_.clear();
    gToCC_.clear();
    neutrinosLepB_.clear();
    neutrinosLepB_C_.clear();
    alltaus_.clear();
    Bhadron_.clear();
    Bhadron_daughter_.clear();

    // Generator-level information (GEN particles)
    gen_particle_pt.clear();
    gen_particle_eta.clear();
    gen_particle_phi.clear();
    gen_particle_mass.clear();
    gen_particle_id.clear();
    gen_particle_status.clear();
    gen_particle_daughters_id.clear();
    gen_particle_daughters_igen.clear();
    gen_particle_daughters_status.clear();
    gen_particle_daughters_pt.clear();
    gen_particle_daughters_eta.clear();
    gen_particle_daughters_phi.clear();
    gen_particle_daughters_mass.clear();
    gen_particle_daughters_charge.clear();

    genLepFromResonance4V_.clear();
    genMuonsFromResonance4V_.clear();
    genElectronsFromResonance4V_.clear();
    tau_gen_visible_.clear();
    tau_gen_.clear();
    tau_gen_charge_.clear();
    tau_gen_nch_.clear();
    tau_gen_np0_.clear();
    tau_gen_nnh_.clear();


    if(!iEvent.isRealData())
    {
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken(genparticles_, genParticles);

        for (const reco::Candidate &genC : *genParticles){
            const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
            if((abs(gen.pdgId())>500&&abs(gen.pdgId())<600)||(abs(gen.pdgId())>5000&&abs(gen.pdgId())<6000)){
                Bhadron_.push_back(gen);
                if(gen.numberOfDaughters()>0){
                    if( (abs(gen.daughter(0)->pdgId())>500&&abs(gen.daughter(0)->pdgId())<600)||(abs(gen.daughter(0)->pdgId())>5000&&abs(gen.daughter(0)->pdgId())<6000)){
                        if(gen.daughter(0)->numberOfDaughters()>0){
                            const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*(gen.daughter(0)->daughter(0)));
                            if(daughter_.vx()!=gen.vx()){
                                Bhadron_daughter_.push_back(daughter_);
                            }else Bhadron_daughter_.push_back(gen);
                        }else  Bhadron_daughter_.push_back(gen);
                    }else{
                        const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*gen.daughter(0));
                        Bhadron_daughter_.push_back(daughter_);
                    }
                }else {
                    Bhadron_daughter_.push_back(gen);
                }
            }
        }

        for (const reco::Candidate &genC : *genParticles){
            const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
            if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16){
                const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (gen.mother());
                if(mother!=NULL){
                    if((abs(mother->pdgId())>500&&abs(mother->pdgId())<600)||(abs(mother->pdgId())>5000&&abs(mother->pdgId())<6000)){
                        neutrinosLepB_.emplace_back(gen);
                    }
                    if((abs(mother->pdgId())>400&&abs(mother->pdgId())<500)||(abs(mother->pdgId())>4000&&abs(mother->pdgId())<5000)){
                        neutrinosLepB_C_.emplace_back(gen);
                    }
                }else {
                    std::cout << "No mother" << std::endl;
                }
            }
            int id(std::abs(gen.pdgId()));
            int status(gen.status());
            if (id == 21 && status >= 21 && status <= 59){ //// Pythia8 hard scatter, ISR, or FSR
                if ( gen.numberOfDaughters() == 2 ){
                    const reco::Candidate* d0 = gen.daughter(0);
                    const reco::Candidate* d1 = gen.daughter(1);
                    if ( std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5 && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4){
                        gToBB_.push_back(gen);
                    }
                    if ( std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4 && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4){
                        gToCC_.push_back(gen);
                    }
                }
            }
            if(id == 15 && false){
                alltaus_.push_back(gen);
            }
        }

        // ---------------------------------------
        // from PNET

  // GEN particle informatio
    if(genParticles.isValid()){
        unsigned int igen = 0;
        for (auto gens_iter = genParticles->begin(); gens_iter != genParticles->end(); ++gens_iter) {      

            if((abs(gens_iter->pdgId()) == 25 or abs(gens_iter->pdgId()) == 24 or abs(gens_iter->pdgId()) == 23) and
            gens_iter->isLastCopy() and 
            gens_iter->statusFlags().fromHardProcess()){ 

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
                    gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
                    gen_particle_daughters_igen.push_back(igen);
                    gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
                    gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
                    gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
                    gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());
                    gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
                    gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
                }
                igen++;
            }

            // Final states Leptons (e,mu) and Neutrinos --> exclude taus. They need to be prompt or from Tau decay      
            if (abs(gens_iter->pdgId()) > 10 and abs(gens_iter->pdgId()) < 17 and abs(gens_iter->pdgId()) != 15  and 
            (gens_iter->isPromptFinalState() or 
            gens_iter->isDirectPromptTauDecayProductFinalState())) { 

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                // No need to save daughters here
                igen++;
            }

            // Final state quarks or gluons from the hard process before the shower --> partons in which H/Z/W/top decay into
            if (((abs(gens_iter->pdgId()) >= 1 and abs(gens_iter->pdgId()) <= 5) or abs(gens_iter->pdgId()) == 21) and 
            gens_iter->statusFlags().fromHardProcess() and 
            gens_iter->statusFlags().isFirstCopy()){
                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());
                igen++;
                // no need to save daughters
            }

            // Special case of taus: last-copy, from hard process and, prompt and decayed
            if(abs(gens_iter->pdgId()) == 15 and 
            gens_iter->isLastCopy() and
            gens_iter->statusFlags().fromHardProcess() and
            gens_iter->isPromptDecayed()){ // hadronic taus

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                // only store the final decay particles
                for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
                    if(not dynamic_cast<const reco::GenParticle*>(gens_iter->daughter(idau))->statusFlags().isPromptTauDecayProduct()) continue;
                    gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
                    gen_particle_daughters_igen.push_back(igen);
                    gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
                    gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
                    gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
                    gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());    
                    gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
                    gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
                }
                igen++;
            }  
        }
    }



        for(size_t igen = 0; igen < gen_particle_pt.size(); igen++){
            // select resonances like Higgs, W, Z, taus
            if(abs(gen_particle_id.at(igen)) == 25 or
            abs(gen_particle_id.at(igen)) == 23 or
            abs(gen_particle_id.at(igen)) == 24 or
            abs(gen_particle_id.at(igen)) == 15){	
            for(size_t idau = 0; idau < gen_particle_daughters_id.size(); idau++){
                // select electrons or muons from the resonance / tau decay
                if(gen_particle_daughters_igen.at(idau) == igen and
                (abs(gen_particle_daughters_id.at(idau)) == 11 or
                abs(gen_particle_daughters_id.at(idau)) == 13)){
                    TLorentzVector gen4V;
                    gen4V.SetPtEtaPhiM(gen_particle_daughters_pt.at(idau),gen_particle_daughters_eta.at(idau),gen_particle_daughters_phi.at(idau),gen_particle_daughters_mass.at(idau));
                    if(std::find(genLepFromResonance4V_.begin(),genLepFromResonance4V_.end(),gen4V) == genLepFromResonance4V_.end())
                    genLepFromResonance4V_.push_back(gen4V);
                    if(abs(gen_particle_daughters_id.at(idau)) == 13 and 
                    std::find(genMuonsFromResonance4V_.begin(),genMuonsFromResonance4V_.end(),gen4V) == genMuonsFromResonance4V_.end()){
                    genMuonsFromResonance4V_.push_back(gen4V);
                    }
                    if(abs(gen_particle_daughters_id.at(idau)) == 11 and 
                    std::find(genElectronsFromResonance4V_.begin(),genElectronsFromResonance4V_.end(),gen4V) == genElectronsFromResonance4V_.end()){
                        genElectronsFromResonance4V_.push_back(gen4V);		
                    }
                }
            }
            }
        }   

        
        // Gen hadronic taus	
        for(size_t igen = 0; igen < gen_particle_pt.size(); igen++){
            if(abs(gen_particle_id.at(igen)) == 15){ // hadronic or leptonic tau
                TLorentzVector tau_gen_tmp;
                unsigned int tau_gen_nch_tmp(0);
                unsigned int tau_gen_np0_tmp(0);
                unsigned int tau_gen_nnh_tmp(0);
                for(size_t idau = 0; idau < gen_particle_daughters_pt.size(); idau++){
                    if(gen_particle_daughters_igen.at(idau) == igen and
                    abs(gen_particle_daughters_id.at(idau)) != 11 and // no mu
                    abs(gen_particle_daughters_id.at(idau)) != 13 and // no el
                    abs(gen_particle_daughters_id.at(idau)) != 12 and // no neutrinos
                    abs(gen_particle_daughters_id.at(idau)) != 14 and
                    abs(gen_particle_daughters_id.at(idau)) != 16){
                    TLorentzVector tmp4V; 
                    tmp4V.SetPtEtaPhiM(gen_particle_daughters_pt.at(idau),gen_particle_daughters_eta.at(idau),gen_particle_daughters_phi.at(idau),gen_particle_daughters_mass.at(idau));
                    tau_gen_tmp += tmp4V;
                    if (gen_particle_daughters_charge.at(idau) != 0 and gen_particle_daughters_status.at(idau) == 1) tau_gen_nch_tmp ++; // charged particles
                    else if(gen_particle_daughters_charge.at(idau) == 0 and gen_particle_daughters_id.at(idau) == 111) tau_gen_np0_tmp++;
                    else if(gen_particle_daughters_charge.at(idau) == 0 and gen_particle_daughters_id.at(idau) != 111) tau_gen_nnh_tmp++;
                    }
                }	
                if(tau_gen_tmp.Pt() > 0){ // good hadronic tau
                    tau_gen_visible_.push_back(tau_gen_tmp);
                    tau_gen_tmp.SetPtEtaPhiM(gen_particle_pt.at(igen),gen_particle_eta.at(igen),gen_particle_phi.at(igen),gen_particle_mass.at(igen));
                    tau_gen_charge_.push_back((gen_particle_id.at(igen) > 0) ? -1 : 1);
                    // std::cout<<((gen_particle_id.at(igen) > 0) ? -1 : 1)<<std::endl;
                    tau_gen_.push_back(tau_gen_tmp);
                    tau_gen_nch_.push_back(tau_gen_nch_tmp);
                    tau_gen_np0_.push_back(tau_gen_np0_tmp);
                    tau_gen_nnh_.push_back(tau_gen_nnh_tmp);
                }
            }
        }       
    }
}










//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MyGenTauNTuplizer);
