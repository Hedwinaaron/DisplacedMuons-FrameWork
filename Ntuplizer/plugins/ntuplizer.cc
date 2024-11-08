#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

float dxy_value(const reco::GenParticle &p, const reco::Vertex &pv){
    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();
    float pv_x = pv.x();
    float pv_y = pv.y();
  
    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;
}


class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntuplizer(const edm::ParameterSet&);
      ~ntuplizer();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      bool isData = true;
      bool isAOD  = false;
      //
      // --- Tokens and Handles
      //

      // trigger bits
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::Handle<edm::TriggerResults> triggerBits;

      // displacedGlobalMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dglToken;
      edm::Handle<edm::View<reco::Track> > dgls;
      // displacedStandAloneMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dsaToken;
      edm::Handle<edm::View<reco::Track> > dsas;
      // displacedMuons (reco::Muon // pat::Muon)
      edm::EDGetTokenT<edm::View<reco::Muon> > dmuToken;
      edm::Handle<edm::View<reco::Muon> > dmuons;
    
      // PrimaryVertices
      edm::EDGetTokenT<edm::View<reco::Vertex> > thePrimaryVertexCollection;
      edm::Handle<edm::View<reco::Vertex> > primaryvertices;
      // GenParticles
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;

      //
      // --- Variables
      //


      // Trigger tags
      std::vector<std::string> HLTPaths_;
      bool triggerPass[200] = {false};
      
      // Tag trigger paths
      std::vector<std::string> Tag_HLTPaths_;
      bool Tag_triggerPass[200] = {false};
      // Probe trigger paths
      std::vector<std::string> Probe_HLTPaths_;
      bool Probe_triggerPass[200] = {false};
      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      // ----------------------------------
      // displacedGlobalMuons
      // ----------------------------------
      Int_t ndgl = 0;
      Float_t dgl_pt[200] = {0.};
      Float_t dgl_eta[200] = {0.};
      Float_t dgl_phi[200] = {0.};
      Float_t dgl_ptError[200] = {0.};
      Float_t dgl_dxy[200] = {0.};
      Float_t dgl_dz[200] = {0.};
      Float_t dgl_normalizedChi2[200] = {0.};
      Float_t dgl_charge[200] = {0.};
      Int_t dgl_nMuonHits[200] = {0};
      Int_t dgl_nValidMuonHits[200] = {0};
      Int_t dgl_nValidMuonDTHits[200] = {0};
      Int_t dgl_nValidMuonCSCHits[200] = {0};
      Int_t dgl_nValidMuonRPCHits[200] = {0};
      Int_t dgl_nValidStripHits[200] = {0};
      Int_t dgl_nhits[200] = {0};
      Int_t dgl_nLostMuonHits[200] = {0};
      Int_t dgl_nLostMuonDTHits[200] = {0};
      Int_t dgl_nLostMuonCSCHits[200] = {0};
      Int_t dgl_nLostMuonRPCHits[200] = {0};

      // ----------------------------------
      // displacedStandAloneMuons
      // ----------------------------------
      Int_t ndsa = 0;
      Float_t dsa_pt[200] = {0.};
      Float_t dsa_eta[200] = {0.};
      Float_t dsa_phi[200] = {0.};
      Float_t dsa_ptError[200] = {0.};
      Float_t dsa_dxy[200] = {0.};
      Float_t dsa_dz[200] = {0.};
      Float_t dsa_normalizedChi2[200] = {0.};
      Float_t dsa_charge[200] = {0.};
      Int_t dsa_nMuonHits[200] = {0};
      Int_t dsa_nValidMuonHits[200] = {0};
      Int_t dsa_nValidMuonDTHits[200] = {0};
      Int_t dsa_nValidMuonCSCHits[200] = {0};
      Int_t dsa_nValidMuonRPCHits[200] = {0};
      Int_t dsa_nValidStripHits[200] = {0};
      Int_t dsa_nhits[200] = {0};
      Int_t dsa_nLostMuonHits[200] = {0};
      Int_t dsa_nLostMuonDTHits[200] = {0};
      Int_t dsa_nLostMuonCSCHits[200] = {0};
      Int_t dsa_nLostMuonRPCHits[200] = {0};
      Int_t dsa_dtStationsWithValidHits[200] = {0};
      Int_t dsa_cscStationsWithValidHits[200] = {0};

      // ----------------------------------
      // displacedMuons
      // ----------------------------------
      Int_t ndmu = 0;
      Int_t dmu_isDSA[200] = {0};
      Int_t dmu_isDGL[200] = {0};
      Int_t dmu_isDTK[200] = {0};
      Int_t dmu_isMatchesValid[200] = {0};
      Int_t dmu_numberOfMatches[200] = {0};
      Int_t dmu_numberOfChambers[200] = {0};
      Int_t dmu_numberOfChambersCSCorDT[200] = {0};
      Int_t dmu_numberOfMatchedStations[200] = {0};
      Int_t dmu_numberOfMatchedRPCLayers[200] = {0};
      Float_t dmu_dsa_pt[200] = {0.};
      Float_t dmu_dsa_eta[200] = {0.};
      Float_t dmu_dsa_phi[200] = {0.};
      Float_t dmu_dsa_ptError[200] = {0.};
      Float_t dmu_dsa_dxy[200] = {0.};
      Float_t dmu_dsa_dz[200] = {0.};
      Float_t dmu_dsa_normalizedChi2[200] = {0.};
      Float_t dmu_dsa_charge[200] = {0.};
      Int_t dmu_dsa_nMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonDTHits[200] = {0};
      Int_t dmu_dsa_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dsa_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dsa_nValidStripHits[200] = {0};
      Int_t dmu_dsa_nhits[200] = {0};
      Int_t dmu_dsa_dtStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_cscStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_nsegments[200] = {0};
      Float_t dmu_dgl_pt[200] = {0.};
      Float_t dmu_dgl_eta[200] = {0.};
      Float_t dmu_dgl_phi[200] = {0.};
      Float_t dmu_dgl_ptError[200] = {0.};
      Float_t dmu_dgl_dxy[200] = {0.};
      Float_t dmu_dgl_dz[200] = {0.};
      Float_t dmu_dgl_normalizedChi2[200] = {0.};
      Float_t dmu_dgl_charge[200] = {0.};
      Int_t dmu_dgl_nMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonDTHits[200] = {0};
      Int_t dmu_dgl_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dgl_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dgl_nValidStripHits[200] = {0};
      Int_t dmu_dgl_nhits[200] = {0};
      Float_t dmu_dtk_pt[200] = {0.};
      Float_t dmu_dtk_eta[200] = {0.};
      Float_t dmu_dtk_phi[200] = {0.};
      Float_t dmu_dtk_ptError[200] = {0.};
      Float_t dmu_dtk_dxy[200] = {0.};
      Float_t dmu_dtk_dz[200] = {0.};
      Float_t dmu_dtk_normalizedChi2[200] = {0.};
      Float_t dmu_dtk_charge[200] = {0.};
      Int_t dmu_dtk_nMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonDTHits[200] = {0};
      Int_t dmu_dtk_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dtk_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dtk_nValidStripHits[200] = {0};
      Int_t dmu_dtk_nhits[200] = {0};

      // ----------------------------------
      // PrimaryVertices
      // ----------------------------------
      Int_t nPV;
      Int_t nTruePV;
      Int_t PV_passAcceptance;
      Float_t PV_vx;
      Float_t PV_vy;
      Float_t PV_vz;

      // ----------------------------------
      // GenParticles
      // ----------------------------------
      Int_t nGenMuon;
      Int_t nGenMuon_PFS;
      Int_t nGenMuon_HPFS;
      Int_t nGenMuon_PTDP;
      Int_t nGenMuon_HDP;
      Float_t GenLeptonSel_pt[30];
      Float_t GenLeptonSel_E[30];
      Float_t GenLeptonSel_et[30];
      Float_t GenLeptonSel_eta[30];
      Float_t GenLeptonSel_phi[30];
      Int_t GenLeptonSel_pdgId[30];
      Float_t GenLeptonSel_dxy[30];
      Float_t GenLeptonSel_vx[30];
      Float_t GenLeptonSel_vy[30];
      Float_t GenLeptonSel_vz[30];
      Int_t GenLeptonSel_motherPdgId[30];
      Int_t GenLeptonSel_fromHardProcessFinalState[30];
      Int_t GenLeptonSel_isPromptFinalState[30];
      Int_t GenLeptonSel_isDirectPromptTauDecayProductFinalState[30];
      Int_t GenLeptonSel_isDirectHadronDecayProduct[30];
      
      Int_t nHardProcessParticle;
      Float_t HardProcessParticle_pt[30];
      Float_t HardProcessParticle_E[30];
      Float_t HardProcessParticle_eta[30];
      Float_t HardProcessParticle_phi[30];
      Float_t HardProcessParticle_vx[30];
      Float_t HardProcessParticle_vy[30];
      Float_t HardProcessParticle_vz[30];
      Int_t HardProcessParticle_pdgId[30];
      //---------------------------------
      // Variables for hemispheres
      //---------------------------------
      Int_t ndgl_up;
      Int_t ndgl_down;
      Int_t ndmu_upper[200]= {0};
      Int_t ndmu_lower[200]= {0};

      Float_t upper_muon_eta[200]={0.};
      Float_t upper_muon_phi[200]={0.};
      Float_t upper_muon_pts[200]={0.};
      Float_t lower_muon_eta[200]={0.};
      Float_t lower_muon_phi[200]={0.};
      Float_t lower_muon_pts[200]={0.};


      //---------------------------------
      // Variables for hemispheres Tag
      //---------------------------------
      Int_t tag_dgl=0;
      Float_t tag_dgl_pt[200]= {0.};
      Float_t tag_dgl_eta[200]= {0.};
      Float_t tag_dgl_phi[200]= {0.};
      Int_t   is_tag_dgl[200]={0}; 
      Float_t tag_dgl_ptError[200] = {0.};
      Float_t tag_dgl_dxy[200] = {0.};
      Float_t tag_dgl_dz[200] = {0.};
      Float_t tag_dgl_normalizedChi2[200] = {0.};
      Float_t tag_dgl_charge[200] = {0.};
      Int_t tag_dgl_nMuonHits[200] = {0};
      Int_t tag_dgl_nValidMuonHits[200] = {0};
      Int_t tag_dgl_nValidMuonDTHits[200] = {0};
      Int_t tag_dgl_nValidMuonCSCHits[200] = {0};
      Int_t tag_dgl_nValidMuonRPCHits[200] = {0};
      Int_t tag_dgl_nValidStripHits[200] = {0};
      Int_t tag_dgl_nhits[200] = {0};
      
      
      //---------------------------------
      // Variables for hemispheres probe
      //---------------------------------
      Int_t probe_dgl=0;
      Float_t probe_dgl_pt[200]= {0.};
      Float_t probe_dgl_eta[200]= {0.};
      Float_t probe_dgl_phi[200]= {0.};
      Int_t    is_probe_dgl[200]={0};
      Float_t probe_dgl_ptError[200] = {0.};
      Float_t probe_dgl_dxy[200] = {0.};
      Float_t probe_dgl_dz[200] = {0.};
      Float_t probe_dgl_normalizedChi2[200] = {0.};
      Float_t probe_dgl_charge[200] = {0.};
      Int_t probe_dgl_nMuonHits[200] = {0};
      Int_t probe_dgl_nValidMuonHits[200] = {0};
      Int_t probe_dgl_nValidMuonDTHits[200] = {0};
      Int_t probe_dgl_nValidMuonCSCHits[200] = {0};
      Int_t probe_dgl_nValidMuonRPCHits[200] = {0};
      Int_t probe_dgl_nValidStripHits[200] = {0};
      Int_t probe_dgl_nhits[200] = {0}; 
            
      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) {

   usesResource("TFileService");

   parameters = iConfig;

   // Analyzer parameters
   isData = parameters.getParameter<bool>("isData");
   isAOD = parameters.getParameter<bool>("isAOD");  
 
   counts = new TH1F("counts", "", 1, 0, 1);

   dglToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedGlobalCollection"));
   dsaToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedStandAloneCollection"));
   dmuToken = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("displacedMuonCollection"));
   if (!isData) {
     theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));
     thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   }

   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
}


// Destructor
ntuplizer::~ntuplizer() {
}


// beginJob (Before first event)
void ntuplizer::beginJob() {

   std::cout << "Begin Job" << std::endl;

   // Init the file and the TTree
   output_filename = parameters.getParameter<std::string>("nameOfOutput");
   file_out = new TFile(output_filename.c_str(), "RECREATE");
   tree_out = new TTree("Events", "Events");

   // Load HLT paths
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX");
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX");
   // Load Tag filter 
   Tag_HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX_v");
   Tag_HLTPaths_.push_back("hltL1fL1sMuOpenNotBptxORNoHaloMu3BXL1Filtered0");
   Tag_HLTPaths_.push_back("hltL2fL1sMuOpenNotBptxORNoHaloMu3BXL1f0NoVtxCosmicSeedMeanTimerL2Filtered10");
   // Load Tag filter
   Probe_HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX_v");

   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");

   // ----------------------------------
   // displacedGlobalMuons
   // ----------------------------------
   tree_out->Branch("ndgl", &ndgl, "ndgl/I");
   tree_out->Branch("dgl_pt", dgl_pt, "dgl_pt[ndgl]/F");
   tree_out->Branch("dgl_eta", dgl_eta, "dgl_eta[ndgl]/F");
   tree_out->Branch("dgl_phi", dgl_phi, "dgl_phi[ndgl]/F");
   tree_out->Branch("dgl_ptError", dgl_ptError, "dgl_ptError[ndgl]/F");
   tree_out->Branch("dgl_dxy", dgl_dxy, "dgl_dxy[ndgl]/F");
   tree_out->Branch("dgl_dz", dgl_dz, "dgl_dz[ndgl]/F");
   tree_out->Branch("dgl_normalizedChi2", dgl_normalizedChi2, "dgl_normalizedChi2[ndgl]/F");
   tree_out->Branch("dgl_charge", dgl_charge, "dgl_charge[ndgl]/F");
   tree_out->Branch("dgl_nMuonHits", dgl_nMuonHits, "dgl_nMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonHits", dgl_nValidMuonHits, "dgl_nValidMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonDTHits", dgl_nValidMuonDTHits, "dgl_nValidMuonDTHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonCSCHits", dgl_nValidMuonCSCHits, "dgl_nValidMuonCSCHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonRPCHits", dgl_nValidMuonRPCHits, "dgl_nValidMuonRPCHits[ndgl]/I");
   tree_out->Branch("dgl_nValidStripHits", dgl_nValidStripHits, "dgl_nValidStripHits[ndgl]/I");
   tree_out->Branch("dgl_nhits", dgl_nhits, "dgl_nhits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonHits", dgl_nLostMuonHits, "dgl_nLostMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonDTHits", dgl_nLostMuonDTHits, "dgl_nLostMuonDTHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonCSCHits", dgl_nLostMuonCSCHits, "dgl_nLostMuonCSCHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonRPCHits", dgl_nLostMuonRPCHits, "dgl_nLostMuonRPCHits[ndgl]/I");

   // ----------------------------------
   // displacedStandAloneMuons
   // ----------------------------------
   tree_out->Branch("ndsa", &ndsa, "ndsa/I");
   tree_out->Branch("dsa_pt", dsa_pt, "dsa_pt[ndsa]/F");
   tree_out->Branch("dsa_eta", dsa_eta, "dsa_eta[ndsa]/F");
   tree_out->Branch("dsa_phi", dsa_phi, "dsa_phi[ndsa]/F");
   tree_out->Branch("dsa_ptError", dsa_ptError, "dsa_ptError[ndsa]/F");
   tree_out->Branch("dsa_dxy", dsa_dxy, "dsa_dxy[ndsa]/F");
   tree_out->Branch("dsa_dz", dsa_dz, "dsa_dz[ndsa]/F");
   tree_out->Branch("dsa_normalizedChi2", dsa_normalizedChi2, "dsa_normalizedChi2[ndsa]/F");
   tree_out->Branch("dsa_charge", dsa_charge, "dsa_charge[ndsa]/F");
   tree_out->Branch("dsa_nMuonHits", dsa_nMuonHits, "dsa_nMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonHits", dsa_nValidMuonHits, "dsa_nValidMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonDTHits", dsa_nValidMuonDTHits, "dsa_nValidMuonDTHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonCSCHits", dsa_nValidMuonCSCHits, "dsa_nValidMuonCSCHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonRPCHits", dsa_nValidMuonRPCHits, "dsa_nValidMuonRPCHits[ndsa]/I");
   tree_out->Branch("dsa_nValidStripHits", dsa_nValidStripHits, "dsa_nValidStripHits[ndsa]/I");
   tree_out->Branch("dsa_nhits", dsa_nhits, "dsa_nhits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonHits", dsa_nLostMuonHits, "dsa_nLostMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonDTHits", dsa_nLostMuonDTHits, "dsa_nLostMuonDTHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonCSCHits", dsa_nLostMuonCSCHits, "dsa_nLostMuonCSCHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonRPCHits", dsa_nLostMuonRPCHits, "dsa_nLostMuonRPCHits[ndsa]/I");
   tree_out->Branch("dsa_dtStationsWithValidHits", dsa_dtStationsWithValidHits, "dsa_dtStationsWithValidHits[ndsa]/I");
   tree_out->Branch("dsa_cscStationsWithValidHits", dsa_cscStationsWithValidHits, "dsa_cscStationsWithValidHits[ndsa]/I");

   // ----------------------------------
   // displacedMuons
   // ----------------------------------
   tree_out->Branch("ndmu", &ndmu, "ndmu/I");
   tree_out->Branch("dmu_isDSA", dmu_isDSA, "dmu_isDSA[ndmu]/I");
   tree_out->Branch("dmu_isDGL", dmu_isDGL, "dmu_isDGL[ndmu]/I");
   tree_out->Branch("dmu_isDTK", dmu_isDTK, "dmu_isDTK[ndmu]/I");
   tree_out->Branch("dmu_isMatchesValid", dmu_isMatchesValid, "dmu_isMatchesValid[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatches", dmu_numberOfMatches, "dmu_numberOfMatches[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambers", dmu_numberOfChambers, "dmu_numberOfChambers[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambersCSCorDT", dmu_numberOfChambersCSCorDT, "dmu_numberOfChambersCSCorDT[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedStations", dmu_numberOfMatchedStations, "dmu_numberOfMatchedStations[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedRPCLayers", dmu_numberOfMatchedRPCLayers, "dmu_numberOfMatchedRPCLayers[ndmu]/I");
   // dmu_dsa
   tree_out->Branch("dmu_dsa_pt", dmu_dsa_pt, "dmu_dsa_pt[ndmu]/F");
   tree_out->Branch("dmu_dsa_eta", dmu_dsa_eta, "dmu_dsa_eta[ndmu]/F");
   tree_out->Branch("dmu_dsa_phi", dmu_dsa_phi, "dmu_dsa_phi[ndmu]/F");
   tree_out->Branch("dmu_dsa_ptError", dmu_dsa_ptError, "dmu_dsa_ptError[ndmu]/F");
   tree_out->Branch("dmu_dsa_dxy", dmu_dsa_dxy, "dmu_dsa_dxy[ndmu]/F");
   tree_out->Branch("dmu_dsa_dz", dmu_dsa_dz, "dmu_dsa_dz[ndmu]/F");
   tree_out->Branch("dmu_dsa_normalizedChi2", dmu_dsa_normalizedChi2, "dmu_dsa_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dsa_charge", dmu_dsa_charge, "dmu_dsa_charge[ndmu]/F");
   tree_out->Branch("dmu_dsa_nMuonHits", dmu_dsa_nMuonHits, "dmu_dsa_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonHits", dmu_dsa_nValidMuonHits, "dmu_dsa_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonDTHits", dmu_dsa_nValidMuonDTHits, "dmu_dsa_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonCSCHits", dmu_dsa_nValidMuonCSCHits, "dmu_dsa_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonRPCHits", dmu_dsa_nValidMuonRPCHits, "dmu_dsa_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidStripHits", dmu_dsa_nValidStripHits, "dmu_dsa_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nhits", dmu_dsa_nhits, "dmu_dsa_nhits[ndmu]/I");
   tree_out->Branch("dmu_dsa_dtStationsWithValidHits", dmu_dsa_dtStationsWithValidHits, "dmu_dsa_dtStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_cscStationsWithValidHits", dmu_dsa_cscStationsWithValidHits, "dmu_dsa_cscStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nsegments", dmu_dsa_nsegments, "dmu_dsa_nsegments[ndmu]/I");
   // dmu_dgl
   tree_out->Branch("dmu_dgl_pt", dmu_dgl_pt, "dmu_dgl_pt[ndmu]/F");
   tree_out->Branch("dmu_dgl_eta", dmu_dgl_eta, "dmu_dgl_eta[ndmu]/F");
   tree_out->Branch("dmu_dgl_phi", dmu_dgl_phi, "dmu_dgl_phi[ndmu]/F");
   tree_out->Branch("dmu_dgl_ptError", dmu_dgl_ptError, "dmu_dgl_ptError[ndmu]/F");
   tree_out->Branch("dmu_dgl_dxy", dmu_dgl_dxy, "dmu_dgl_dxy[ndmu]/F");
   tree_out->Branch("dmu_dgl_dz", dmu_dgl_dz, "dmu_dgl_dz[ndmu]/F");
   tree_out->Branch("dmu_dgl_normalizedChi2", dmu_dgl_normalizedChi2, "dmu_dgl_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dgl_charge", dmu_dgl_charge, "dmu_dgl_charge[ndmu]/F");
   tree_out->Branch("dmu_dgl_nMuonHits", dmu_dgl_nMuonHits, "dmu_dgl_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonHits", dmu_dgl_nValidMuonHits, "dmu_dgl_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonDTHits", dmu_dgl_nValidMuonDTHits, "dmu_dgl_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonCSCHits", dmu_dgl_nValidMuonCSCHits, "dmu_dgl_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonRPCHits", dmu_dgl_nValidMuonRPCHits, "dmu_dgl_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidStripHits", dmu_dgl_nValidStripHits, "dmu_dgl_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nhits", dmu_dgl_nhits, "dmu_dgl_nhits[ndmu]/I");
   // dmu_dtk
   tree_out->Branch("dmu_dtk_pt", dmu_dtk_pt, "dmu_dtk_pt[ndmu]/F");
   tree_out->Branch("dmu_dtk_eta", dmu_dtk_eta, "dmu_dtk_eta[ndmu]/F");
   tree_out->Branch("dmu_dtk_phi", dmu_dtk_phi, "dmu_dtk_phi[ndmu]/F");
   tree_out->Branch("dmu_dtk_ptError", dmu_dtk_ptError, "dmu_dtk_ptError[ndmu]/F");
   tree_out->Branch("dmu_dtk_dxy", dmu_dtk_dxy, "dmu_dtk_dxy[ndmu]/F");
   tree_out->Branch("dmu_dtk_dz", dmu_dtk_dz, "dmu_dtk_dz[ndmu]/F");
   tree_out->Branch("dmu_dtk_normalizedChi2", dmu_dtk_normalizedChi2, "dmu_dtk_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dtk_charge", dmu_dtk_charge, "dmu_dtk_charge[ndmu]/F");
   tree_out->Branch("dmu_dtk_nMuonHits", dmu_dtk_nMuonHits, "dmu_dtk_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonHits", dmu_dtk_nValidMuonHits, "dmu_dtk_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonDTHits", dmu_dtk_nValidMuonDTHits, "dmu_dtk_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonCSCHits", dmu_dtk_nValidMuonCSCHits, "dmu_dtk_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonRPCHits", dmu_dtk_nValidMuonRPCHits, "dmu_dtk_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidStripHits", dmu_dtk_nValidStripHits, "dmu_dtk_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nhits", dmu_dtk_nhits, "dmu_dtk_nhits[ndmu]/I");

    // Store values in the tree or output (e.g., add branches for upper/lower muons)
    
    
    tree_out->Branch("ndmu_upper", &ndmu_upper, "ndmu_upper/I");
    tree_out->Branch("upper_muon_pts", &upper_muon_pts, "dgl_pt[ndgl]/F");
    tree_out->Branch("upper_muon_eta", &upper_muon_eta, "dgl_eta[ndgl]/F");
    tree_out->Branch("upper_muon", &upper_muon_phi, "dgl_phi[ndgl]/F");
    
    tree_out->Branch("ndmu_lower", &ndmu_lower, "ndmu_lower/I");
    tree_out->Branch("lower_muon_pts", &lower_muon_pts, "dgl_pt[ndgl]/F");
    tree_out->Branch("lower_muon_eta", &lower_muon_eta, "dgl_eta[ndgl]/F");
    tree_out->Branch("lower_muon_phi", &lower_muon_phi, "dgl_phi[ndgl]/F");
  
    //store tag and probe hist
    // Tag
    tree_out->Branch("tag_dgl", &tag_dgl, "tag_dgl/I");
    tree_out->Branch("tag_dgl_pt", tag_dgl_pt , "tag_dgl_pt[tag_dgl]/F");
    tree_out->Branch("tag_dgl_eta", tag_dgl_eta , "tag_dgl_eta[tag_dgl]/F");
    tree_out->Branch("tag_dgl_phi", tag_dgl_phi , "tag_dgl_phi[tag_dgl]/F");
    tree_out->Branch("is_tag_dgl", is_tag_dgl , "is_tag_dgl[tag_dgl]/I");
    tree_out->Branch("tag_dgl_ptError", tag_dgl_ptError, "tag_dgl_ptError[tag_dgl]/F");
    tree_out->Branch("tag_dgl_dxy", tag_dgl_dxy, "tag_dgl_dxy[tag_dgl]/F");
    tree_out->Branch("tag_dgl_dz", tag_dgl_dz, "tag_dgl_dz[tag_dgl]/F");
    tree_out->Branch("tag_dgl_normalizedChi2", tag_dgl_normalizedChi2, "tag_dgl_normalizedChi2[tag_dgl]/F");
    tree_out->Branch("tag_dgl_charge", tag_dgl_charge, "tag_dgl_charge[tag_dgl]/F");
    tree_out->Branch("tag_dgl_nMuonHits", tag_dgl_nMuonHits, "tag_dgl_nMuonHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nValidMuonHits", tag_dgl_nValidMuonHits, "tag_dgl_nValidMuonHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nValidMuonDTHits", tag_dgl_nValidMuonDTHits, "tag_dgl_nValidMuonDTHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nValidMuonCSCHits", tag_dgl_nValidMuonCSCHits, "tag_dgl_nValidMuonCSCHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nValidMuonRPCHits", tag_dgl_nValidMuonRPCHits, "tag_dgl_nValidMuonRPCHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nValidStripHits", tag_dgl_nValidStripHits, "tag_dgl_nValidStripHits[tag_dgl]/I");
    tree_out->Branch("tag_dgl_nhits", tag_dgl_nhits, "tag_dgl_nhits[tag_dgl]/I");
    
    //probe
    tree_out->Branch("probe_dgl", &probe_dgl, "probe_dgl/I");
    tree_out->Branch("probe_dgl_pt", probe_dgl_pt , "probe_dgl_pt[probe_dgl]/F");
    tree_out->Branch("probe_dgl_eta", probe_dgl_eta , "probe_dgl_eta[probe_dgl]/F");
    tree_out->Branch("probe_dgl_phi", probe_dgl_phi , "probe_dgl_phi[probe_dgl]/F");
    tree_out->Branch("is_probe_dgl", is_probe_dgl , "is_probe_dgl[probe_dgl]/I");
    tree_out->Branch("probe_dgl_ptError", probe_dgl_ptError, "probe_dgl_ptError[probe_dgl]/F");
    tree_out->Branch("probe_dgl_dxy", probe_dgl_dxy, "probe_dgl_dxy[probe_dgl]/F");
    tree_out->Branch("probe_dgl_dz", probe_dgl_dz, "probe_dgl_dz[probe_dgl]/F");
    tree_out->Branch("probe_dgl_normalizedChi2", probe_dgl_normalizedChi2, "probe_dgl_normalizedChi2[probe_dgl]/F");
    tree_out->Branch("probe_dgl_charge", probe_dgl_charge, "probe_dgl_charge[probe_dgl]/F");
    tree_out->Branch("probe_dgl_nMuonHits", probe_dgl_nMuonHits, "probe_dgl_nMuonHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nValidMuonHits", probe_dgl_nValidMuonHits, "probe_dgl_nValidMuonHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nValidMuonDTHits", probe_dgl_nValidMuonDTHits, "probe_dgl_nValidMuonDTHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nValidMuonCSCHits", probe_dgl_nValidMuonCSCHits, "probe_dgl_nValidMuonCSCHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nValidMuonRPCHits", probe_dgl_nValidMuonRPCHits, "probe_dgl_nValidMuonRPCHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nValidStripHits", probe_dgl_nValidStripHits, "probe_dgl_nValidStripHits[probe_dgl]/I");
    tree_out->Branch("probe_dgl_nhits", probe_dgl_nhits, "probe_dgl_nhits[probe_dgl]/I");
    
    
   if (!isData) {
     // ----------------------------------
     // PrimaryVertices
     // ----------------------------------
     tree_out->Branch("nPV", &nPV, "nPV/I");
     tree_out->Branch("nTruePV", &nTruePV, "nTruePV/I");
     tree_out->Branch("PV_passAcceptance", &PV_passAcceptance, "PV_passAcceptance/I");
     tree_out->Branch("PV_vx", &PV_vx, "PV_vx/F");
     tree_out->Branch("PV_vy", &PV_vy, "PV_vy/F");
     tree_out->Branch("PV_vz", &PV_vz, "PV_vz/F");   

     // ----------------------------------
     // GenParticles
     // ----------------------------------
     tree_out->Branch("nGenMuon", &nGenMuon, "nGenMuon/I");
     tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_E", GenLeptonSel_E, "GenLeptonSel_E[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_phi", GenLeptonSel_phi, "GenLeptonSel_phi[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_dxy", GenLeptonSel_dxy, "GenLeptonSel_dxy[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vx", GenLeptonSel_vx, "GenLeptonSel_vx[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vy", GenLeptonSel_vy, "GenLeptonSel_vy[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vz", GenLeptonSel_vz, "GenLeptonSel_vz[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_pdgId", GenLeptonSel_pdgId, "GenLeptonSel_pdgId[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_motherPdgId", GenLeptonSel_motherPdgId, "GenLeptonSel_motherPdgId[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isPromptFinalState", GenLeptonSel_isPromptFinalState, "GenLeptonSel_isPromptFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_fromHardProcessFinalState", GenLeptonSel_fromHardProcessFinalState, "GenLeptonSel_fromHardProcessFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isDirectPromptTauDecayProductFinalState", GenLeptonSel_isDirectPromptTauDecayProductFinalState, "GenLeptonSel_isDirectPromptTauDecayProductFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isDirectHadronDecayProduct", GenLeptonSel_isDirectHadronDecayProduct, "GenLeptonSel_isDirectHadronDecayProduct[nGenMuon]/I");
 
     tree_out->Branch("nHardProcessParticle", &nHardProcessParticle, "nHardProcessParticle/I");
     tree_out->Branch("HardProcessParticle_E", HardProcessParticle_E, "HardProcessParticle_E[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_pt", HardProcessParticle_pt, "HardProcessParticle_pt[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_eta", HardProcessParticle_eta, "HardProcessParticle_eta[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_phi", HardProcessParticle_phi, "HardProcessParticle_phi[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vx", HardProcessParticle_vx, "HardProcessParticle_vx[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vy", HardProcessParticle_vy, "HardProcessParticle_vy[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vz", HardProcessParticle_vz, "HardProcessParticle_vz[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_pdgId", HardProcessParticle_pdgId, "HardProcessParticle_pdgId[nHardProcessParticle]/I");
   }

   // Trigger branches
   for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
     tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
   }
   
   // Tag Trigger branches
   for (unsigned int ihlt = 0; ihlt < Tag_HLTPaths_.size(); ihlt++) {
     tree_out->Branch(TString(Tag_HLTPaths_[ihlt]), &Tag_triggerPass[ihlt]);
   }
   // Tag Trigger branches
   for (unsigned int ihlt = 0; ihlt < Probe_HLTPaths_.size(); ihlt++) {
     tree_out->Branch(TString(Probe_HLTPaths_[ihlt]), &Probe_triggerPass[ihlt]);
   }

}

// endJob (After event loop has finished)
void ntuplizer::endJob()
{

    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();

}


// fillDescriptions
void ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   iEvent.getByToken(dglToken, dgls);
   iEvent.getByToken(dsaToken, dsas);
   iEvent.getByToken(dmuToken, dmuons);
   iEvent.getByToken(triggerBits_, triggerBits);
   if (!isData){
      iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
      iEvent.getByToken(theGenParticleCollection, genParticles);
   }

   // Count number of events read
   counts->Fill(0);
   

   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   // ----------------------------------
   // displacedGlobalMuons Collection
   // ----------------------------------
   ndgl = 0;
   for (unsigned int i = 0; i < dgls->size(); i++) {
     const reco::Track& dgl(dgls->at(i));
     dgl_pt[ndgl] = dgl.pt();
     dgl_eta[ndgl] = dgl.eta();
     dgl_phi[ndgl] = dgl.phi();
     dgl_ptError[ndgl] = dgl.ptError();
     dgl_dxy[ndgl] = dgl.dxy();
     dgl_dz[ndgl] = dgl.dz();
     dgl_normalizedChi2[ndgl] = dgl.normalizedChi2();
     dgl_charge[ndgl] = dgl.charge();

     dgl_nMuonHits[ndgl] = dgl.hitPattern().numberOfMuonHits();
     dgl_nValidMuonDTHits[ndgl] = dgl.hitPattern().numberOfValidMuonDTHits();
     dgl_nValidMuonCSCHits[ndgl] = dgl.hitPattern().numberOfValidMuonCSCHits();
     dgl_nValidMuonRPCHits[ndgl] = dgl.hitPattern().numberOfValidMuonRPCHits();
     dgl_nValidMuonHits[ndgl] = dgl.hitPattern().numberOfValidMuonHits();
     dgl_nValidStripHits[ndgl] = dgl.hitPattern().numberOfValidStripHits();
     dgl_nhits[ndgl] = dgl.hitPattern().numberOfValidHits();

     dgl_nLostMuonHits[ndgl] = dgl.hitPattern().numberOfLostMuonHits();
     dgl_nLostMuonDTHits[ndgl] = dgl.hitPattern().numberOfLostMuonDTHits();
     dgl_nLostMuonCSCHits[ndgl] = dgl.hitPattern().numberOfLostMuonCSCHits();
     dgl_nLostMuonRPCHits[ndgl] = dgl.hitPattern().numberOfLostMuonRPCHits();
     ndgl++;
   }
  
  
  // Loop over displacedMuonCollection
  ndgl_up=0;
  ndgl_down=0;
  for (unsigned int i = 0; i < dmuons->size(); i++) {
      const reco::Muon& muon = dmuons->at(i);

      // Get the muon eta value to sort into hemispheres
      float muon_eta = muon.eta();
      float muon_charge = muon.charge();

      if (muon_eta > 0 && muon_charge < 0 ) {
          // Muon in the upper hemisphere
          ndmu_upper[ndgl_up]=muon.numberOfMatches(reco::Muon::SegmentArbitration);
          upper_muon_pts[ndgl_up]=muon.pt(); 
          upper_muon_eta[ndgl_up]=muon.eta();
          upper_muon_phi[ndgl_up]=muon.phi();
          ndgl_up++;
      } else if (muon_eta < 0 && muon_charge < 0) {
          // Muon in the lower hemisphere
          ndmu_lower[ndgl_down]=muon.numberOfMatches(reco::Muon::SegmentArbitration);
          lower_muon_pts[ndgl_down]=muon.pt();
          lower_muon_eta[ndgl_up]=muon.eta();
          lower_muon_phi[ndgl_up]=muon.phi();
          ndgl_down++;
      }
  }
   
   
   
   // ----------------------------------
   // displacedStandAloneMuons Collection
   // ----------------------------------
   ndsa = 0;
   for (unsigned int i = 0; i < dsas->size(); i++) {
     const reco::Track& dsa(dsas->at(i));
     dsa_pt[ndsa] = dsa.pt();
     dsa_eta[ndsa] = dsa.eta();
     dsa_phi[ndsa] = dsa.phi();
     dsa_ptError[ndsa] = dsa.ptError();
     dsa_dxy[ndsa] = dsa.dxy();
     dsa_dz[ndsa] = dsa.dz();
     dsa_normalizedChi2[ndsa] = dsa.normalizedChi2();
     dsa_charge[ndsa] = dsa.charge();

     dsa_nMuonHits[ndsa] = dsa.hitPattern().numberOfMuonHits();
     dsa_nValidMuonHits[ndsa] = dsa.hitPattern().numberOfValidMuonHits();
     dsa_nValidMuonDTHits[ndsa] = dsa.hitPattern().numberOfValidMuonDTHits();
     dsa_nValidMuonCSCHits[ndsa] = dsa.hitPattern().numberOfValidMuonCSCHits();
     dsa_nValidMuonRPCHits[ndsa] = dsa.hitPattern().numberOfValidMuonRPCHits();
     dsa_nValidStripHits[ndsa] = dsa.hitPattern().numberOfValidStripHits();
     dsa_nhits[ndsa] = dsa.hitPattern().numberOfValidHits();

     dsa_nLostMuonHits[ndsa] = dsa.hitPattern().numberOfLostMuonHits();
     dsa_nLostMuonDTHits[ndsa] = dsa.hitPattern().numberOfLostMuonDTHits();
     dsa_nLostMuonCSCHits[ndsa] = dsa.hitPattern().numberOfLostMuonCSCHits();
     dsa_nLostMuonRPCHits[ndsa] = dsa.hitPattern().numberOfLostMuonRPCHits();

     dsa_dtStationsWithValidHits[ndsa] = dsa.hitPattern().dtStationsWithValidHits();
     dsa_cscStationsWithValidHits[ndsa] = dsa.hitPattern().cscStationsWithValidHits();

     ndsa++;
   }
    
        // ----------------------------------
    // displacedMuons Collection tag
    // ----------------------------------
    tag_dgl = 0;
    // Check if trigger fired:
    const edm::TriggerNames &Tag_names = iEvent.triggerNames(*triggerBits);
    unsigned int Tag_ipath = 0;
    std::vector<reco::Muon> tag_muon;  // Vector to store selected tag muons
    for (unsigned int i = 0; i < dmuons->size(); i++) {
        std::cout << " - - tag_dgl" << tag_dgl << std::endl;
        const reco::Muon& dglmuon(dmuons->at(i));
        
        if (dglmuon.isGlobalMuon()) {
            const reco::Track* globalTrack = (dglmuon.combinedMuon()).get();
            if (globalTrack == nullptr) {
              continue;  // Skip if probe global track is null
            }
            
            
            // Get the transverse momentum (pT), eta, and phi
            float pt = globalTrack->pt();
            float eta = globalTrack->eta();
            float phi = globalTrack->phi();  
            
            // Get the uncertainty in transverse momentum 
            float ptError = globalTrack->ptError();
            float relativePtError = ptError / pt;  
            
            // Apply Tag kinematic cuts:
            if (pt > 10 && phi < 0) {
                for(auto path: Tag_HLTPaths_){
                  bool fired=false;
                  for(unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg){
                    TString Tag_TrigPath = Tag_names.triggerName(itrg);
                    if(!triggerBits->accept(itrg))
                      continue;
                    if (!Tag_TrigPath.Contains(path)){
                     continue;
                   }
                   fired = true;
                 }
                 Tag_triggerPass[Tag_ipath] = fired;
                 tag_dgl_pt[tag_dgl] = pt;
                 tag_dgl_eta[tag_dgl] = eta;
                 tag_dgl_phi[tag_dgl] = phi;  // Storing the φ value
                 is_tag_dgl[tag_dgl] = dglmuon.isGlobalMuon();
                 tag_dgl_ptError[tag_dgl] = globalTrack->ptError();
                 tag_dgl_dxy[tag_dgl] = globalTrack->dxy();
                 tag_dgl_dz[tag_dgl] = globalTrack->dz();
                 tag_dgl_normalizedChi2[tag_dgl] = globalTrack->normalizedChi2();
                 tag_dgl_charge[tag_dgl] = globalTrack->charge();
                 tag_dgl_nMuonHits[tag_dgl] = globalTrack->hitPattern().numberOfMuonHits();
                 tag_dgl_nValidMuonHits[tag_dgl] = globalTrack->hitPattern().numberOfValidMuonHits();
                 tag_dgl_nValidMuonDTHits[tag_dgl] = globalTrack->hitPattern().numberOfValidMuonDTHits();
                 tag_dgl_nValidMuonCSCHits[tag_dgl] = globalTrack->hitPattern().numberOfValidMuonCSCHits();
                 tag_dgl_nValidMuonRPCHits[tag_dgl] = globalTrack->hitPattern().numberOfValidMuonRPCHits();
                 tag_dgl_nValidStripHits[tag_dgl] = globalTrack->hitPattern().numberOfValidStripHits();
                 tag_dgl_nhits[tag_dgl] = globalTrack->hitPattern().numberOfValidHits();
                 tag_muon.push_back(dglmuon);  // Store the tag muon
                 tag_dgl++;
                 Tag_ipath++;
                 } 
                
               } else {
                // If muon doesn't pass the cuts, set the variables to 0
                tag_dgl_pt[tag_dgl] = 0;
                tag_dgl_eta[tag_dgl] = 0;
                tag_dgl_phi[tag_dgl] = 0;
                is_tag_dgl[tag_dgl] = 0;
                tag_dgl_ptError[tag_dgl] = 0;
                tag_dgl_dxy[tag_dgl] = 0;
                tag_dgl_dz[tag_dgl] = 0;
                tag_dgl_normalizedChi2[tag_dgl] = 0;
                tag_dgl_charge[tag_dgl] = 0;
                tag_dgl_nMuonHits[tag_dgl] = 0;
                tag_dgl_nValidMuonHits[tag_dgl] = 0;
                tag_dgl_nValidMuonDTHits[tag_dgl] = 0;
                tag_dgl_nValidMuonCSCHits[tag_dgl] = 0;
                tag_dgl_nValidMuonRPCHits[tag_dgl] = 0;
                tag_dgl_nValidStripHits[tag_dgl] = 0;
                tag_dgl_nhits[tag_dgl] = 0;
            }//end Apply all cuts
        }//end is global 
    }//end tag selection
    
    // ----------------------------------
    // displacedMuons Collection probe
    // ----------------------------------
    probe_dgl = 0;
    
    for (unsigned int i = 0; i < tag_muon.size(); i++) {
        const reco::Muon& tag_dglmuon(tag_muon[i]);
        
        for (unsigned int l = 0; l < dmuons->size(); l++) {  
            std::cout << "- - probe selection - -" << std::endl;
            const reco::Muon& probe_dglmuon(dmuons->at(l));
            
            // Calculate pT difference between Tag and Probe muons
            //float tag_Muon_pT = tag_dglmuon.pt();
            //float probe_Muon_pt = probe_dglmuon.pt();
            
            // Check if Probe muon is back-to-back with Tag muon 
            float dPhi = fabs(deltaPhi(probe_dglmuon.phi(),tag_dglmuon.phi() ));
            
            if (dPhi > 3.14) {
                
                dPhi = 2 * 3.14 - dPhi;  // Adjust to handle wrapping around
            }
            
            // Calculate the azimuthal angle difference
            float alpha = dPhi;  // For muons nearly back-to-back
         
            const reco::Track* probe_globalTrack = (probe_dglmuon.combinedMuon()).get();
           
            if (probe_globalTrack == nullptr) {
              continue;  // Skip if probe global track is null
            }
           
            // Apply Probe cuts: back-to-back 
            if (alpha >1.0) {
                std::cout << "- - probe selection pass - -" << std::endl;
                probe_dgl_pt[probe_dgl] = probe_globalTrack->pt();
                probe_dgl_eta[probe_dgl] = probe_globalTrack->eta();
                probe_dgl_phi[probe_dgl] = probe_globalTrack->phi();
                is_probe_dgl[probe_dgl] = probe_dglmuon.isGlobalMuon();
                probe_dgl_ptError[probe_dgl] = probe_globalTrack->ptError();
                probe_dgl_dxy[probe_dgl] = probe_globalTrack->dxy();
                probe_dgl_dz[probe_dgl] = probe_globalTrack->dz();
                probe_dgl_normalizedChi2[probe_dgl] = probe_globalTrack->normalizedChi2();
                probe_dgl_charge[probe_dgl] = probe_globalTrack->charge();
                probe_dgl_nMuonHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfMuonHits();
                probe_dgl_nValidMuonHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidMuonHits();
                probe_dgl_nValidMuonDTHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidMuonDTHits();
                probe_dgl_nValidMuonCSCHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidMuonCSCHits();
                probe_dgl_nValidMuonRPCHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidMuonRPCHits();
                probe_dgl_nValidStripHits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidStripHits();
                probe_dgl_nhits[probe_dgl] = probe_globalTrack->hitPattern().numberOfValidHits();
                probe_dgl++;
            } else {
                // If probe muon doesn't pass the cuts, set the variables to 0
                probe_dgl_pt[probe_dgl] = 0;
                probe_dgl_eta[probe_dgl] = 0;
                probe_dgl_phi[probe_dgl] = 0;
                probe_dgl_ptError[probe_dgl] = 0;
                probe_dgl_dxy[probe_dgl] = 0;
                probe_dgl_dz[probe_dgl] = 0;
                probe_dgl_normalizedChi2[probe_dgl] = 0;
                probe_dgl_charge[probe_dgl] = 0;
                probe_dgl_nMuonHits[probe_dgl] = 0;
                probe_dgl_nValidMuonHits[probe_dgl] = 0;
                probe_dgl_nValidMuonDTHits[probe_dgl] = 0;
                probe_dgl_nValidMuonCSCHits[probe_dgl] = 0;
                probe_dgl_nValidMuonRPCHits[probe_dgl] = 0;
                probe_dgl_nValidStripHits[probe_dgl] = 0;
                probe_dgl_nhits[probe_dgl] = 0;
            }
        }
    }

    
    
    
    
    
    
    
    
    
    
    
    
    /*// ----------------------------------
    // displacedMuons Collection tag
    // ----------------------------------
    tag_dgl = 0;
    
    std::vector<reco::Muon> tag_muon;  // Vector to store selected tag muons
    for (unsigned int i = 0; i < dmuons->size(); i++) {
        std::cout << " - - tag_dgl" << tag_dgl << std::endl;
        const reco::Muon& dglmuon(dmuons->at(i));
        
        if (dglmuon.isGlobalMuon()) {
            const reco::Track* globalTrack = (dglmuon.combinedMuon()).get();
            if (globalTrack == nullptr) {
              continue;  // Skip if probe global track is null
            }
            
            
            // Get the transverse momentum (pT), eta, and phi
            float pt = globalTrack->pt();
            float eta = globalTrack->eta();
            float phi = globalTrack->phi();  
            
            // Get the uncertainty in transverse momentum 
            float ptError = globalTrack->ptError();
            float relativePtError = ptError / pt;  
            
            // Apply all Tag cuts:
            if (pt > 20 && fabs(eta) < 0.9 && relativePtError < 0.3 && phi > -2.6 && phi < -0.6 && globalTrack->hitPattern().numberOfValidHits() > 12 && globalTrack->hitPattern().numberOfValidStripHits() > 5) {
                tag_dgl_pt[tag_dgl] = pt;
                tag_dgl_eta[tag_dgl] = eta;
                tag_dgl_phi[tag_dgl] = phi;  // Storing the φ value
                is_tag_dgl[tag_dgl] = dglmuon.isGlobalMuon();
                tag_muon.push_back(dglmuon);  // Store the tag muon
                tag_dgl++;
            } else {
                // If muon doesn't pass the cuts, set the variables to 0
                tag_dgl_pt[tag_dgl] = 0;
                tag_dgl_eta[tag_dgl] = 0;
                tag_dgl_phi[tag_dgl] = 0;
            }//end Apply all cuts
        }//end is global 
    }//end tag selection
    
    // ----------------------------------
    // displacedMuons Collection probe
    // ----------------------------------
    probe_dgl = 0;
    
    for (unsigned int i = 0; i < tag_muon.size(); i++) {
        const reco::Muon& tag_dglmuon(tag_muon[i]);
        
        for (unsigned int l = 0; l < dmuons->size(); l++) {  
            std::cout << "- - probe selection - -" << std::endl;
            const reco::Muon& probe_dglmuon(dmuons->at(l));
            
            // Calculate pT difference between Tag and Probe muons
            //float tag_Muon_pT = tag_dglmuon.pt();
            //float probe_Muon_pt = probe_dglmuon.pt();
            
            // Check if Probe muon is back-to-back with Tag muon 
            float deltaPhi = fabs(tag_dglmuon.phi() - probe_dglmuon.phi());
            
            if (deltaPhi > 3.14) {
                
                deltaPhi = 2 * 3.14 - deltaPhi;  // Adjust to handle wrapping around
            }
            
            // Calculate the azimuthal angle difference
            float alpha = deltaPhi;  // For muons nearly back-to-back
         
            const reco::Track* probe_globalTrack = (probe_dglmuon.combinedMuon()).get();
           
            if (probe_globalTrack == nullptr) {
              continue;  // Skip if probe global track is null
            }
           
            // Apply Probe cuts: back-to-back 
            if (alpha > 2.8 && probe_globalTrack->pt() > 20) {
                std::cout << "- - probe selection pass - -" << std::endl;
                probe_dgl_pt[probe_dgl] = probe_globalTrack->pt();
                probe_dgl_eta[probe_dgl] = probe_globalTrack->eta();
                probe_dgl_phi[probe_dgl] = probe_globalTrack->phi();
                is_probe_dgl[probe_dgl] = probe_dglmuon.isGlobalMuon();
                probe_dgl++;
            } else {
                // If probe muon doesn't pass the cuts, set the variables to 0
                probe_dgl_pt[probe_dgl] = 0;
                probe_dgl_eta[probe_dgl] = 0;
                probe_dgl_phi[probe_dgl] = 0;
            }
        }
    }*/

   
   
   // ----------------------------------
   // displacedMuons Collection
   // ----------------------------------
   ndmu = 0;
   for (unsigned int i = 0; i < dmuons->size(); i++) {
     std::cout << " - - ndmu: " << ndmu << std::endl;
     const reco::Muon& dmuon(dmuons->at(i));
     dmu_isDGL[ndmu] = dmuon.isGlobalMuon();
     dmu_isDSA[ndmu] = dmuon.isStandAloneMuon();
     dmu_isDTK[ndmu] = dmuon.isTrackerMuon();
     dmu_isMatchesValid[ndmu] = dmuon.isMatchesValid();
     dmu_numberOfMatches[ndmu] = dmuon.numberOfMatches();
     dmu_numberOfChambers[ndmu] = dmuon.numberOfChambers();
     dmu_numberOfChambersCSCorDT[ndmu] = dmuon.numberOfChambersCSCorDT();
     dmu_numberOfMatchedStations[ndmu] = dmuon.numberOfMatchedStations();
     dmu_numberOfMatchedRPCLayers[ndmu] = dmuon.numberOfMatchedRPCLayers();

     // Access the DGL track associated to the displacedMuon
     std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
     if ( dmuon.isGlobalMuon() ) {
       const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
       dmu_dgl_pt[ndmu] = globalTrack->pt();
       dmu_dgl_eta[ndmu] = globalTrack->eta();
       dmu_dgl_phi[ndmu] = globalTrack->phi();
       dmu_dgl_ptError[ndmu] = globalTrack->ptError();
       dmu_dgl_dxy[ndmu] = globalTrack->dxy();
       dmu_dgl_dz[ndmu] = globalTrack->dz();
       dmu_dgl_normalizedChi2[ndmu] = globalTrack->normalizedChi2();
       dmu_dgl_charge[ndmu] = globalTrack->charge();
       dmu_dgl_nMuonHits[ndmu] = globalTrack->hitPattern().numberOfMuonHits();
       dmu_dgl_nValidMuonHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonHits();
       dmu_dgl_nValidMuonDTHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dgl_nValidMuonCSCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dgl_nValidMuonRPCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dgl_nValidStripHits[ndmu] = globalTrack->hitPattern().numberOfValidStripHits();
       dmu_dgl_nhits[ndmu] = globalTrack->hitPattern().numberOfValidHits();
       
       

     } else {
       dmu_dgl_pt[ndmu] = 0;
       dmu_dgl_eta[ndmu] = 0;
       dmu_dgl_phi[ndmu] = 0;
       dmu_dgl_ptError[ndmu] = 0;
       dmu_dgl_dxy[ndmu] = 0;
       dmu_dgl_dz[ndmu] = 0;
       dmu_dgl_normalizedChi2[ndmu] = 0;
       dmu_dgl_charge[ndmu] = 0;
       dmu_dgl_nMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonDTHits[ndmu] = 0;
       dmu_dgl_nValidMuonCSCHits[ndmu] = 0;
       dmu_dgl_nValidMuonRPCHits[ndmu] = 0;
       dmu_dgl_nValidStripHits[ndmu] = 0;
       dmu_dgl_nhits[ndmu] = 0;
     }     

     // Access the DSA track associated to the displacedMuon
     std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
     if ( dmuon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
       dmu_dsa_pt[ndmu] = outerTrack->pt();
       dmu_dsa_eta[ndmu] = outerTrack->eta();
       dmu_dsa_phi[ndmu] = outerTrack->phi();
       dmu_dsa_ptError[ndmu] = outerTrack->ptError();
       dmu_dsa_dxy[ndmu] = outerTrack->dxy();
       dmu_dsa_dz[ndmu] = outerTrack->dz();
       dmu_dsa_normalizedChi2[ndmu] = outerTrack->normalizedChi2();
       dmu_dsa_charge[ndmu] = outerTrack->charge();
       dmu_dsa_nMuonHits[ndmu] = outerTrack->hitPattern().numberOfMuonHits();
       dmu_dsa_nValidMuonHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonHits();
       dmu_dsa_nValidMuonDTHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dsa_nValidMuonCSCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dsa_nValidMuonRPCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dsa_nValidStripHits[ndmu] = outerTrack->hitPattern().numberOfValidStripHits();
       dmu_dsa_nhits[ndmu] = outerTrack->hitPattern().numberOfValidHits();
       dmu_dsa_dtStationsWithValidHits[ndmu] = outerTrack->hitPattern().dtStationsWithValidHits();
       dmu_dsa_cscStationsWithValidHits[ndmu] = outerTrack->hitPattern().cscStationsWithValidHits();
       if (isAOD) {
         // Number of DT+CSC segments
         unsigned int nsegments = 0;
         for (trackingRecHit_iterator hit = outerTrack->recHitsBegin(); hit != outerTrack->recHitsEnd(); ++hit) {
           if (!(*hit)->isValid()) continue;
           DetId id = (*hit)->geographicalId();
           if (id.det() != DetId::Muon) continue;
           if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
             nsegments++;
           }
         }
         dmu_dsa_nsegments[ndmu] = nsegments;
       }
     } else {
       dmu_dsa_pt[ndmu] = 0;
       dmu_dsa_eta[ndmu] = 0;
       dmu_dsa_phi[ndmu] = 0;
       dmu_dsa_ptError[ndmu] = 0;
       dmu_dsa_dxy[ndmu] = 0;
       dmu_dsa_dz[ndmu] = 0;
       dmu_dsa_normalizedChi2[ndmu] = 0;
       dmu_dsa_charge[ndmu] = 0;
       dmu_dsa_nMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonDTHits[ndmu] = 0;
       dmu_dsa_nValidMuonCSCHits[ndmu] = 0;
       dmu_dsa_nValidMuonRPCHits[ndmu] = 0;
       dmu_dsa_nValidStripHits[ndmu] = 0;
       dmu_dsa_nhits[ndmu] = 0;
       dmu_dsa_dtStationsWithValidHits[ndsa] = 0;
       dmu_dsa_cscStationsWithValidHits[ndsa] = 0;
       dmu_dsa_nsegments[ndmu] = 0;
     }

     // Access the DTK track associated to the displacedMuon
     std::cout << "isTrackerMuon: " << dmuon.isTrackerMuon() << std::endl;
     if ( dmuon.isTrackerMuon() ) {
       const reco::Track* innerTrack = (dmuon.track()).get();
       dmu_dtk_pt[ndmu] = innerTrack->pt();
       dmu_dtk_eta[ndmu] = innerTrack->eta();
       dmu_dtk_phi[ndmu] = innerTrack->phi();
       dmu_dtk_ptError[ndmu] = innerTrack->ptError();
       dmu_dtk_dxy[ndmu] = innerTrack->dxy();
       dmu_dtk_dz[ndmu] = innerTrack->dz();
       dmu_dtk_normalizedChi2[ndmu] = innerTrack->normalizedChi2();
       dmu_dtk_charge[ndmu] = innerTrack->charge();
       dmu_dtk_nMuonHits[ndmu] = innerTrack->hitPattern().numberOfMuonHits();
       dmu_dtk_nValidMuonHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonHits();
       dmu_dtk_nValidMuonDTHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dtk_nValidMuonCSCHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dtk_nValidMuonRPCHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dtk_nValidStripHits[ndmu] = innerTrack->hitPattern().numberOfValidStripHits();
       dmu_dtk_nhits[ndmu] = innerTrack->hitPattern().numberOfValidHits();
     } else {
       dmu_dtk_pt[ndmu] = 0;
       dmu_dtk_eta[ndmu] = 0;
       dmu_dtk_phi[ndmu] = 0;
       dmu_dtk_ptError[ndmu] = 0;
       dmu_dtk_dxy[ndmu] = 0;
       dmu_dtk_dz[ndmu] = 0;
       dmu_dtk_normalizedChi2[ndmu] = 0;
       dmu_dtk_charge[ndmu] = 0;
       dmu_dtk_nMuonHits[ndmu] = 0;
       dmu_dtk_nValidMuonHits[ndmu] = 0;
       dmu_dtk_nValidMuonDTHits[ndmu] = 0;
       dmu_dtk_nValidMuonCSCHits[ndmu] = 0;
       dmu_dtk_nValidMuonRPCHits[ndmu] = 0;
       dmu_dtk_nValidStripHits[ndmu] = 0;
       dmu_dtk_nhits[ndmu] = 0;
     }

     ndmu++;
     std::cout << "End muon" << std::endl;
   }

   if (!isData) {
     // ----------------------------------
     // PrimaryVertices collection
     // ----------------------------------
     nTruePV = 0;
     nPV = primaryvertices->size();

     for (size_t i = 0; i < primaryvertices->size(); i ++){
         const reco::Vertex &current_vertex = (*primaryvertices)[i];
         if(current_vertex.isValid()){ nTruePV++; }
     }

     // The PV information:
     const reco::Vertex &thePrimaryVertex = (*primaryvertices)[0];
     
     PV_vx = thePrimaryVertex.x();
     PV_vy = thePrimaryVertex.y();
     PV_vz = thePrimaryVertex.z();
     PV_passAcceptance = false;
     
     if (!thePrimaryVertex.isFake() && thePrimaryVertex.ndof() > 4 && fabs(thePrimaryVertex.z()) < 25 && thePrimaryVertex.position().rho() <= 2) {
        PV_passAcceptance = true;
     }   
     GlobalPoint _PVpoint(thePrimaryVertex.x(), thePrimaryVertex.y(), thePrimaryVertex.z());

     // ----------------------------------
     // GenParticle collection
     // ----------------------------------
     std::vector<int> iGM; // gen muons count

     reco::GenParticleRef mref;
     reco::GenParticle m;
     // Get the muons that can be reconstructed
     for (size_t i = 0; i < genParticles->size(); i++) {
       const reco::GenParticle &genparticle = (*genParticles)[i];
       if ( abs(genparticle.pdgId()) == 13 && genparticle.status() == 1) {
         iGM.push_back(i);
       }
     }
     nGenMuon = iGM.size();
     std::cout << "Number of gen muons = " << nGenMuon << std::endl;

     for (size_t i = 0; i < iGM.size(); i++) {
       const reco::GenParticle &genparticle = (*genParticles)[iGM.at(i)];
         
       GenLeptonSel_pt[i] = genparticle.pt();
       GenLeptonSel_E[i] = genparticle.energy();
       GenLeptonSel_et[i] = genparticle.et();
       GenLeptonSel_eta[i] = genparticle.eta();
       GenLeptonSel_phi[i] = genparticle.phi();
       GenLeptonSel_pdgId[i] = genparticle.pdgId();

       // Bottom-up to get the real decaying particle:
       if (genparticle.mother()->pdgId() == genparticle.pdgId()) {
         mref = genparticle.motherRef();
         m = *mref;
         while (m.pdgId() == m.mother()->pdgId()) {
           mref = m.motherRef();
           m = *mref;
         }

         GenLeptonSel_vx[i] = m.vx();
         GenLeptonSel_vy[i] = m.vy();
         GenLeptonSel_vz[i] = m.vz();
         GenLeptonSel_dxy[i] = dxy_value(m, thePrimaryVertex); // should be computed here or before?

         if (m.numberOfMothers() != 0) {
           GenLeptonSel_motherPdgId[i] = m.motherRef()->pdgId();
         } else {
           GenLeptonSel_motherPdgId[i] = 0; 
         }
       } else {
         GenLeptonSel_vx[i] = genparticle.vx();
         GenLeptonSel_vy[i] = genparticle.vy();
         GenLeptonSel_vz[i] = genparticle.vz();
         GenLeptonSel_dxy[i] = dxy_value(genparticle, thePrimaryVertex); // should be computed here or before?

         GenLeptonSel_motherPdgId[i] = genparticle.motherRef()->pdgId();
       }
            
       // Flags
       GenLeptonSel_isPromptFinalState[i] = genparticle.isPromptFinalState();
       GenLeptonSel_fromHardProcessFinalState[i] = genparticle.fromHardProcessFinalState(); // has to be done with the last one
       GenLeptonSel_isDirectPromptTauDecayProductFinalState[i] = genparticle.isDirectPromptTauDecayProductFinalState(); 
       GenLeptonSel_isDirectHadronDecayProduct[i] = genparticle.statusFlags().isDirectHadronDecayProduct(); 
     }
         
     // Counters initialization
     nGenMuon_PFS = 0; 
     nGenMuon_HPFS = 0; 
     nGenMuon_PTDP = 0; 
     nGenMuon_HDP = 0; 
     for (size_t i = 0; i < iGM.size(); i++) {
       if (GenLeptonSel_isPromptFinalState[i]) { nGenMuon_PFS++; }
       if (GenLeptonSel_fromHardProcessFinalState[i]) { nGenMuon_HPFS++; }
       if (GenLeptonSel_isDirectPromptTauDecayProductFinalState[i]) { nGenMuon_PTDP++; }
       if (GenLeptonSel_isDirectHadronDecayProduct[i]) { nGenMuon_HDP++; }
     }

     // ----------------------------------
     // Hard Process Collection
     // ----------------------------------
     nHardProcessParticle = 0;
     for (size_t i = 0; i < genParticles->size(); i++) {
       const reco::GenParticle &genparticle = (*genParticles)[i];
       if (genparticle.isHardProcess()){
         HardProcessParticle_pt[nHardProcessParticle] = genparticle.pt();
         HardProcessParticle_E[nHardProcessParticle] = genparticle.energy();
         HardProcessParticle_eta[nHardProcessParticle] = genparticle.eta();
         HardProcessParticle_phi[nHardProcessParticle] = genparticle.phi();
         HardProcessParticle_vx[nHardProcessParticle] = genparticle.vx();
         HardProcessParticle_vy[nHardProcessParticle] = genparticle.vy();
         HardProcessParticle_vz[nHardProcessParticle] = genparticle.vz();
         HardProcessParticle_pdgId[nHardProcessParticle] = genparticle.pdgId();
         nHardProcessParticle++;
       }       
     }
   }

   // Check if trigger fired:
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   unsigned int ipath = 0;
   for (auto path : HLTPaths_) {
     std::string path_v = path + "_v";
     // std::cout << path << "\t" << std::endl;
     bool fired = false;
     for (unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
       TString TrigPath = names.triggerName(itrg);
       if (!triggerBits->accept(itrg))
         continue;
       if (!TrigPath.Contains(path_v)){
         continue;
       }
       fired = true;
     }
     triggerPass[ipath] = fired;
     ipath++;
   } 

   //-> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(ntuplizer);
