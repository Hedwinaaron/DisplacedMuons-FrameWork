#include <sys/stat.h>
#include <iostream>
#include <filesystem>
#include <TROOT.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <string>
#include <cmath> // For std::abs
#include <vector> 
namespace fs = std::filesystem;

void Dmuons_plots_tnp() {
    // Path to the directory containing the ROOT files
    std::string path = "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2022F_MiniAOD-Ntuples_TnP_CMSSW_13_2_0_pre1_9_23_test2/241029_185245/0000";
    // Vector to store the file paths
    std::vector<std::string> files;

    // Check if the provided path is a directory or a single file
    if (fs::is_directory(path)) {
        // Loop through the directory and collect all ROOT file paths
        for (const auto& entry : fs::directory_iterator(path)) {
            if (entry.path().extension() == ".root") {
                files.push_back(entry.path());
                std::cout << entry.path() << std::endl;
            }
        }
    } else if (fs::is_regular_file(path) && fs::path(path).extension() == ".root") {
        // If the path is a single ROOT file, add it to the list
        files.push_back(path);
        std::cout << path << std::endl;
    } else {
        std::cerr << "Invalid path: " << path << std::endl;
        return;
    }
    // Create a histograms
    TH1F* dSamu_pt_his = new TH1F("dSamu_pt_his", "dSamu_pt_his", 100, 0, 100);
    TH1F* dSamu_pt_his_probe = new TH1F("dSamu_pt_probe", "dSamu_pt_probe", 100, 0, 100);

    std::vector<int> Tag_Probe;
    std::vector<int> only_tag;
    std::vector<int> pT;
    // Loop through each file
    for (const auto& file : files) {

        std::cout << "Reading file: " << file << std::endl;

        // Open the ROOT file
        TFile* f = TFile::Open(file.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << file << std::endl;
            continue;
        }
        // Initialize a TTreeReader for the tree
        TTreeReader myReader("Events", f);

        // TTreeReaderArray to read the branches
        TTreeReaderArray<Int_t> isDGL(myReader, "dmu_isDGL");
        TTreeReaderArray<Int_t> isDGL_tag(myReader, "is_tag_dgl");
        TTreeReaderArray<Int_t> isDGL_probe(myReader, "is_probe_dgl");
        //pT
        TTreeReaderArray<Float_t> mupt_dGL_tag(myReader, "tag_dgl_pt");
        TTreeReaderArray<Float_t> mupt_dGL_probe(myReader, "probe_dgl_pt");
        //Eta
        TTreeReaderArray<Float_t> mueta_dGL_tag(myReader, "tag_dgl_eta");
        TTreeReaderArray<Float_t> mueta_dGL_probe(myReader, "probe_dgl_eta");
        //Phi
        TTreeReaderArray<Float_t> muphi_dGL_tag(myReader, "tag_dgl_phi");
        TTreeReaderArray<Float_t> muphi_dGL_probe(myReader, "probe_dgl_phi");
        
        //trigger match
        TTreeReaderArray<bool> HLT_L2Mu10_NoVertex_NoBPTX3BX(myReader, "HLT_L2Mu10_NoVertex_NoBPTX3BX");
        TTreeReaderArray<bool> HLT_L2Mu10_NoVertex_NoBPTX(myReader, "HLT_L2Mu10_NoVertex_NoBPTX");

     
        // Loop over all entries in the tree
        while (myReader.Next()) {
            
           
            //std::cout<<"/////////////////////////////////////Vegining of the event///////////////////////////////////////////"<<std::endl;


            


            // Loop through the elements of the arrays
            
            for (size_t i = 0; i < isDGL.GetSize(); ++i) {
                // Fill the histogram
                if(isDGL_tag[i]==1 && isDGL_probe[i]==1 && HLT_L2Mu10_NoVertex_NoBPTX3BX[i]==1 && HLT_L2Mu10_NoVertex_NoBPTX[i]==1)
                dSamu_pt_his->Fill(mupt_dGL_tag[i]);
                dSamu_pt_his_probe->Fill(mupt_dGL_probe[i]);
                   
                    
            }



        }

        // Close the ROOT file
        f->Close();
        delete f;

    }

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
    dSamu_pt_his->Draw();
    c1->SaveAs("dmu_dgl_pt_tag.png");
    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 800);
    dSamu_pt_his_probe->Draw();
    c2->SaveAs("dmu_dgl_pt_probe.png");

    // Clean up
    delete c1;
    delete dSamu_pt_his;
}
