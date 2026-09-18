#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <deque>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "configuration.h"
#include "trees.h"
#include "utilities.h"

struct GenericBranchBuffer {
  std::string name;
  std::string type; // "int", "float", "double"
  int val_i = 0;
  float val_f = 0.0f;
  double val_d = 0.0;

  void bind_input(TTree *tree, TBranch *branch) {
    name = branch->GetName();
    TLeaf *leaf = (TLeaf *)branch->GetListOfLeaves()->At(0);
    TString tname = leaf ? leaf->GetTypeName() : "Double_t";
    if (tname == "Int_t" || tname == "int" || tname == "Short_t" ||
        tname == "UShort_t" || tname == "UInt_t" || tname == "Long64_t" ||
        tname == "ULong64_t" || tname == "Bool_t" || tname == "Char_t" ||
        tname == "UChar_t") {
      type = "int";
      tree->SetBranchAddress(name.c_str(), &val_i);
    } else if (tname == "Float_t" || tname == "float") {
      type = "float";
      tree->SetBranchAddress(name.c_str(), &val_f);
    } else {
      type = "double";
      tree->SetBranchAddress(name.c_str(), &val_d);
    }
  }

  void bind_output(TTree *tree) {
    if (type == "int") {
      tree->Branch(name.c_str(), &val_i, (name + "/I").c_str());
    } else if (type == "float") {
      tree->Branch(name.c_str(), &val_f, (name + "/F").c_str());
    } else {
      tree->Branch(name.c_str(), &val_d, (name + "/D").c_str());
    }
  }

  double get_as_double() const {
    if (type == "int")
      return (double)val_i;
    if (type == "float")
      return (double)val_f;
    return val_d;
  }
};

void copy_no_syst(cfg::ConfigurationTable table, TTree *out_tree,
                  TTree *in_tree);
void copy_with_syst(cfg::ConfigurationTable config,
                    cfg::ConfigurationTable table, TTree *out_tree,
                    TTree *in_tree, TTree *syst_in_tree, std::string syst_type);

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <config_file.toml> [input_file.root]"
              << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////
  /// Config
  /////////////////////////////////////////////////////////////
  cfg::ConfigurationTable config;
  config.set_config(argv[1]);
  std::vector<cfg::ConfigurationTable> tables;
  tables = config.get_subtables("tree");

  /////////////////////////////////////////////////////////////
  /// Input
  /////////////////////////////////////////////////////////////
  std::string input_filename;
  if (argc > 2)
    input_filename = argv[2];
  else
    input_filename = config.get_string_field("output.path");

  TFile *input = TFile::Open(input_filename.c_str(), "READ");
  if (!input || input->IsZombie()) {
    std::cerr << "Error: Could not open input ROOT file " << input_filename
              << std::endl;
    return 1;
  }

  /// Syst trees
  TTree *multisigma_tree =
      (TTree *)input->Get("events/cvext/selected_multisigmaTree");
  TTree *variation_tree =
      (TTree *)input->Get("events/cvext/selected_variationTree");
  TTree *multisim_tree =
      (TTree *)input->Get("events/full/selected_multisimTree");
  TTree *NuMIflux_tree =
      (TTree *)input->Get("events/cvext/selected_NuMIfluxsimTree");

  /////////////////////////////////////////////////////////////
  /// Output
  /////////////////////////////////////////////////////////////
  TFile *output = TFile::Open(
      "icarus_numi_numu_mc_onbeam_offbeam_syst_gundam.root", "RECREATE");
  if (!output || output->IsZombie()) {
    std::cerr << "Error: Could not create output ROOT file." << std::endl;
    input->Close();
    return 1;
  }

  /////////////////////////////////////////////////////////////
  // Loop over every tree in TOML config
  /////////////////////////////////////////////////////////////
  for (cfg::ConfigurationTable &table : tables) {
    // Grab input TTree
    TTree *in_tree =
        (TTree *)input->Get(table.get_string_field("origin").c_str());
    if (!in_tree) {
      std::cerr << "Warning: Tree " << table.get_string_field("origin")
                << " not found in input file." << std::endl;
      continue;
    }

    std::cout << "Processing " << table.get_string_field("origin") << std::endl;

    // Create output TTree
    TDirectory *directory = (TDirectory *)output;
    directory =
        create_directory(directory, table.get_string_field("destination"));
    directory->cd();

    TDirectory *parent = (TDirectory *)input;
    parent =
        get_parent_directory(parent, table.get_string_field("destination"));
    if (parent) {
      TH1D *pot = (TH1D *)parent->Get("POT");
      TH1D *livetime = (TH1D *)parent->Get("Livetime");
      if (pot)
        directory->WriteObject(pot, "POT");
      if (livetime)
        directory->WriteObject(livetime, "Livetime");
    }

    TTree *out_tree = new TTree(table.get_string_field("name").c_str(),
                                table.get_string_field("name").c_str());
    copy_no_syst(table, out_tree, in_tree);
    if (table.get_bool_field("gundam_store_syst") == true) {
      if (multisigma_tree)
        copy_with_syst(config, table, out_tree, in_tree, multisigma_tree,
                       "multisigma");
      if (variation_tree)
        copy_with_syst(config, table, out_tree, in_tree, variation_tree,
                       "variation");
      if (NuMIflux_tree)
        copy_with_syst(config, table, out_tree, in_tree, NuMIflux_tree,
                       "NuMIfluxsim");
    }
    out_tree->Write();
  }

  if (multisim_tree) {
    output->cd();
    TTree *multisim_tree_copy = multisim_tree->CloneTree(-1);
    if (multisim_tree_copy)
      multisim_tree_copy->Write();
  }

  output->Close();
  input->Close();

  return 0;
}

void copy_no_syst(cfg::ConfigurationTable table, TTree *out_tree,
                  TTree *in_tree) {
  if (!in_tree || !out_tree) {
    std::cerr << "Error: in_tree or out_tree is null!" << std::endl;
    return;
  }

  TObjArray *branches = in_tree->GetListOfBranches();
  if (!branches)
    return;

  int n_branches = branches->GetEntries();
  std::cout << "Number of branches: " << n_branches << std::endl;

  // Variables for special branches
  int run(0), subrun(0), event(0);
  int cut_type(0), is_nu(0), is_data(0), category(0);

  GenericBranchBuffer run_buf, subrun_buf, evt_buf;
  GenericBranchBuffer cut_type_buf, is_nu_buf, is_data_buf, category_buf;

  bool has_run = false, has_subrun = false, has_evt = false;
  bool has_cut_type = false, has_is_nu = false, has_is_data = false,
       has_category = false;

  // deque: buffer addresses are bound to the trees, so elements must not move
  std::deque<GenericBranchBuffer> other_branches;

  for (int i = 0; i < n_branches; ++i) {
    TBranch *branch = (TBranch *)branches->At(i);
    std::string bname = branch->GetName();

    if (bname == "Run" || bname == "run") {
      run_buf.bind_input(in_tree, branch);
      has_run = true;
    } else if (bname == "Subrun" || bname == "subrun") {
      subrun_buf.bind_input(in_tree, branch);
      has_subrun = true;
    } else if (bname == "Evt" || bname == "evt" || bname == "Event" ||
               bname == "event") {
      evt_buf.bind_input(in_tree, branch);
      has_evt = true;
    } else if (bname == "reco_cut_type" || bname == "cut_type") {
      cut_type_buf.bind_input(in_tree, branch);
      has_cut_type = true;
    } else if (bname == "reco_is_nu" || bname == "is_nu") {
      is_nu_buf.bind_input(in_tree, branch);
      has_is_nu = true;
    } else if (bname == "reco_is_data" || bname == "is_data") {
      is_data_buf.bind_input(in_tree, branch);
      has_is_data = true;
    } else if (bname == "true_category" || bname == "category") {
      category_buf.bind_input(in_tree, branch);
      has_category = true;
    } else {
      other_branches.emplace_back();
      GenericBranchBuffer& gbuf = other_branches.back();
      gbuf.bind_input(in_tree, branch);
      gbuf.bind_output(out_tree);
    }
  }

  // Create standard required output branches
  out_tree->Branch("Run", &run, "Run/I");
  out_tree->Branch("Subrun", &subrun, "Subrun/I");
  out_tree->Branch("Evt", &event, "Evt/I");
  out_tree->Branch("cut_type", &cut_type, "cut_type/I");
  out_tree->Branch("is_nu", &is_nu, "is_nu/I");
  out_tree->Branch("is_data", &is_data, "is_data/I");
  out_tree->Branch("category", &category, "category/I");

  // Copy entries
  Long64_t nentries = in_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    in_tree->GetEntry(i);

    if (has_run)
      run = (int)run_buf.get_as_double();
    if (has_subrun)
      subrun = (int)subrun_buf.get_as_double();
    if (has_evt)
      event = (int)evt_buf.get_as_double();

    if (has_cut_type)
      cut_type = (int)cut_type_buf.get_as_double();
    else
      cut_type = 0;

    if (table.get_bool_field("is_nu"))
      is_nu = 1;
    else if (has_is_nu)
      is_nu = (int)is_nu_buf.get_as_double();
    else
      is_nu = 0;

    if (table.get_bool_field("is_data")) {
      is_data = 1;
      category = 10;
    } else {
      if (has_is_data)
        is_data = (int)is_data_buf.get_as_double();
      else
        is_data = 0;

      if (has_category)
        category = (int)category_buf.get_as_double();
      else
        category = 0;
    }

    out_tree->Fill();
  }
}

void copy_with_syst(cfg::ConfigurationTable config,
                    cfg::ConfigurationTable table, TTree *out_tree,
                    TTree *in_tree, TTree *syst_in_tree,
                    std::string syst_type) {
  if (!syst_in_tree || !out_tree) {
    std::cerr << "CRITICAL ERROR: syst_in_tree or out_tree is null!"
              << std::endl;
    return;
  }

  // Get relevant syst branches from TOML
  int num_syst_branches(0);
  std::vector<std::string> syst_branch_names;

  for (cfg::ConfigurationTable &t : config.get_subtables("sys")) {
    if (!strcmp(t.get_string_field("type").c_str(), syst_type.c_str())) {
      syst_branch_names.push_back(t.get_string_field("name"));
      num_syst_branches++;
    }
  }

  if (num_syst_branches == 0) {
    std::cout << "No systematic branches found matching type: " << syst_type
              << std::endl;
    return;
  }

  // Initialize TClonesArray
  std::vector<TClonesArray *> arrSyst(num_syst_branches, nullptr);

  // Get the list of branches from input syst tree
  TObjArray *branches = syst_in_tree->GetListOfBranches();
  if (!branches)
    return;

  // Create a map to store the addresses of the vectors
  std::map<std::string, std::vector<Float_t> *> branch_addresses_nsigmas;
  std::map<std::string, std::vector<Float_t> *> branch_addresses_weights;
  std::map<std::string, int> sys_idxs;
  std::vector<TString> relevant_branch_names;

  // Loop through each branch and set its address
  int sysIdx(0);
  for (int i = 0; i < branches->GetEntries(); ++i) {
    TBranch *branch = (TBranch *)branches->At(i);
    TString branch_name = branch->GetName();

    if (std::find(syst_branch_names.begin(), syst_branch_names.end(),
                  std::string(branch_name.Data())) != syst_branch_names.end()) {
      sys_idxs[branch_name.Data()] = sysIdx;
      relevant_branch_names.push_back(branch_name);

      std::vector<Float_t> *vec2 = nullptr;
      branch_addresses_weights[branch_name.Data()] = vec2;
      syst_in_tree->SetBranchAddress(
          branch_name, &(branch_addresses_weights[branch_name.Data()]));

      TString temp(branch_name);
      TString dummy("_sigma");
      TString branch_name_nsigma = temp + dummy;

      std::vector<Float_t> *vec1 = nullptr;
      branch_addresses_nsigmas[branch_name.Data()] = vec1;
      syst_in_tree->SetBranchAddress(
          branch_name_nsigma, &(branch_addresses_nsigmas[branch_name.Data()]));

      arrSyst[sysIdx] = new TClonesArray("TGraph", 1);
      out_tree->Branch(branch_name, &arrSyst[sysIdx], 32000, -1);
      sysIdx++;
    }
  }

  // Define morph dial names for special handling
  std::vector<std::string> morph_names = {
      "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
      "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
      "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
      "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"};

  // Loop through the entries of the tree and access the data
  Long64_t nEntries = syst_in_tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    syst_in_tree->GetEntry(i);

    // Loop over each systematic parameter
    for (auto const &[name, y_vec_ptr] : branch_addresses_weights) {
      std::vector<Float_t> *x_vec_ptr = branch_addresses_nsigmas.at(name);
      std::vector<Float_t> nsigmas;
      std::vector<Float_t> weights;

      // Check if this is a morph dial
      bool is_morph = (!strcmp(syst_type.c_str(), "multisigma")) &&
                      (std::find(morph_names.begin(), morph_names.end(),
                                 name) != morph_names.end());

      if (is_morph) {
        // Build 5-point symmetric spline for morph dials
        Float_t w1 =
            (y_vec_ptr && !y_vec_ptr->empty()) ? y_vec_ptr->front() : 1.f;
        if (!std::isfinite(w1))
          w1 = 1.f;
        if (table.get_bool_field("is_nu") == false)
          w1 = 1.f;

        Float_t w05 = 1.f + 0.5f * (w1 - 1.f);

        nsigmas.clear();
        weights.clear();

        nsigmas.push_back(-1.f);
        nsigmas.push_back(-0.5f);
        nsigmas.push_back(0.f);
        nsigmas.push_back(0.5f);
        nsigmas.push_back(1.f);

        weights.push_back(w1);
        weights.push_back(w05);
        weights.push_back(1.f);
        weights.push_back(w05);
        weights.push_back(1.f);
      } else {
        // Original behavior for standard multisigma dials
        if (x_vec_ptr) {
          for (Float_t val : *x_vec_ptr)
            nsigmas.push_back(val);
        }

        if (y_vec_ptr) {
          for (Float_t val : *y_vec_ptr) {
            if (table.get_bool_field("is_nu") == false)
              weights.push_back(-5);
            else
              weights.push_back(val);
          }
        }

        // weight = 1 for nsigma = 0
        if (!strcmp(syst_type.c_str(), "multisigma")) {
          nsigmas.push_back(0);
          if (table.get_bool_field("is_nu") == false)
            weights.push_back(-5);
          else
            weights.push_back(1);
        }
      }

      // Create a TGraph for every event
      if (!nsigmas.empty() && !weights.empty()) {
        TGraph graph(nsigmas.size(), &nsigmas[0], &weights[0]);
        graph.Sort();
        new ((*arrSyst[sys_idxs[name]])[0])
            TGraph(graph.GetN(), graph.GetX(), graph.GetY());
      }
    } // end systematic parameter loop

    // Fill specified systematic branches
    TObjArray *out_branches = out_tree->GetListOfBranches();
    if (out_branches) {
      for (int j = 0; j < out_branches->GetEntries(); ++j) {
        TBranch *branch = (TBranch *)out_branches->At(j);
        TString branch_name = branch->GetName();
        if (std::find(relevant_branch_names.begin(),
                      relevant_branch_names.end(),
                      branch_name) != relevant_branch_names.end()) {
          branch->Fill();
        }
      }
    }
  } // end entry loop

} // end function
