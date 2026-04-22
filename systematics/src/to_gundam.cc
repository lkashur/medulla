#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraph.h>
#include <vector>
#include "configuration.h"
#include "trees.h"

void copy_no_syst(cfg::ConfigurationTable table, TTree * out_tree, TTree * in_tree);
void copy_with_syst(cfg::ConfigurationTable config, cfg::ConfigurationTable table, TTree * out_tree, TTree * in_tree, TTree* syst_in_tree, std::string syst_type, std::string sel_or_sig);

int main(int argc, char * argv[])
{

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
  input_filename = argv[2];
  TFile *input = TFile::Open(input_filename.c_str(), "READ");

  /////////////////////////////////////////////////////////////
  /// Test
  /////////////////////////////////////////////////////////////
  std::string sel_or_sig;
  sel_or_sig = argv[3];
  
  TTree* multisigma_tree = nullptr;
  TTree* variation_tree = nullptr;
  TTree* multisim_tree = nullptr;
  
  if(sel_or_sig == "sel")
  {
    input->GetObject("events/cvfull/selected_nu_multisigmaTree", multisigma_tree);
    input->GetObject("events/cvfull/selected_nu_variationTree", variation_tree);
    input->GetObject("events/cvfull/selected_nu_multisimTree",multisim_tree);
  }
  else
  {
    input->GetObject("events/cvfull/signal_multisigmaTree", multisigma_tree);
    input->GetObject("events/cvfull/signal_variationTree", variation_tree);
    input->GetObject("events/cvfull/signal_multisimTree", multisim_tree);
  }

  /// Syst trees
  /// To-do: Automate this
  //TTree* multisigma_tree = (TTree*)input->Get("events/cvfull/selected_nu_multisigmaTree");
  //TTree* variation_tree = (TTree*)input->Get("events/cvfull/selected_nu_variationTree");
  //TTree* multisim_tree = (TTree*)input->Get("events/cvfull/selected_nu_multisimTree"); // TEST...delete if anything breaks
  
  //TTree* multisigma_tree = (TTree*)input->Get("events/cvfull/signal_multisigmaTree");
  //TTree* variation_tree = (TTree*)input->Get("events/cvfull/signal_variationTree");
  //TTree* multisim_tree = (TTree*)input->Get("events/cvfull/signal_multisimTree");

  /////////////////////////////////////////////////////////////
  /// Output
  /////////////////////////////////////////////////////////////
  TFile* output = TFile::Open("output.root", "RECREATE");

  /////////////////////////////////////////////////////////////
  // Loop over every tree in TOML config
  /////////////////////////////////////////////////////////////
  for(cfg::ConfigurationTable & table : tables)
    {
      // Grab input TTree
      TTree* in_tree = (TTree*)input->Get(table.get_string_field("origin").c_str());                                                                                                  
      
      std::cout << "Processing " <<table.get_string_field("origin").c_str() << std::endl;
                   
      // Create output TTree
      TDirectory * directory = (TDirectory *) output;
      directory = create_directory(directory, table.get_string_field("destination").c_str());
      directory->cd();
      TDirectory * parent = (TDirectory *) input;
      parent = get_parent_directory(parent, table.get_string_field("destination").c_str());
      TH1D * pot = (TH1D *) parent->Get("POT");
      TH1D * livetime = (TH1D *) parent->Get("Livetime");
      directory->WriteObject(pot, "POT");
      directory->WriteObject(livetime, "Livetime");
      TTree* out_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());
      copy_no_syst(table, out_tree, in_tree);
      
      if(table.get_bool_field("gundam_store_syst") == true)
        {
          copy_with_syst(config, table, out_tree, in_tree, multisigma_tree, "multisigma", sel_or_sig);
          copy_with_syst(config, table, out_tree, in_tree, variation_tree, "variation", sel_or_sig);
        }
      
      out_tree->Write();
    }

  // TEST: Copy multisim tree as-is...delete if anything breaks 
  output->cd();
  TTree* multisim_tree_copy = multisim_tree->CloneTree(-1);
  multisim_tree_copy->Write();
    

  output->Close();
  input->Close();

  return 0;
}

void copy_no_syst(cfg::ConfigurationTable table, TTree * out_tree, TTree * in_tree)
{


  TObjArray* branches = in_tree->GetListOfBranches();
  
  // Input tree
  int run, subrun, event;
  double br[in_tree->GetNbranches()-3];

  //for (int i = 0; i < in_tree->GetNbranches()-3; i++)
  //  in_tree->SetBranchAddress(in_tree->GetListOfBranches()->At(i)->GetName(), br+i);

  // Test
  int temp_i(0);
  for (int i = 0; i < in_tree->GetNbranches(); i++)
    {
      TBranch* branch = static_cast<TBranch*>(branches->At(i));
      const std::string branchName = branch->GetName();
      
      if(branchName != "Run" && branchName != "Subrun" && branchName != "Evt")
	{
	  //std::cout << branchName << std::endl;
	  in_tree->SetBranchAddress(branch->GetName(), br+temp_i);
	  temp_i++;
	}
    }


  in_tree->SetBranchAddress("Run", &run);
  in_tree->SetBranchAddress("Subrun", &subrun);
  in_tree->SetBranchAddress("Evt", &event);

  // These are branches we wish to modify
  double _cut_type, _is_nu, _is_data, _category_topology_ccpi0_xsec1, _interaction_mode;
  in_tree->SetBranchAddress("reco_cut_type", &_cut_type);
  in_tree->SetBranchAddress("reco_is_nu", &_is_nu);
  in_tree->SetBranchAddress("reco_is_data", &_is_data);
  //in_tree->SetBranchAddress("true_category_topology_ccpi0_simple1", &_category_topology_ccpi0_simple1);
  in_tree->SetBranchAddress("true_category_topology_ccpi0_xsec1", &_category_topology_ccpi0_xsec1);
  //in_tree->SetBranchAddress("true_category_topology_ccpi0_complete", &_category_topology_ccpi0_complete);
  in_tree->SetBranchAddress("true_interaction_mode", &_interaction_mode);

  
  // Output tree
  //for (int i = 0; i < in_tree->GetNbranches()-3; i++)
  //  out_tree->Branch(in_tree->GetListOfBranches()->At(i)->GetName(), br+i);

  // Test
  temp_i = 0;
  for (int i = 0; i < in_tree->GetNbranches(); i++)
    {
      TBranch* branch = static_cast<TBranch*>(branches->At(i));
      const std::string branchName = branch->GetName();

      if(branchName != "Run" && branchName != "Subrun" && branchName != "Evt")
        {
          out_tree->Branch(branch->GetName(), br+temp_i);
          temp_i++;
        }
    }

  out_tree->Branch("Run", &run);
  out_tree->Branch("Subrun", &subrun);
  out_tree->Branch("Evt", &event);

  int cut_type, is_nu, is_data, category_topology_ccpi0_xsec1, interaction_mode;
  out_tree->Branch("cut_type", &cut_type, "cut_type/I");
  out_tree->Branch("is_nu", &is_nu, "is_nu/I");
  out_tree->Branch("is_data", &is_data, "is_data/I");
  //out_tree->Branch("category_topology_ccpi0_simple1", &category_topology_ccpi0_simple1, "category_topology_ccpi0_simple1/I");
  out_tree->Branch("category_topology_ccpi0_xsec1", &category_topology_ccpi0_xsec1, "category_topology_ccpi0_xsec1/I");
  //out_tree->Branch("category_topology_ccpi0_complete", &category_topology_ccpi0_complete, "category_topology_ccpi0_complete/I");
  out_tree->Branch("interaction_mode", &interaction_mode, "interaction_mode/I");

  // Copy entries from input tree to output tree
  for(int i(0); i < in_tree->GetEntries(); ++i)
    {
      in_tree->GetEntry(i);
      
      // Any modifications happen here...

      // Make "cut_type" an int
      cut_type = (int)_cut_type;

      // Make "interaction mode" an int
      interaction_mode = (int)_interaction_mode;

      // Replace dummy value for "is_nu" (currently -5)
      if(table.get_bool_field("is_nu") == true)
	is_nu = 1;
      else
	is_nu = 0;

      // Replace dummy value for "is_data" (currently -5)
      if(table.get_bool_field("is_data") == true)
	{
	  is_data = 1;
	  //category_topology_ccpi0_simple1 = (int)10;
	  category_topology_ccpi0_xsec1 = (int)10;
	  //category_topology_ccpi0_complete = (int)10;
	}
      else
	{
	  is_data = 0;
	  //category_topology_ccpi0_simple1 = (int)_category_topology_ccpi0_simple1;
	  category_topology_ccpi0_xsec1 = (int)_category_topology_ccpi0_xsec1;
	  //category_topology_ccpi0_complete = (int)_category_topology_ccpi0_complete;
	}
      
      out_tree->Fill();
    }

}

void copy_with_syst(cfg::ConfigurationTable config, cfg::ConfigurationTable table, TTree * out_tree, TTree * in_tree, TTree* syst_in_tree, std::string syst_type, std::string sel_or_sig)
{

  // Get relevant syst branches from TOML
  int num_syst_branches(0);
  std::vector<std::string> syst_branch_names;
  for(cfg::ConfigurationTable & t : config.get_subtables("sys"))
    {
      if(!strcmp(t.get_string_field("type").c_str(), syst_type.c_str()))
        {
          syst_branch_names.push_back(t.get_string_field("name").c_str());
          num_syst_branches++;
        }
    }

  // Specify morph knobs as these are handled differently
  std::vector<std::string> morph_branch_names = {"GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"};

  // FSI NuSyst knobs are also handled differently
  std::vector<std::string> fsi_nusyst_branch_names = {"FSIReweight_SBNNuSyst_FSI_hNReweight_multisigma_FSIReweight", "FSIReweight_SBNNuSyst_FSI_INCLReweight_multisigma_FSIReweight", "FSIReweight_SBNNuSyst_FSI_G4BCReweight_multisigma_FSIReweight"};

  // Initialize TClonesArray 
  TClonesArray *arrSyst[num_syst_branches];

  // Get the list of branches from input syst tree
  TObjArray *branches = syst_in_tree->GetListOfBranches();

  // Create a map to store the addresses of the vectors
  std::map<std::string, std::vector<Float_t>*> branch_addresses_nsigmas;
  std::map<std::string, std::vector<Float_t>*> branch_addresses_weights;
  std::map<std::string, int> sys_idxs;
  std::vector<TString> relevant_branch_names;

  // Loop through each branch and set its address 
  int sysIdx(0);
  for (int i = 0; i < branches->GetEntries(); ++i)
    {
      TBranch *branch = (TBranch*)branches->At(i);
      TString branch_name = branch->GetName();

      if(std::find(syst_branch_names.begin(), syst_branch_names.end(), branch_name) != syst_branch_names.end())
        {
          sys_idxs[branch_name.Data()] = sysIdx;
	  relevant_branch_names.push_back(branch_name);

	  std::vector<Float_t> *vec2 = nullptr;
          branch_addresses_weights[branch_name.Data()] = vec2;
          syst_in_tree->SetBranchAddress(branch_name, &(branch_addresses_weights[branch_name.Data()]));

          TString temp(branch_name);
          TString dummy("_sigma");
          TString branch_name_nsigma = temp + dummy;

	  std::vector<Float_t> *vec1 = nullptr;
          branch_addresses_nsigmas[branch_name.Data()] = vec1;
          syst_in_tree->SetBranchAddress(branch_name_nsigma, &(branch_addresses_nsigmas[branch_name.Data()]));

	  //std::cout << branch_name_nsigma << " " << branch_name << std::endl;

          arrSyst[sysIdx] = new TClonesArray("TGraph", 1);
          out_tree->Branch(branch_name, &arrSyst[sysIdx], 32000, -1);
          sysIdx++;
        }
    }  

  // Loop through the entries of the tree and access the data 
  Long64_t nEntries = syst_in_tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    syst_in_tree->GetEntry(i);

    // Loop over each systematic parameter
    for (auto const& [name, y_vec_ptr] : branch_addresses_weights)
      {
	std::vector<Float_t>* x_vec_ptr = branch_addresses_nsigmas.at(name);
	std::vector<Float_t> nsigmas;
        for (Float_t val : *x_vec_ptr)
          {
            nsigmas.push_back(val);
          }
	std::vector<Float_t> weights;
        for (Float_t val : *y_vec_ptr)
          {
	    if(table.get_bool_field("is_nu") == false)
	      weights.push_back(-5);
	    else if(!strcmp(syst_type.c_str(), "variation") && sel_or_sig == "sig")
	      weights.push_back(1);
	    else
	      weights.push_back(val);
          }

	// weight = 1 for nsigma = 0
	if(!strcmp(syst_type.c_str(), "multisigma") && std::find(fsi_nusyst_branch_names.begin(), fsi_nusyst_branch_names.end(), name) == fsi_nusyst_branch_names.end())
	  {
	    nsigmas.push_back(0);
	    if(table.get_bool_field("is_nu") == false)
	      weights.push_back(-5);
	    else
	      weights.push_back(1);
	  }

	// Morph knobs
	if(!strcmp(syst_type.c_str(), "multisigma") && std::find(morph_branch_names.begin(), morph_branch_names.end(), name) != morph_branch_names.end())
	  {
	    for (const auto& z : {-1.0, -0.5, 0.0, 0.5, 1.0})
	      {
		// Values we already handled
		if(z == 0.0 || z == 1.0) continue;
		
		// Other values
		nsigmas.push_back(z);
		weights.push_back( 1 + fabs(z)*( y_vec_ptr->front() - 1 ) );

	      }
	  }

        // Create a TGraph for every event
        TGraph *graph = new TGraph(nsigmas.size(), &nsigmas[0], &weights[0]);
        graph->Sort();
        new( (*arrSyst[sys_idxs[name]])[0]) TGraph(nsigmas.size(),graph->GetX(),graph->GetY());
      } // end systematic parameter loop
    
    // This is needed in order to not fill with output TTree with duplicates
    // (only fill specified systematic branches)
    TObjArray *branches = out_tree->GetListOfBranches();
    for (int j = 0; j < branches->GetEntries(); ++j)
      {
        TBranch *branch = (TBranch*)branches->At(j);
	TString branch_name = branch->GetName();
	if(std::find(relevant_branch_names.begin(), relevant_branch_names.end(), branch_name) != relevant_branch_names.end())
	  {
	    branch->Fill();
	  }
      }
  } // end entry loop

} // end function

