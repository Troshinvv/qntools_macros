//
// Created by Misha on 3/7/2023.
//

void k_short_correct(std::string list/*, std::string str_efficiency_file*/,std::string calib_in_file="qa.root"){
  std::vector<int> f1_modules = {34, 35,
                                 36, 37,
                                 38, 39,
                                 40, 41,
                                 42, 43};
  std::vector<int> f2_modules = {0, 1, 2, 3, 4,
                                 5, 6, 7, 8, 9,
                                10,11,12,13,14,
                                15,16,   17,18,
                                19,20,21,22,23,
                                24,25,26,27,28,
                                29,30,31,32,33};
  std::vector<int> f3_modules = {44, 45,
                                 46, 47,
                                 48, 49,
                                 50, 51,
                                 52, 53};
//  TFilePtr efficiency_file{ str_efficiency_file };
//  TH2D* efficiency_histo{nullptr};
//  efficiency_file->GetObject( "h2_efficiency", efficiency_histo );
//  if( !efficiency_histo ){ std::cout << "Efficiency histogram cannot be retrieved from file" << std::endl; }
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("candidate_pT", "std::vector<float> pT; for( auto mom : candidate_momenta ){ pT.push_back( mom.Pt() ); } return pT;")
          .Define("candidate_phi", "std::vector<float> phi; for( auto mom : candidate_momenta ){ phi.push_back( mom.Phi() ); } return phi;")
          .Define("candidate_rapidity", "std::vector<float> rapidity; for( auto mom : candidate_momenta ){ rapidity.push_back( mom.Rapidity() - 1.15141 ); } return rapidity;")
          .Define("daughter1_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(0).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(0).at(i) ); } return chi2;")
          .Define("daughter2_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(1).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(1).at(i) ); } return chi2;")
          .Define("candidate_good",
                  "std::vector<int> good_candidate;\n"
                  "for(int i=0; i<daughter_cosines.at(0).size(); ++i){\n"
                  "if( daughter1_chi2_prim.at(i) < 80 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( daughter2_chi2_prim.at(i) < 80 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_L.at(i) < 0.25 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_LdL.at(i) < 6 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_cos_topo.at(i) < 0.99){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_chi2_topo.at(i) > 50 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_chi2_geo.at(i) > 30  ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( daughter_dca.at(i) > 1.5  ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  " good_candidate.push_back( 1 );\n"
                  "}\n"
                  "return good_candidate;\n"
          )
          .Filter("vtxChi2>0.0001"); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("centrality"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("candidate_(pT|rapidity|phi|mass|good)")
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 3, 10, 30} );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f1.AddCut( "fhcalModId", [f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f2.AddCut( "fhcalModId", [f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f3.AddCut( "fhcalModId", [f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  std::vector<Qn::AxisD> rec_k_short_axes{
          { "candidate_rapidity", 3, 0, 1.5 },
          { "candidate_pT", 3, 0.0, 1.5 },
          { "candidate_mass", 20, 0.4, 0.6 },
  };

  VectorConfig k_short_good( "k_short_good", "candidate_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  k_short_good.SetHarmonicArray( {1, 2} );
  k_short_good.SetCorrections( {CORRECTION::PLAIN, CORRECTION::TWIST_RESCALING} );
  k_short_good.SetCorrectionAxes( rec_k_short_axes );
  k_short_good.AddCut( "candidate_good", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is good candidate" );
  k_short_good.AddHisto1D({"candidate_phi",70,-3.5,3.5});
  correction_task.AddVector(k_short_good);


  correction_task.Run();
}
