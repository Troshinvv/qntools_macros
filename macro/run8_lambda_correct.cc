//
// Created by Misha on 3/7/2023.
//

void run8_lambda_correct(std::string list, std::string str_efficiency_file){
  std::vector<int> f1_modules = {11, 12, 13, 16, 17, 20, 21, 22};
  std::vector<int> f2_modules = {5, 6, 7, 8, 9, 10, 14, 15, 18, 19, 23, 24, 25, 26, 27, 28};
  std::vector<int> f3_modules = {0, 1, 2, 3, 4, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
  const float FHCAL_Z = 980;
const float PROTON_M = 0.938;
  const float Y_CM = 1.15141;
  const float LAMBDA_M = 1.11568;

  auto centrality_function =
  []
  (ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >> vec_mom){
      float centrality;
      std::vector<float> centrality_percentage{ 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100 };
      std::vector<int> multiplicity_edges{ 248, 150, 122, 100, 82, 67, 54, 43, 34, 21, 13, 8, 4, 3, 1 };
      auto multiplicity = vec_mom.size();
      if( multiplicity > multiplicity_edges[0] )
        return -1.0f;
      int idx = 0;
      float bin_edge = multiplicity_edges[idx];
      while( multiplicity < bin_edge &&
        idx < multiplicity_edges.size()-1 ){
        idx++;
        bin_edge = multiplicity_edges[idx];
      }
      centrality = (centrality_percentage[idx-1] + centrality_percentage[idx])/2.0f;
      return centrality;
    };
  auto function_fhcal_x =
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<std::vector<float>> vec_param ){
      std::vector<float> vec_x{};
      vec_x.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto x = par.at(0);
        auto z = par.at(2);
        auto tx = par.at(3);
        auto dz = FHCAL_Z - z;
        auto dx = tx * dz;
        vec_x.push_back( x+dx );
      }
      return vec_x;
    };
  auto function_fhcal_y =
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<vector<float>> vec_param ){
      std::vector<float> vec_y{};
      vec_y.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto y = par.at(1);
        auto z = par.at(2);
        auto ty = par.at(4);
        auto dz = FHCAL_Z - z;
        auto dy = ty * dz;
        vec_y.push_back( y+dy );
      }
      return vec_y;
    };
  TFilePtr efficiency_file{ str_efficiency_file };
  TH2D* efficiency_histo{nullptr};
  efficiency_file->GetObject( "h2_efficiency", efficiency_histo );
  if( !efficiency_histo ){ std::cout << "Efficiency histogram cannot be retrieved from file" << std::endl; }
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
	  .Define("track_multiplicity", "return trMom.size();")
          .Define( "stsNdigits","return stsDigits.size()" )
          .Define("centrality", centrality_function, {"trMom"} )
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
	  .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
          .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
          .Define("scwallModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
	  .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
	  .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
	  .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("candidate_pT", "std::vector<float> pT; for( auto mom : candidate_momenta ){ pT.push_back( mom.Pt() ); } return pT;")
          .Define("candidate_phi", "std::vector<float> phi; for( auto mom : candidate_momenta ){ phi.push_back( mom.Phi() ); } return phi;")
          .Define("candidate_rapidity", "std::vector<float> rapidity; for( auto mom : candidate_momenta ){ rapidity.push_back( mom.Rapidity() - 1.0 ); } return rapidity;")
          .Define("candidate_weight",
                  [efficiency_histo]( ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >> momenta ){
            std::vector<float> weights;
            for( auto mom : momenta ){
              if( !efficiency_histo ){
                weights.push_back(1.f);
              }
              auto y = mom.Rapidity();
              auto pT = mom.Pt();
              auto y_bin = efficiency_histo->GetXaxis()->FindBin( y );
              auto pT_bin = efficiency_histo->GetYaxis()->FindBin( pT );
              auto efficiency = efficiency_histo->GetBinContent( y_bin, pT_bin );
              auto weight = efficiency > 0.0001 ? 1.0/efficiency : 0.0;
              weights.push_back( weight );
            }
            return weights;
          },
                  {"candidate_momenta"})
          .Define("m_err", "std::vector<float> err; for( auto mom : candidate_momentum_errors ){ err.push_back( mom.at(3) ); } return err;")
          .Define("daughter1_cos", "std::vector<float> cosine; for( int i=0; i<daughter_cosines.at(0).size(); ++i ){ cosine.push_back( daughter_cosines.at(0).at(i) ); } return cosine;")
          .Define("daughter2_cos", "std::vector<float> cosine; for( int i=0; i<daughter_cosines.at(1).size(); ++i ){ cosine.push_back( daughter_cosines.at(1).at(i) ); } return cosine;")
          .Define("daughter1_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(0).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(0).at(i) ); } return chi2;")
          .Define("daughter2_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(1).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(1).at(i) ); } return chi2;")
          .Define("candidate_good",
                  "std::vector<int> good_candidate;\n"
                  "for(int i=0; i<daughter_cosines.at(0).size(); ++i){\n"
		  "if( candidate_cos_topo.at(i) < 0.998 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( daughter2_cos.at(i) < 0.997 || daughter2_cos.at(i) > 0.9998 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( daughter1_chi2_prim.at(i) < 400 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( daughter2_chi2_prim.at(i) < 10 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( candidate_L.at(i) < 2.25 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( candidate_LdL.at(i) < 6.25 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( candidate_chi2_topo.at(i) > 50 ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                 "if( candidate_chi2_geo.at(i) > 20  ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"
                                                "if( daughter_dca.at(i) > 0.7  ){"
                                                        "good_candidate.push_back(0);"
                                                        "continue;"
                                                "}"

                                                "good_candidate.push_back( 1 );"
                                        "} return good_candidate;"
          )
	  .Alias("trStsNhits", "stsTrackNhits")
          .Alias("trStsChi2", "stsTrackChi2Ndf")
	  .Filter([]( ROOT::VecOps::RVec<float> bc1_int, ROOT::VecOps::RVec<float> fd_int ){
            float trigger = fd_int[0]-bc1_int[0]*0.602;
            return -25900 < trigger && trigger < -6030;
          }, {"bc1sIntegral", "fdIntegral"} )
	  .Filter( []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }, {"triggerMapAR"} )
          .Filter([]( unsigned long sts_digits, unsigned long n_tracks ){
            double sts_min = sts_digits-n_tracks*(4.81632+0.0332792*n_tracks-9.62078e-05*n_tracks*n_tracks);
            double sts_max = sts_digits-n_tracks*(19.4203-0.0518774*n_tracks+4.56033e-05*n_tracks*n_tracks);
            return -74.0087 < sts_min && sts_max < 188.248;
          }, {"stsNdigits", "track_multiplicity"} )
          .Filter("vtxChi2/vtxNdf > 0.1")
  ;


  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("centrality|psiRP"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("candidate_(pT|rapidity|phi|weight|mass|good)"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 1, 10, 40} );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f1.AddCut( "fhcalModId", [&f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f2.AddCut( "fhcalModId", [&f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f3.AddCut( "fhcalModId", [&f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  VectorConfig Tp( "Tp", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tp.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  Tp.AddCut( "trEta", [](double eta){
    return 2.0 < eta && eta < 3.0;
  }, "eta cut" );
  Tp.AddCut( "trPt", [](double pT){
    return pT > 0.2;
  }, "pT cut" );
  Tp.AddCut( "trFhcalX", [](double pos){
    return pos < 10.0 || pos > 120;
    }, "cut on x-pos in fhcal plane" );
  Tp.AddCut( "trFhcalY", [](double pos){
    return pos < -50.0 || pos > 50;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tp);

  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "trCharge", [](double charge){
    return charge < 0.0;
    }, "charge" );
  Tneg.AddCut( "trEta", [](double eta){
    return 1.5 < eta && eta < 4.0;
    }, "eta cut" );
  Tneg.AddCut( "trPt", [](double pT){
    return pT > 0.2;
    }, "pT cut" );
  Tneg.AddCut( "trFhcalX", [](double pos){
    return pos < 10.0 || pos > 120;
    }, "cut on x-pos in fhcal plane" );
  Tneg.AddCut( "trFhcalY", [](double pos){
    return pos < -50.0 || pos > 50;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tneg);

  std::vector<Qn::AxisD> rec_lamda_axes{
          { "candidate_rapidity", 6, -0.2, 1.0 },
          { "candidate_pT", 7, 0.0, 1.4 },
          { "candidate_mass", 10, 1.10, 1.13 },
  };

  VectorConfig lambda_good( "lambda_good", "candidate_phi", "candidate_weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_good.SetHarmonicArray( {1, 2} );
  lambda_good.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, } );
  lambda_good.SetCorrectionAxes( rec_lamda_axes );
  lambda_good.AddCut( "candidate_good", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is good candidate" );
  correction_task.AddVector(lambda_good);


  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}
