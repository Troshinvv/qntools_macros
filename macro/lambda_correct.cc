//
// Created by Misha on 3/7/2023.
//
#include <cmath>
#include <math.h>
#include <random>
#include <vector>
void lambda_correct(std::string list, std::string str_efficiency_file,std::string calib_in_file="qa.root"){
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
const float FHCAL_Z = 980;

/*auto function_fhcal_x =
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
*/
std::random_device rd{};
std::mt19937 random_engine( rd() );
std::uniform_real_distribution<double> phi_distribution_{ -M_PI, M_PI };
    const auto function_random_shuffle_phi = 
  [&random_engine, &phi_distribution_]( ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> vec_mom ) mutable {
    ROOT::VecOps::RVec< std::array<double, 2> > vec_shuffled(vec_mom.size());
    auto old_res_x = double{};
    auto old_res_y = double{};
    auto shuffled_rex_x = double{};
    auto shuffled_rex_y = double{};

    std::for_each( vec_shuffled.begin(), vec_shuffled.end(), [ &, i=0 ]( auto& new_mom ) mutable {
      new_mom = std::array<double, 2>{vec_mom.at(i).Px(), vec_mom.at(i).Py()};
      auto dphi = phi_distribution_(random_engine);
      auto px = new_mom[0];
      auto py = new_mom[1];
      old_res_x+=px;
      old_res_y+=px;
      auto new_px = px*cos(dphi) - py*sin(dphi);
      auto new_py = px*sin(dphi) + py*cos(dphi);
      shuffled_rex_x+=new_px;
      shuffled_rex_y+=new_py;
      new_mom[0] = new_px;
      new_mom[1] = new_py;
      ++i;
    } );
    auto mov_x = (shuffled_rex_x - old_res_x) / vec_shuffled.size();
    auto mov_y = (shuffled_rex_y - old_res_y) / vec_shuffled.size();
    std::vector<double> vec_shuffled_phi{};
    std::for_each( vec_shuffled.begin(), vec_shuffled.end(), [ & ]( auto& new_mom ) mutable {
      new_mom[0] = new_mom[0] - mov_x;
      new_mom[1] = new_mom[1] - mov_y;
      vec_shuffled_phi.push_back( atan2( new_mom[1], new_mom[0] ) );
    });
    return vec_shuffled_phi;
  };
std::unique_ptr<TFile> efficiency_file{TFile::Open( str_efficiency_file.c_str(), "READ" )};
  TH2D* efficiency_histo{nullptr};
  efficiency_file->GetObject( "h2_pT_y_signal", efficiency_histo );
  if( !efficiency_histo ){ std::cout << "Efficiency histogram cannot be retrieved from file" << std::endl; }
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
std::vector<int> physical_runs{7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7941, 7942, 7944, 7948, 7949, 7950, 7952, 7954,7955, 7957, 7958, 7960, 7961, 7962, 7963, 7965, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8026, 8027, 8028, 8029, 8030, 8031, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8217, 8219, 8220, 8221, 8228, 8229, 8230, 8231, 8235, 8236, 8238, 8239, 8240, 8242, 8244, 8245, 8246, 8247, 8248, 8250, 8251, 8253, 8254, 8255, 8256, 8257, 8258, 8265, 8266, 8267, 8268, 8270, 8271, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8281, 8284, 8286, 8287, 8288, 8289, 8290, 8292, 8293, 8294, 8295, 8297, 8298, 8299, 8300, 8305, 8306};
  auto dd=d/*.Filter([&physical_runs]( UInt_t run_id ){
            if( std::find( physical_runs.begin(), physical_runs.end(), run_id) == physical_runs.end() )
              return false;
            return true;
          }, {"runId"} )*/
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("candidate_pT", "std::vector<float> pT; for( auto mom : candidate_momenta ){ pT.push_back( mom.Pt() ); } return pT;")
          .Define("candidate_phi", "std::vector<float> phi; for( auto mom : candidate_momenta ){ phi.push_back( mom.Phi() ); } return phi;")
          .Define("candidate_rapidity", "std::vector<float> rapidity; for( auto mom : candidate_momenta ){ rapidity.push_back( mom.Rapidity() - 1.15141 ); } return rapidity;")
          .Define("daughter1_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(0).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(0).at(i) ); } return chi2;")
          .Define("daughter2_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(1).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(1).at(i) ); } return chi2;")
          .Define( "candidate_Weight", [efficiency_histo](std::vector<float> vec_y, std::vector<float> vec_pT){
                  if( !efficiency_histo ){
                      return std::vector<float>(vec_y.size(), 1);
                    }
                  std::vector<float> vec_weight{};
                  vec_weight.reserve(vec_y.size());
                  for( int i=0; i<vec_y.size(); ++i ){
                    auto pT = vec_pT.at(i);
                    auto y = vec_y.at(i);
                    auto y_bin = efficiency_histo->GetXaxis()->FindBin( y );
                    auto pT_bin = efficiency_histo->GetYaxis()->FindBin( pT );
                    auto efficiency = efficiency_histo->GetBinContent( y_bin, pT_bin );
                    auto weight = efficiency > 1e-2 ? 1.0 / efficiency : 0.0;
                    vec_weight.push_back( weight );
                  }
                  return vec_weight;
          }, {"candidate_rapidity", "candidate_pT"} )
          .Define("candidate_good",
                  "std::vector<int> good_candidate;\n"
                  "for(int i=0; i<daughter_cosines.at(0).size(); ++i){\n"
                  "if( daughter1_chi2_prim.at(i) < 200 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( daughter2_chi2_prim.at(i) < 10 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_L.at(i) < 1.5 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_LdL.at(i) < 6 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_cos_topo.at(i) < 0.999){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_chi2_topo.at(i) > 40 ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( candidate_chi2_geo.at(i) > 30  ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  "if( daughter_dca.at(i) > 0.7  ){\n"
                          "good_candidate.push_back(0);\n"
                          "continue;\n"
                  "}\n"
                  " good_candidate.push_back( 1 );\n"
                  "}\n"
                  "return good_candidate;\n"
          )
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
    //      .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
    //      .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
          .Define("candidate_ShuffledPhi", function_random_shuffle_phi, {"candidate_momenta"})
          .Filter("vtxChi2>0.0001"); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("centrality"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("candidate_(pT|rapidity|phi|mass|good|Weight)")/*,
					    std::regex("tr(Pt|Eta|Phi|FhcalX|FhcalY|Charge)")*/
                                    });
  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 3, 10, 40} );

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
  /*
  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
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
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tneg.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tneg);

  VectorConfig Tpos( "Tpos", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tpos.SetHarmonicArray( {1, 2} );
  Tpos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  Tpos.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  Tpos.AddCut( "trEta", [](double eta){
    return 2.0 < eta && eta < 3.0;
  }, "eta cut" );
  Tpos.AddCut( "trPt", [](double pT){
    return pT > 0.2;
  }, "pT cut" );
  Tpos.AddCut( "trFhcalX", [](double pos){
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tpos.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tpos);
*/
//{ "candidate_rapidity", {-0.2,0,0.3,0.5,0.8,1.0} },
// { "candidate_rapidity", {-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0} },
  std::vector<Qn::AxisD> rec_lamda_axes{
          { "candidate_rapidity", {-0.2,0,0.3,0.5,0.8,1.0} },
          { "candidate_pT", {0,0.2,0.4,0.6,0.8,1.0,1.2,2} },
          { "candidate_mass", 20, 1.09, 1.15 },
  };
  VectorConfig lambda_good( "lambda_good", "candidate_phi", "candidate_Weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_good.SetHarmonicArray( {1, 2} );
  lambda_good.SetCorrections( {CORRECTION::PLAIN,CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING} );
  lambda_good.SetCorrectionAxes( rec_lamda_axes );
  lambda_good.AddCut( "candidate_good", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is good candidate" );
  lambda_good.AddHisto1D({"candidate_phi",70,-3.5,3.5});
  correction_task.AddVector(lambda_good);


/*VectorConfig nonflow_lambda( "lambda_nonflow", "candidate_ShuffledPhi", "candidate_Weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  nonflow_lambda.SetHarmonicArray( {1, 2} );
  nonflow_lambda.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  nonflow_lambda.SetCorrectionAxes( rec_lamda_axes );
  nonflow_lambda.AddCut( "candidate_good", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is good candidate" );
  nonflow_lambda.AddHisto1D({"candidate_phi",70,-3.5,3.5});
  correction_task.AddVector(nonflow_lambda);*/

  correction_task.Run();
}
