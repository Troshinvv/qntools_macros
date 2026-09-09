//
// Created by Misha on 3/7/2023.
//

#include <cmath>
#include <vector>
#include <cassert>

void run8_proton_qa( std::string list, 
                          std::string str_effieciency_file,
			  std::string str_run_id_efficiency_file,
                          std::string centrality_calib_file,
                          std::string calib_in_file="qa.root" ){

  std::cout << "starting execution" << std::endl;
 gInterpreter->GenerateDictionary("ROOT::RVec<ROOT::RVec<float>>", "ROOT/RVec.hxx"); 
  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

	auto f1_2212_m_400 = new TF1("2212_mean_400", "pol1");
  f1_2212_m_400->SetParameter(0, 0.9666099);
  f1_2212_m_400->SetParameter(1, -0.04425819);

  auto f1_2212_s_400 = new TF1("2212_sigma_400", "pol2");
  f1_2212_s_400->SetParameter(0, 0.08997859);
  f1_2212_s_400->SetParameter(1, -0.03365721);
  f1_2212_s_400->SetParameter(2, 0.01650391);

  auto f1_2212_m_700 = new TF1("2212_mean_700", "pol1");
  f1_2212_m_700->SetParameter(0, 0.9562788);
  f1_2212_m_700->SetParameter(1, -0.03858193);

  auto f1_2212_s_700 = new TF1("2212_sigma_700", "pol2");
  f1_2212_s_700->SetParameter(0, 0.05129182);
  f1_2212_s_700->SetParameter(1, -0.0120711);
  f1_2212_s_700->SetParameter(2, 0.01661445);
  auto file_fit = TFile::Open( centrality_calib_file.c_str(), "READ" );
	file_fit->cd();
	auto g1_FitRunIdFactor_1 = file_fit->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_RefMult_gt_mpd_8150_8200");

auto file_run_id_eff = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( str_run_id_efficiency_file.c_str(), "READ" ), [](auto f){f->Close(); } };
  assert(file_run_id_eff);
  THn* efficiency_eta_pT_phi_run_id{nullptr};
        file_run_id_eff->GetObject("hn_efficiency", efficiency_eta_pT_phi_run_id);
  assert(efficiency_eta_pT_phi_run_id);

const auto trWeightFunction = [efficiency_eta_pT_phi_run_id](
     ROOT::VecOps::RVec<float> pt_vec, ROOT::VecOps::RVec<float> eta_vec, ROOT::VecOps::RVec<float> phi_vec, UInt_t run_id
  ){
    std::vector<float> vec_efficiency(pt_vec.size(), 0.f);
                for( int i=0; i<pt_vec.size(); i++ ){
                        auto pT = pt_vec.at(i);
                        auto eta = eta_vec.at(i);
                        auto phi = phi_vec.at(i);
                        auto coord = std::vector<double>{ eta, pT, phi, static_cast<double>(run_id) };
                        auto bin = efficiency_eta_pT_phi_run_id->GetBin( coord.data() );
                        if( bin <= 0 )
                                continue;
                        if( bin > efficiency_eta_pT_phi_run_id->GetNbins() )
                                continue;
                        auto weight = efficiency_eta_pT_phi_run_id->GetBinContent(bin);
                        vec_efficiency[i] = weight;
                }
                return vec_efficiency;
  };

  auto ref_mult_generator =
  []( TGraphErrors* g1_calib ){
    return [g1_calib](unsigned long _mult, UInt_t _runId){ return (_mult * g1_calib->Eval( static_cast<double>(_runId) )); };
  };
  auto m2_function = 
  []
  ( ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >> vec_mom, 
    ROOT::VecOps::RVec<double> vec_beta){
      std::vector<float> vec_m2;
      vec_m2.reserve( vec_beta.size() );
      for( size_t i=0; i<vec_mom.size(); i++ ){
        auto p = vec_mom.at(i).P();
        auto p2 = p*p;
        auto beta = vec_beta.at(i);
        auto beta2 = beta*beta;
        auto gamma2 = 1 - beta2;
        auto m2 = beta > -990. ? p2 / beta2 * gamma2 : -999.0;
        vec_m2.push_back( m2 );
      }
      return vec_m2;
    };
  const auto n_sigma_generator = []( auto f1_mean, auto f1_sigma ){
    return 
    [ f1_mean, f1_sigma ]
    ( ROOT::VecOps::RVec<float> vec_pq, ROOT::VecOps::RVec<float> vec_m2 ){
        auto vec_n_sigma = std::vector<float>{};
        vec_n_sigma.reserve( vec_pq.size() );
        for( size_t i=0; i < vec_pq.size(); ++i ){
          auto m2 = vec_m2.at(i);
          auto pq = vec_pq.at(i);
          auto mean = f1_mean->Eval(pq);
          auto sigma = f1_sigma->Eval(pq);
          auto n_sigma = fabs( m2 - mean ) / sigma;
          vec_n_sigma.push_back( pq > 0 ? n_sigma : 999. );
        }
      return vec_n_sigma;
    };
  };
   
	const auto n_sigma_particle_function = 
  []
  ( ROOT::VecOps::RVec<float> n_sigma_400, 
   ROOT::VecOps::RVec<float> n_sigma_700 ){
      std::vector<float> vec_n_simga{};
      vec_n_simga.reserve( n_sigma_400.size() );
      for( int i=0; i<n_sigma_400.size(); ++i ){ 
	
        vec_n_simga.push_back( std::min( TMath::Abs(n_sigma_400.at(i)), TMath::Abs(n_sigma_700.at(i)) ) ); }
      return vec_n_simga;
  };
  
	const auto rapidity_generator = []( auto particle_m, auto y_cm ){
    return 
    [particle_m, y_cm]( ROOT::VecOps::RVec<float> vec_pz, ROOT::VecOps::RVec<float> vec_pq ){
      std::vector<float> vec_y{};
      vec_y.reserve( vec_pz.size() );
      for( int i=0; i<vec_pz.size(); ++i ){
        auto pz = vec_pz.at(i);
        auto p = vec_pq.at(i);
        auto E = sqrt( p*p + particle_m*particle_m );
        auto y = 0.5 * log( ( E + pz ) / ( E - pz ) ) - y_cm;
        vec_y.push_back( y );
      }
      return vec_y;
    };
  };
  const auto function_fhcal_x = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> vec_param ){
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
  const auto function_fhcal_y = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> vec_param ){
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
//10%
      auto GetCentrBin_8150_8200 = [](Double_t _refMult){ //RunId_corr_factor_h2_RunId_nTracks_8150_8200
	Double_t binMult = _refMult+1;
        if(_refMult < 0) return -1.;//6
        else if ( _refMult < 5 ) return 85.; // 70-100%
        else if ( _refMult < 8 ) return 65.; // 60-70%
        else if ( _refMult < 13 ) return 55.; // 50-60%
        else if ( _refMult < 20 ) return 45.; // 40-50%
        else if ( _refMult < 28 ) return 35.; // 30-40%
        else if ( _refMult < 40 ) return 25.; // 20-30%
        else if ( _refMult < 56) return 15.; // 10-20%
        else if ( _refMult < 114) return 5.; //  0-10%
        else if ( _refMult >=114) return -1.;
        return -1.;
    };
//5%
/*auto GetCentrBin_8150_8200 = [](Double_t _refMult){ //RunId_corr_factor_h2_RunId_nTracks_8150_8200
	Double_t binMult = _refMult+1;
        if( binMult < 0) return -1.;//6
	else if ( binMult < 5 ) return 85.; // 70-100%
        else if ( binMult < 8 ) return 65.; // 60-70%
        else if ( binMult < 13 ) return 55.; // 50-60%
        else if ( binMult < 20 ) return 45.; // 40-50%
        else if ( binMult < 24 ) return 37.5; // 35-40%
        else if ( binMult < 28 ) return 32.5; // 30-35%
        else if ( binMult < 34 ) return 27.5; // 25-30%
        else if ( binMult < 40 ) return 22.5; // 20-25%
        else if ( binMult < 48 ) return 17.5; // 15-20%
        else if ( binMult < 56 ) return 12.5; // 10-15%
        else if ( binMult < 67 ) return 7.5; //  5-10%
        else if ( binMult < 114) return 2.5; //  0-5%
        else if ( binMult >=114) return -1.;
        return -1.;
    };*/
  
  const auto weight_generator = []( auto efficiency_map ){
    return [efficiency_map](std::vector<float> vec_y, ROOT::VecOps::RVec<float> vec_pT){
      if( !efficiency_map ){
          return std::vector<float>(vec_y.size(), 1);
        }
      std::vector<float> vec_weight{};
      vec_weight.reserve(vec_y.size());
      for( int i=0; i<vec_y.size(); ++i ){
        auto pT = vec_pT.at(i);
        auto y = vec_y.at(i);
        auto y_bin = efficiency_map->GetXaxis()->FindBin( y );
        auto pT_bin = efficiency_map->GetYaxis()->FindBin( pT );
        auto efficiency = efficiency_map->GetBinContent( y_bin, pT_bin );
        auto weight = efficiency > 1e-2 ? 1.0 / efficiency : 0.0;
        vec_weight.push_back( weight );
      }
      return vec_weight;
    };
  };

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2D* efficiency_histo{nullptr};
  TH2D* efficiency_tof400{nullptr};
  TH2D* efficiency_tof700{nullptr};

  effieciency_file->GetObject("efficiency_2212_tof", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No efficiency for both tof was found in file " << str_effieciency_file << "\n";
  effieciency_file->GetObject("efficiency_2212_tof400", efficiency_tof400);
  if( !efficiency_tof400 )
    std::cerr << "Warning: No efficiency for tof-400 was found in file " << str_effieciency_file << "\n";
  effieciency_file->GetObject("efficiency_2212_tof700", efficiency_tof700);
  if( !efficiency_tof700 )
    std::cerr << "Warning: No efficiency for tof-700 was found in file " << str_effieciency_file << "\n";
  
  std::vector<int> physical_runs{6667, 6668, 6669, 6670, 6671, 6672, 6673, 6674, 6675, 6676, 6677, 6678, 6679, 6680, 6681, 6683, 6684, 6685, 666, 6687, 6689, 6690, 6691, 6692, 6694, 6695, 6696, 6698, 6699, 6732, 6733, 6734, 6737, 6738, 6739, 6740, 6745, 6752, 6753, 6760, 6761, 6765, 6766, 6767, 6768, 6769, 6771, 6772, 6773, 6774, 6779, 6780, 6782, 6783, 6785, 6786, 6788, 6794, 6795, 6797, 6799, 6800, 6803, 6815, 6816, 6817, 6818, 6819, 6820, 6821, 6822, 6879, 6882, 6883, 6884, 6886, 6887, 6889, 6891, 6900, 6901, 6902, 6903, 6904, 6905, 6906, 6907, 6908, 6909, 6910, 6911, 6915, 6916, 6918, 6919, 6920, 6921, 6923, 6924, 6926, 6927, 6928, 6929, 6930, 6931, 6932, 6933, 6934, 6935, 6936, 6937, 6939, 6940, 6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6990, 6991, 6992, 6993, 6994, 6995, 6997, 6998, 6999, 7000, 7002, 7003, 7004, 7005, 7006, 7008, 7009, 7010, 7011, 7012, 7030, 7031, 7032, 7033, 7034, 7035, 7037, 7038, 7040, 7041, 7042, 7043, 7044, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7075, 7076, 7077, 7078, 7081, 7082, 7083, 7084, 7086, 7087, 7091, 7092, 7093, 7094, 7096, 7097, 7098, 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7842, 7843, 7845, 7846, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7941, 7942, 7944, 7948, 7949, 7950, 7952, 7954, 7955, 7957, 7958, 7960, 7961, 7962, 7963, 7965, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8026, 8027, 8028, 8029, 8030, 8031, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8217, 8219, 8220, 8221, 8228, 8229, 8230, 8231, 8235, 8236, 8238, 8239, 8240, 8242, 8244, 8245, 8246, 8247, 8248, 8250, 8251, 8253, 8254, 8255, 8256, 8257, 8258, 8265, 8266, 8267, 8268, 8270, 8271, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8281, 8284, 8286, 8287, 8288, 8289, 8290, 8292, 8293, 8294, 8295, 8297, 8298, 8299, 8300, 8305, 8306};
  std::vector<int> bad_runs{7313, 7415, 7417, 7435, 7469, 7517, 7519, 7520, 7537, 7575, 7604, 7630, 7657, 7659, 7679, 7681, 7705, 7735, 7843, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7907, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7954, 7955, 8031, 8032, 8033, 8115, 8121, 8167, 8201, 8204, 8205, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8247, 8265, 8266, 8267, 8281, 8289};

  std::vector<int> f1_modules = {
    35, 36, 
    37, 38, 
    39, 40, 
    41, 42, 
    43, 44, 
  };
  std::vector<int> f2_modules = {
     1,  2,  3,  4,  5,
     6,  7,  8,  9,  10,
    11, 12, 13, 14, 15,
    16, 17,     18, 19, 
    20, 21, 22, 23, 24,
    25, 26, 27, 28, 29,
    30, 31, 32, 33, 34,
  };
  std::vector<int> f3_modules = {
    45, 46,
    47, 48,
    49, 50,
    51, 52, 
    53, 54,
  };

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  auto dd=d
	  .Define( "track_multiplicity_gt_corr", ref_mult_generator( g1_FitRunIdFactor_1 ), {"track_multiplicity_gt", "runId"} )
          .Define("centrality",  [GetCentrBin_8150_8200] (Double_t _refMult) { return GetCentrBin_8150_8200(_refMult);},{"track_multiplicity_gt_corr"})
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
          .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
//	  .Define( "pz", [&](ROOT::VecOps::RVec<float> mom,ROOT::VecOps::RVec<float> momt){std::vector<float> pz;for(int i=0;i<mom.size();i++){float pz1 = TMath::Sqrt(mom.at(i)*mom.at(i)-momt.at(i)*momt.at(i));pz.push_back(pz1);}return pz;},{"trP","trPt"})
          // .Define( "trM2Tof700", m2_function, { "trMom", "trBetaTof700" } )
          // .Define( "trM2Tof400", m2_function, { "trMom", "trBetaTof400" } )
  //        .Define( "trNsigmaProton400", n_sigma_generator(f1_2212_m_400, f1_2212_s_400), { "trPq", "trM2Tof400_corr" } )
 //         .Define( "trNsigmaProton700", n_sigma_generator(f1_2212_m_700, f1_2212_s_700), { "trPq", "trM2Tof700_corr"  } )
          .Define( "trNsigmaProton", n_sigma_particle_function, {"trNsigma_2212_400", "trNsigma_2212_700"} )
	  .Define( "trNsigmaDeit", n_sigma_particle_function, {"trNsigma_1000010020_400", "trNsigma_1000010020_700"} )
          .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"trPz", "trPq"} )
	  .Define( "trWeight", trWeightFunction, {"trPt","trEta","trPhi", "runId"} )
          .Define( "trProtonEfficiency", weight_generator(efficiency_histo), {"trProtonY", "trPt"} )
          .Define( "trProtonEfficiencyTof400", weight_generator(efficiency_tof400), {"trProtonY", "trPt"} )
          .Define( "trProtonEfficiencyTof700", weight_generator(efficiency_tof700), {"trProtonY", "trPt"} )
          .Define( "trProtonWeight", "std::vector<double> weights{}; for( auto i=size_t{0}; i<trWeight.size(); ++i ){ weights.push_back( trWeight[i]*trProtonEfficiency[i] ); } return weights;" )
	  .Define( "Weight", []( std::vector<float> gen_w,ROOT::VecOps::RVec<float> chi2,ROOT::VecOps::RVec<int> nhits, ROOT::VecOps::RVec<float> dcaR){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(gen_w.size());
                        for( int i=0; i<gen_w.size(); ++i ){
                        if(gen_w.at(i)==0){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(chi2.at(i)>=5){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(dcaR.at(i)>=5){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(nhits.at(i)<6){
                                vec_weight.push_back(0);
                                continue;
                        }
                        vec_weight.push_back(gen_w.at(i));
                        }
                  return vec_weight;
                }, {"trWeight","stsTrackChi2Ndf","stsTrackNhits","trDcaR"} )
	  .Define( "Weightone", []( std::vector<float> gen_w,ROOT::VecOps::RVec<float> chi2,ROOT::VecOps::RVec<int> nhits, ROOT::VecOps::RVec<float> dcaR){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(gen_w.size());
                        for( int i=0; i<gen_w.size(); ++i ){
                        if(chi2.at(i)>=5){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(dcaR.at(i)>=5){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(nhits.at(i)<6){
                                vec_weight.push_back(0);
                                continue;
                        }
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"trWeight","stsTrackChi2Ndf","stsTrackNhits","trDcaR"} )
	 .Define( "isProt400", []( std::vector<float> gen_w,ROOT::VecOps::RVec<float> nsigma_400,std::vector<float> nsigma_d){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(gen_w.size());
                        for( int i=0; i<gen_w.size(); ++i ){
                        if(gen_w.at(i)==0){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(TMath::Abs(nsigma_400.at(i))>3){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(nsigma_d.at(i)<=3){
				vec_weight.push_back(0);
                                continue;
                        }
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"Weightone","trNsigma_2212_400","trNsigmaDeit"} )
	.Define( "isProt700", []( std::vector<float> gen_w,ROOT::VecOps::RVec<float> nsigma_700,std::vector<float> nsigma_d){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(gen_w.size());
                        for( int i=0; i<gen_w.size(); ++i ){
                        if(gen_w.at(i)==0){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(TMath::Abs(nsigma_700.at(i))>3){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(nsigma_d.at(i)<=3){
                                vec_weight.push_back(0);
                                continue;
                        }
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"Weightone","trNsigma_2212_700","trNsigmaDeit"} )
	.Define( "isProt", []( std::vector<float> gen_w,std::vector<float> nsigma,std::vector<float> nsigma_d){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(gen_w.size());
                        for( int i=0; i<gen_w.size(); ++i ){
                        if(gen_w.at(i)==0){
                                vec_weight.push_back(0);
                                continue;
                        }
                        if(nsigma.at(i)>3){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(nsigma_d.at(i)<=3){
                                vec_weight.push_back(0);
                                continue;
                        }
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"Weightone","trNsigmaProton","trNsigmaDeit"} )
	 .Define( "isPositive", "ROOT::VecOps::RVec<float> x; for(auto& pos:trCharge){ if(pos>=0){x.push_back(1);}else{x.push_back(0);}}; return x;")
       	 .Define( "isNegative", "ROOT::VecOps::RVec<float> x; for(auto& pos:trCharge){ if(pos<0){x.push_back(1);}else{x.push_back(0);}}; return x;")
	 .Define( "Tpos", []( ROOT::VecOps::RVec<short> vec_charge,std::vector<float> vec_x,std::vector<float> vec_y,ROOT::VecOps::RVec<float> eta_vec){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(vec_charge.size());
                        for( int i=0; i<vec_charge.size(); ++i ){
                        if(vec_charge.at(i)<0){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(eta_vec.at(i)>3){
                                vec_weight.push_back(0);
                                continue;
                        }
			
                      /*  if(vec_y.at(i)>-100 && vec_y.at(i)<100){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(vec_x.at(i)>-100 && vec_x.at(i)<250){
                                vec_weight.push_back(0);
                                continue;
                        }*/
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"trCharge","trFhcalX","trFhcalY","trEta"} )
	.Define( "Tneg", []( ROOT::VecOps::RVec<short> vec_charge,std::vector<float> vec_x,std::vector<float> vec_y,ROOT::VecOps::RVec<float> eta_vec){
                        std::vector<float> vec_weight{};
                  vec_weight.reserve(vec_charge.size());
                        for( int i=0; i<vec_charge.size(); ++i ){
                        if(vec_charge.at(i)>=0){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(eta_vec.at(i)>3){
                                vec_weight.push_back(0);
                                continue;
                        }
                      /*  if(vec_y.at(i)>-100 && vec_y.at(i)<100){
                                vec_weight.push_back(0);
                                continue;
                        }
			if(vec_x.at(i)>-100 && vec_x.at(i)<250){
                                vec_weight.push_back(0);
                                continue;
                        }*/
                        vec_weight.push_back(1);
                        }
                  return vec_weight;
                }, {"trCharge","trFhcalX","trFhcalY","trEta"} )
	  .Define( "trTof700ResX", [](ROOT::VecOps::RVec<float> vec_P,ROOT::VecOps::RVec<float> vec_x,ROOT::VecOps::RVec<int> hit_num){
                  std::vector<float> vec_T_Tp{};
                  vec_T_Tp.reserve(vec_P.size());
                        for( int i=0; i<vec_P.size(); ++i ){
                        if(hit_num.at(i)==-1){
                                vec_T_Tp.push_back(-999);
                                continue;
                        }
                        vec_T_Tp.push_back((vec_x.at(hit_num.at(i))));
                        }
                  return vec_T_Tp;
                }, {"trP", "tof700hitResX","trTof700hit"} )
          .Define( "trTof400ResX", [](ROOT::VecOps::RVec<float> vec_P,ROOT::VecOps::RVec<float> vec_x,ROOT::VecOps::RVec<int> hit_num){
                  std::vector<float> vec_T_Tp{};
                  vec_T_Tp.reserve(vec_P.size());
                        for( int i=0; i<vec_P.size(); ++i ){
                        if(hit_num.at(i)==-1){
                                vec_T_Tp.push_back(-999);
                                continue;
                        }
                        vec_T_Tp.push_back((vec_x.at(hit_num.at(i))));
                        }
                  return vec_T_Tp;
                }, {"trP","tof400hitResX","trTof400hit"} )
          .Define( "trTof700ResY", [](ROOT::VecOps::RVec<float> vec_P,ROOT::VecOps::RVec<float> vec_y,ROOT::VecOps::RVec<int> hit_num){
                  std::vector<float> vec_T_Tp{};
                  vec_T_Tp.reserve(vec_P.size());
                        for( int i=0; i<vec_P.size(); ++i ){
                        if(hit_num.at(i)==-1){
                                vec_T_Tp.push_back(-999);
                                continue;
                        }
                        vec_T_Tp.push_back((vec_y.at(hit_num.at(i))));
                        }
                  return vec_T_Tp;
                }, {"trP", "tof700hitResY","trTof700hit"} )
          .Define( "trTof400ResY", [](ROOT::VecOps::RVec<float> vec_P,ROOT::VecOps::RVec<float> vec_y,ROOT::VecOps::RVec<int> hit_num){
                  std::vector<float> vec_T_Tp{};
                  vec_T_Tp.reserve(vec_P.size());
                        for( int i=0; i<vec_P.size(); ++i ){
                        if(hit_num.at(i)==-1){
                                vec_T_Tp.push_back(-999);
                                continue;
                        }
                        vec_T_Tp.push_back((vec_y.at(hit_num.at(i))));
                        }
                  return vec_T_Tp;
                }, {"trP","tof400hitResY","trTof400hit"} )
          .Alias("trStsNhits", "stsTrackNhits")
          .Alias("trStsChi2", "stsTrackChi2Ndf")
          .Filter([&physical_runs, &bad_runs]( UInt_t run_id ){ 
            if( std::find( physical_runs.begin(), physical_runs.end(), run_id) == physical_runs.end() )
              return false;
            if( std::find( bad_runs.begin(), bad_runs.end(), run_id) != bad_runs.end() )
              return false;
            return true;
          }, {"runId"} )
          .Filter("runId < 8312")
          .Filter( []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }, {"triggerMapAR"} )
          .Filter("vtxNtracks >= 2")
          .Filter("fabs(vtxZcorr)<1")
	  .Filter("fabs(vtxRcorr)<1")
          .Filter("noPileup == 1")
//	  .Filter("centrality<30 && centrality>10")
  ; // at least one filter is mandatory!!!

  std::vector<ROOT::RDF::RResultPtr<::TH2D>> hist2d;
  std::vector<ROOT::RDF::RResultPtr<::TH3D>> hist3d;
  std::vector<ROOT::RDF::RResultPtr<::TH1D>> hist1d;
  std::vector<ROOT::RDF::RResultPtr<::THnD>> hist4d;
std::vector<double> gapy_vec={0,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0};
std::vector<double> gappt_vec={0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,10.0};
std::vector<double> gaprunid_vec;
for(int i=7000;i<8320;i++){
	gaprunid_vec.push_back(i-0.5);
}

double gapphi[9];
std::vector<double> gapphi_vec;
for(int i=0;i<9;i++){
gapphi[i]=-TMath::Pi()+TMath::Pi()*2/8*i;
gapphi_vec.push_back(gapphi[i]);
}
hist1d.push_back( dd.Histo1D( {"h1_z",";;",100,-1,1},"vtxZcorr"));
hist1d.push_back( dd.Histo1D( {"h1_runid",";;",1500,6900,8400},"runId"));
/*hist2d.push_back( dd.Histo2D( {"h2_pt_eta_positive",";#eta;pT",500,-5,20,60,0,3},"trEta","trPt","Tpos"));
hist2d.push_back( dd.Histo2D( {"h2_pt_eta_negative",";#eta;pT",500,-5,20,60,0,3},"trEta","trPt","Tneg"));
hist2d.push_back( dd.Histo2D( {"h2_cent_eta_positive",";#eta;centrality",500,-5,20,20,0,100},"trEta","centrality","Tpos"));
hist2d.push_back( dd.Histo2D( {"h2_cent_eta_negative",";#eta;centrality",500,-5,20,20,0,100},"trEta","centrality","Tneg"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalX_eta_positive",";#eta;trFhcalX",500,-5,20,2000,-2000,2000},"trEta","trFhcalX","Tpos"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalY_eta_positive",";#eta;trFhcalY",500,-5,20,1000,-1000,1000},"trEta","trFhcalY","Tpos"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalX_eta_negative",";#eta;trFhcalX",500,-5,20,2000,-2000,2000},"trEta","trFhcalX","Tneg"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalY_eta_negative",";#eta;trFhcalY",500,-5,20,1000,-1000,1000},"trEta","trFhcalY","Tneg"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalX_Y_positive",";trFhcalX;trFhcalY",2000,-2000,2000,1000,-1000,1000},"trFhcalX","trFhcalY","Tpos"));
hist2d.push_back( dd.Histo2D( {"h2_trFhcalX_Y_negative",";trFhcslX;trFhcalY",2000,-2000,2000,1000,-1000,1000},"trFhcalX","trFhcalY","Tneg"));*/
/*hist1d.push_back( dd.Profile1D( {"p_y_runid_weight",";runId;y_{cm}",1500,6900,8400}, "runId","trProtonY","Weight"));
hist1d.push_back( dd.Profile1D( {"p_pt_runid_weight",";runId;p_{T}, GeV/c",1500,6900,8400}, "runId","trPt","Weight"));
hist1d.push_back( dd.Profile1D( {"p_phi_runid_weight",";runId;#phi, rad",1500,6900,8400}, "runId","trPhi","Weight"));
hist1d.push_back( dd.Profile1D( {"p_y_runid",";runId;y_{cm}",1500,6900,8400}, "runId","trProtonY","Weightone"));
hist1d.push_back( dd.Profile1D( {"p_pt_runid",";runId;p_{T}, GeV/c",1500,6900,8400}, "runId","trPt","Weightone"));
hist1d.push_back( dd.Profile1D( {"p_phi_runid",";runId;#phi, rad",1500,6900,8400}, "runId","trPhi","Weightone"));*/
//hist3d.push_back(dd.Histo3D( {"h3_runid_pt_eta",";#eta;pT;runId",10,0,3,10,0,3,1500,6900,8400},"trEta","trPt","runId","Weightone"));
//hist4d.push_back( dd.HistoND( { "hN_runid_pt_eta_phi",";runId; p_{T} (GeV/c);#eta;phi, rad",4,{1319,8,14,8},{gaprunid_vec,gappt_vec,gapy_vec,gapphi_vec}},{ "runId", "trPt","trEta","trPhi", "Weightone"} ) );
hist2d.push_back( dd.Histo2D( {"h2_nsigma_tof_400_pq",";p/q (GeV/c);#frac{m^{2}-<m^{2}>}{#sigma}",100,0,10,100,-5,5},"trPq","trNsigma_2212_400","Weightone"));
hist2d.push_back( dd.Histo2D( {"h2_nsigma_tof_700_pq",";p/q (GeV/c);#frac{m^{2}-<m^{2}>}{#sigma}",100,0,10,100,-5,5},"trPq","trNsigma_2212_700","Weightone"));
hist2d.push_back( dd.Histo2D( {"h2_m2_tof_400_pq",";p/q (GeV/c);m^{2} (GeV^{2}/c^{4})",100,0,10,1100,-1,10},"trPq","trM2Tof400_corr","isProt400"));
hist2d.push_back( dd.Histo2D( {"h2_m2_tof_700_pq",";p/q (GeV/c);m^{2} (GeV^{2}/c^{4})",100,0,10,1100,-1,10},"trPq","trM2Tof700_corr","isProt700"));
hist2d.push_back( dd.Histo2D( {"h2_400_prot_acceptance",";y_{cm};p_{T} (GeV/c)",400,-1,3,300,0,3},"trProtonY","trPt","isProt400"));
hist2d.push_back( dd.Histo2D( {"h2_700_prot_acceptance",";y_{cm};p_{T} (GeV/c)",400,-1,3,300,0,3},"trProtonY","trPt","isProt700"));
hist2d.push_back( dd.Histo2D( {"h2_prot_acceptance",";y_{cm};p_{T} (GeV/c)",400,-1,3,300,0,3},"trProtonY","trPt","isProt"));
hist2d.push_back( dd.Histo2D( {"h2_prot_y_phi",";#phi (rad);y_{cm}",800,-4,4,400,-1,3},"trPhi","trProtonY","isProt"));
hist2d.push_back( dd.Histo2D( {"h2_resx_prot_400_pq",";pq, GeV/c;ResX, cm",100,0,10,800,-10,10}, "trP","trTof400ResX","isProt400"));
hist2d.push_back( dd.Histo2D( {"h2_resx_prot_700_pq",";pq, GeV/c;ResX, cm",100,0,10,800,-10,10}, "trP","trTof700ResX","isProt700"));
hist2d.push_back( dd.Histo2D( {"h2_resy_prot_400_pq",";pq, GeV/c;ResY, cm",100,0,10,800,-10,10}, "trP","trTof400ResY","isProt400"));
hist2d.push_back( dd.Histo2D( {"h2_resy_prot_700_pq",";pq, GeV/c;ResY, cm",100,0,10,800,-10,10}, "trP","trTof700ResY","isProt700"));
auto file_out = TFile::Open("out_qa.root", "recreate");
        for( auto& p2 : hist2d )
                p2->Write();
	for( auto& p3 : hist3d )
                p3->Write();
	for( auto& p1 : hist1d )
                p1->Write();
	for( auto& p4 : hist4d )
                p4->Write();
        file_out->Close();

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}
