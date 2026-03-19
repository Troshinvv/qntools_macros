//
// Created by Misha on 3/7/2023.
//

#include <cmath>
#include <vector>

void run8_proton_correct_clean_runid( std::string list, 
                          std::string str_effieciency_file,
                          std::string centrality_calib_file,
                          std::string calib_in_file="qa.root" ){

  std::cout << "starting execution" << std::endl;
  
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
	auto g1_FitVtxX = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_x");
	auto g1_FitVtxY = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_y");
	auto g1_FitVtxZ = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_z");

	auto g1_FitRunIdFactor_1 = file_fit->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_nTracks_8120_8170");
	auto g1_FitRunIdFactor_2 = file_fit->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_nTracks_7400_7450");

  auto vtx_correction_generator = 
  []( TGraphErrors* g1_calib ){
    return [g1_calib](double _vtx, UInt_t _runId){return _vtx - g1_calib->Eval( static_cast<double>(_runId) ); };
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
    ( std::vector<float> vec_pq, ROOT::VecOps::RVec<float> vec_m2 ){
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
  ( std::vector<float> n_sigma_400, 
    std::vector<float> n_sigma_700 ){
      std::vector<int> vec_n_simga{};
      vec_n_simga.reserve( n_sigma_400.size() );
      for( int i=0; i<n_sigma_400.size(); ++i ){ 
        vec_n_simga.push_back( std::min( n_sigma_400.at(i), n_sigma_700.at(i) ) ); }
      return vec_n_simga;
  };
  
	const auto rapidity_generator = []( auto particle_m, auto y_cm ){
    return 
    [particle_m, y_cm]( std::vector<float> vec_pz, std::vector<float> vec_pq ){
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
  const auto function_fhcal_y = 
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

  const auto centrality_function = 
  []
  (double multiplicity){
      float centrality;
//Glauber
      std::vector<float> centrality_percentage{ 0, 10, 20, 30, 40, 50, 60, 70, 100 };
      std::vector<int> multiplicity_edges{ 236, 137, 99, 71, 49, 33, 22, 12, 0 };
//      std::vector<float> centrality_percentage{ 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 100 };
//      std::vector<int> multiplicity_edges{ 236, 183, 167, 155, 145, 136, 127, 119, 111, 104, 97, 71, 49, 33, 22, 12, 0  };
//	std::vector<int> multiplicity_edges{236, 180, 170, 162, 155, 149, 143, 138, 132, 127, 105,85, 68, 54, 42, 31, 1};
// - gamma fit 25.04 8120-8170
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
  
  const auto dca_function = [](std::vector<float> vec_x, std::vector<float> vec_y){
    std::vector<float> vec_r{};
    vec_r.reserve(vec_x.size());
    for (int i=0; i<vec_x.size(); ++i) {
      auto x = vec_x.at(i);
      auto y = vec_y.at(i);
      auto r = std::sqrt( x*x + y*y );
      vec_r.push_back(r);
    }
    return vec_r;
  };

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
    34, 35, 
    36, 37, 
    38, 39, 
    40, 41, 
    42, 43, 
  };
  std::vector<int> f2_modules = {
     0,  1,  2,  3,  4,
     5,  6,  7,  8,  9,
    10, 11, 12, 13, 14,
    15, 16,     17, 18, 
    19, 20, 21, 22, 23,
    24, 25, 26, 27, 28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_modules = {
    44, 45,
    46, 47,
    48, 49,
    50, 51, 
    52, 53,
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
          .Define( "vtxXcorr", vtx_correction_generator(g1_FitVtxX), {"vtxX","runId"})
          .Define( "vtxYcorr", vtx_correction_generator(g1_FitVtxY), {"vtxY","runId"})
          .Define( "vtxZcorr", vtx_correction_generator(g1_FitVtxY), {"vtxZ","runId"})
          .Define( "vtxRcorr", "return sqrt(vtxXcorr*vtxXcorr + vtxYcorr*vtxYcorr);" )
	  .Define( "vtxRMpd", "return sqrt(vtxXMpd*vtxXMpd + vtxYMpd*vtxYMpd);" )
          .Define("track_multiplicity", "return trMom.size();")
          .Define( "ref_multiplicity", ref_mult_generator( g1_FitRunIdFactor_1 ), {"track_multiplicity","runId"} )
          .Define("stsNdigits","return stsDigits.size()" )
          .Define("centrality", centrality_function, {"ref_multiplicity"} )
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
          .Define( "trDcaX", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(0) - vtxXMpd ); } return vec_par; " )
		      .Define( "trDcaY", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(1) - vtxYMpd ); } return vec_par; " )
          .Define( "trDcaR", dca_function, {"trDcaX", "trDcaY"} )
          .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
          .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
          .Define( "trChi2Ndf", " std::vector<float> vec_par; for( int i=0; i<trChi2.size(); ++i ){ vec_par.push_back( trChi2.at(i)/trNdf.at(i) ); } return vec_par; " )
          .Define( "trPx", " std::vector<float> px; for( auto mom : trMom ){ px.push_back( mom.Px() ); } return px; " )
          .Define( "trPy", " std::vector<float> py; for( auto mom : trMom ){ py.push_back( mom.Py() ); } return py; " )
          .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
          .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P() / trCharge.at(i) ); } return pq;" )
          // .Define( "trM2Tof700", m2_function, { "trMom", "trBetaTof700" } )
          // .Define( "trM2Tof400", m2_function, { "trMom", "trBetaTof400" } )
          .Define( "trNsigmaProton400", n_sigma_generator(f1_2212_m_400, f1_2212_s_400), { "pq", "trM2Tof400" } )
          .Define( "trNsigmaProton700", n_sigma_generator(f1_2212_m_700, f1_2212_s_700), { "pq", "trM2Tof700"  } )
          .Define( "trNsigmaProton", n_sigma_particle_function, {"trNsigmaProton400", "trNsigmaProton700"} )
          .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"pz", "pq"} )
          .Define( "trWeight", weight_generator(efficiency_histo), {"trProtonY", "trPt"} )
          .Define( "trWeightTof400", weight_generator(efficiency_tof400), {"trProtonY", "trPt"} )
          .Define( "trWeightTof700", weight_generator(efficiency_tof700), {"trProtonY", "trPt"} )
          .Alias("trStsNhits", "stsTrackNhits")
          .Alias("trStsChi2", "stsTrackChi2Ndf")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
          .Filter([&physical_runs, &bad_runs]( UInt_t run_id ){ 
            if( std::find( physical_runs.begin(), physical_runs.end(), run_id) == physical_runs.end() )
              return false;
            if( std::find( bad_runs.begin(), bad_runs.end(), run_id) != bad_runs.end() )
              return false;
            return true;
          }, {"runId"} )
          .Filter("runId < 8312")
          .Filter( []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }, {"triggerMapAR"} )
        //  .Filter([]( unsigned long sts_digits, unsigned long n_tracks ){ 
        //    double sts_min = sts_digits-n_tracks*(4.81632+0.0332792*n_tracks-9.62078e-05*n_tracks*n_tracks);
        //    double sts_max = sts_digits-n_tracks*(19.4203-0.0518774*n_tracks+4.56033e-05*n_tracks*n_tracks);
        //    return -74.0087 < sts_min && sts_max < 188.248; 
        //  }, {"stsNdigits", "track_multiplicity"} )
          .Filter("vtxNtracks > 2")
          .Filter("fabs(vtxZMpd)<1")
	  .Filter("fabs(vtxRMpd)<1")
        //  .Filter("fabs(vtxRcorr)<1")
        //  .Filter("fabs(vtxZcorr)<0.1")
          .Filter("noPileup == 1")
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("centrality|runId"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Px|Py|Eta|Phi|NsigmaProton|NsigmaProton400|NsigmaProton700|Charge|ProtonY|DcaR|Chi2Ndf|Nhits|Weight|WeightTof400|WeightTof700|FhcalX|FhcalY|StsNhits|StsChi2)"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", {0, 10, 20, 30, 40, 50}} );
  correction_task.AddEventAxis( {"runId", {6666.5, 6667.5, 6668.5, 6669.5, 6670.5, 6671.5, 6672.5, 6673.5, 6674.5, 6675.5, 6676.5, 6677.5, 6678.5, 6679.5, 6680.5, 6682.0, 6683.5, 6684.5, 6686.0, 6688.0, 6690.5, 6691.5, 6692.5, 6693.0, 6694.5, 6695.5, 6697.0, 6698.5, 6715.5, 6731.5, 6732.5, 6733.5, 6735.5, 6737.5, 6738.5, 6739.5, 6742.5, 6748.5, 6752.5, 6756.5, 6760.5, 6763.0, 6765.5, 6766.5, 6767.5, 6768.5, 6770.0, 6771.5, 6772.5, 6773.5, 6776.5, 6779.5, 6781.0, 6782.5, 6784.0, 6785.5, 6787.0, 6790.5, 6794.5, 6796.0, 6797.5, 6798.0, 6799.5, 6801.5, 6803.5, 6815.5, 6816.5, 6817.5, 6818.5, 6819.5, 6820.5, 6821.5, 6822.5, 6879.5, 6881.0, 6882.5, 6883.5, 6885.0, 6886.5, 6888.0, 6890.0, 6895.5, 6900.5, 6901.5, 6902.5, 6903.5, 6904.5, 6905.5, 6906.5, 6907.5, 6908.5, 6909.5, 6910.5, 6913.0, 6915.5, 6917.0, 6918.5, 6919.5, 6920.5, 6922.0, 6923.5, 6925.0, 6926.5, 6927.5, 6928.5, 6929.5, 6930.5, 6931.5, 6932.5, 6933.5, 6934.5, 6935.5, 6936.5, 6938.0, 6939.5, 6954.0, 6967.5, 6968.5, 6970.0, 6971.5, 6972.5, 6973.5, 6974.5, 6975.5, 6976.5, 6977.5, 6978.5, 6979.5, 6980.5, 6981.5, 6982.5, 6983.5, 6987.0, 6990.5, 6991.5, 6992.5, 6993.5, 6994.5, 6996.0, 6997.5, 6998.5, 6999.5, 7001.0, 7002.5, 7003.5, 7004.5, 7005.5, 7007.0, 7008.5, 7009.5, 7010.5, 7011.5, 7012.5, 7021.0, 7030.5, 7031.5, 7032.5, 7033.5, 7034.5, 7036.0, 7037.5, 7039.0, 7040.5, 7041.5, 7042.5, 7043.5, 7045.0, 7046.5, 7047.5, 7048.5, 7049.5, 7050.5, 7051.5, 7052.5, 7053.5, 7054.5, 7055.5, 7056.5, 7065.5, 7075.5, 7076.5, 7077.5, 7078.5, 7080.0, 7081.5, 7082.5, 7083.5, 7085.0, 7086.5, 7088.0, 7091.5, 7092.5, 7093.5, 7095.0, 7096.5, 7097.5, 7098.5, 7101.5, 7102.5, 7103.5, 7104.5, 7114.5, 7125.5, 7126.5, 7127.5, 7128.5, 7129.5, 7130.5, 7131.5, 7132.5, 7134.0, 7135.5, 7136.5, 7137.5, 7141.5, 7146.5, 7147.5, 7148.0, 7149.5, 7150.5, 7152.5, 7154.5, 7155.5, 7156.5, 7158.0, 7159.5, 7160.5, 7161.5, 7162.5, 7163.5, 7164.5, 7165.5, 7166.5, 7167.5, 7170.5, 7173.5, 7174.5, 7175.5, 7176.5, 7177.5, 7178.5, 7179.5, 7180.5, 7181.5, 7183.0, 7184.5, 7185.0, 7186.5, 7187.5, 7189.0, 7191.5, 7192.5, 7193.5, 7194.5, 7195.5, 7197.5, 7200.5, 7201.0, 7202.5, 7204.0, 7205.5, 7206.5, 7207.5, 7208.5, 7210.0, 7211.5, 7212.5, 7213.5, 7214.5, 7215.5, 7216.5, 7217.5, 7218.5, 7219.5, 7221.5, 7223.5, 7224.0, 7240.0, 7255.5, 7256.5, 7258.5, 7259.5, 7261.5, 7262.0, 7263.5, 7264.5, 7265.5, 7266.5, 7267.5, 7268.5, 7269.5, 7270.5, 7271.5, 7273.0, 7274.5, 7275.0, 7276.5, 7277.5, 7278.5, 7280.0, 7281.5, 7282.5, 7284.5, 7285.0, 7286.5, 7287.0, 7288.5, 7289.5, 7290.5, 7301.5, 7311.5, 7312.5, 7316.5, 7320.5, 7321.5, 7322.5, 7323.5, 7324.0, 7325.5, 7326.5, 7327.5, 7328.5, 7329.0, 7336.5, 7337.5, 7339.5, 7342.5, 7343.5, 7344.5, 7345.5, 7346.5, 7347.0, 7348.5, 7349.5, 7350.5, 7351.5, 7352.5, 7353.5, 7354.5, 7355.5, 7356.5, 7357.5, 7358.5, 7360.0, 7361.5, 7362.0, 7363.5, 7364.5, 7365.5, 7366.0, 7367.5, 7368.0, 7369.5, 7371.5, 7374.5, 7375.0, 7376.5, 7377.5, 7378.5, 7379.5, 7380.5, 7381.5, 7382.5, 7384.0, 7386.5, 7387.5, 7388.5, 7389.5, 7390.5, 7391.5, 7392.5, 7393.5, 7394.0, 7395.5, 7396.5, 7397.5, 7398.5, 7399.5, 7400.5, 7401.5, 7402.5, 7404.0, 7405.5, 7406.5, 7407.0, 7408.5, 7409.5, 7410.5, 7411.5, 7412.5, 7413.5, 7414.5, 7416.0, 7417.5, 7418.5, 7420.0, 7421.5, 7422.5, 7424.0, 7425.5, 7426.0, 7427.5, 7428.5, 7430.0, 7431.5, 7432.5, 7433.5, 7434.5, 7436.0, 7437.5, 7438.0, 7439.5, 7440.5, 7441.5, 7443.0, 7444.5, 7445.5, 7446.5, 7448.0, 7449.5, 7450.0, 7451.5, 7452.5, 7453.5, 7454.5, 7455.5, 7456.5, 7457.5, 7459.0, 7460.5, 7465.0, 7469.5, 7470.5, 7471.5, 7472.5, 7473.5, 7475.5, 7477.5, 7479.0, 7480.5, 7481.5,7482.5, 7483.5, 7485.5, 7487.5, 7488.5, 7489.5, 7490.5, 7491.5, 7492.5, 7494.0, 7495.5, 7496.0, 7497.5, 7499.0, 7500.5, 7501.5, 7507.5, 7513.5, 7514.5, 7516.0, 7517.5, 7518.0, 7519.5, 7520.5, 7524.5, 7528.5, 7529.5, 7530.5, 7531.5, 7532.5, 7533.5, 7535.5, 7537.5, 7538.5, 7540.5, 7542.5, 7543.5, 7544.0, 7545.5, 7546.5, 7548.0, 7549.5, 7550.5, 7551.5, 7552.5, 7553.5, 7558.5, 7564.5, 7565.5, 7566.5, 7568.0, 7569.5, 7571.0, 7572.5, 7573.5, 7574.5, 7576.0, 7577.5, 7578.0, 7579.5, 7580.5, 7581.5, 7582.5, 7584.5, 7585.5, 7586.5, 7587.5, 7588.5, 7590.5, 7591.5, 7594.0, 7596.5, 7597.5, 7598.0, 7599.5, 7600.5, 7602.0, 7604.5, 7605.5, 7606.5, 7607.5, 7608.5, 7609.5, 7610.0, 7611.5, 7612.5, 7617.5, 7622.5, 7623.5, 7624.0, 7625.5, 7626.5, 7627.5, 7628.5, 7629.0, 7630.5, 7631.5, 7632.0, 7633.5, 7634.5, 7635.5, 7636.5, 7637.0, 7638.5, 7639.5, 7640.5, 7641.5, 7642.0, 7643.5, 7644.5, 7645.5, 7646.5, 7648.0, 7649.5, 7652.0, 7655.5, 7656.5, 7658.0, 7659.5, 7661.0, 7662.5, 7663.5, 7664.5, 7665.5, 7667.0, 7668.5, 7669.5, 7670.5, 7672.0, 7673.5, 7674.5, 7675.5, 7676.5, 7677.5, 7678.5, 7680.0, 7681.5, 7683.0, 7684.5, 7686.0, 7687.5, 7688.5, 7689.5, 7691.0, 7692.5, 7693.5, 7695.0, 7696.5, 7697.0, 7698.5, 7699.0, 7700.5, 7701.5, 7702.5, 7703.5, 7704.5, 7706.0, 7707.5, 7710.5, 7711.0, 7712.5, 7713.5, 7714.5, 7715.5, 7716.5, 7717.5, 7718.5, 7719.5, 7721.5, 7722.0, 7723.5, 7724.5, 7725.5, 7726.5, 7727.5, 7728.5, 7729.5, 7730.5, 7731.0, 7732.5, 7733.5, 7734.5, 7736.0, 7737.5, 7744.0, 7751.5, 7752.5, 7754.0, 7755.5, 7756.5, 7758.5, 7761.5, 7762.5, 7763.5, 7764.5, 7765.0, 7766.5, 7767.5, 7768.5, 7770.0, 7771.5, 7773.5, 7775.5, 7776.5, 7777.0, 7778.5, 7779.5, 7780.5, 7781.5, 7782.0, 7783.5, 7784.5, 7785.5, 7787.0, 7788.5, 7789.5, 7790.5, 7792.5, 7794.5, 7795.5, 7796.5, 7797.5, 7798.5, 7799.5, 7801.5, 7802.5, 7803.5, 7808.5, 7814.5, 7815.0, 7816.5, 7817.5, 7819.5, 7820.0, 7821.5, 7822.5, 7824.5, 7825.5, 7826.5, 7828.5, 7829.5, 7830.5, 7831.5, 7833.0, 7834.5, 7835.5, 7838.5, 7842.5, 7844.0, 7845.5, 7846.5, 7848.0, 7849.5, 7854.5, 7861.5, 7865.5, 7866.5, 7868.5, 7869.5, 7871.0, 7872.5, 7873.5, 7874.5, 7875.0, 7876.5, 7877.5, 7878.5, 7879.0, 7880.5, 7881.0, 7882.5, 7883.5, 7884.5, 7885.5, 7886.5, 7888.0, 7889.5, 7890.5, 7891.5, 7892.5, 7893.5, 7895.0, 7896.5, 7897.5, 7898.5, 7899.5, 7900.5, 7902.0, 7903.5, 7904.5, 7905.5, 7906.5, 7908.0, 7909.5, 7910.5, 7911.5, 7912.5, 7913.5, 7922.5, 7931.5, 7932.5, 7933.5, 7934.0, 7935.5, 7936.0, 7937.5, 7938.5, 7940.0, 7941.5, 7942.5, 7943.0, 7944.5, 7946.0, 7948.5, 7949.5, 7951.0, 7952.5, 7953.0, 7954.5, 7955.5, 7956.5, 7957.5, 7958.5, 7959.0, 7960.5, 7961.5, 7962.5, 7964.0, 7965.5, 7966.5, 7967.5, 7971.0, 7975.5, 7976.0, 7977.5, 7978.5, 7980.0, 7981.5, 7982.5, 7984.0, 7986.5, 7987.0, 7988.5, 7989.5, 7990.5, 7991.5, 7993.5, 7995.5, 7996.5, 7997.5, 7998.5, 7999.5, 8000.5, 8001.5, 8003.0, 8004.5, 8005.5, 8006.5, 8007.5, 8008.5, 8010.5, 8013.5, 8014.5, 8015.5, 8017.0, 8018.5, 8019.0, 8020.5, 8021.5, 8022.5, 8024.5, 8026.5, 8027.5, 8028.5, 8029.5, 8030.5, 8031.5, 8032.5, 8035.5, 8038.5, 8039.5, 8040.5, 8041.5, 8043.0, 8044.5, 8045.5, 8046.5, 8047.5, 8049.0, 8050.5, 8051.5, 8052.5, 8054.0, 8055.5, 8056.5, 8057.5, 8058.5, 8060.0, 8061.5, 8062.0, 8063.5, 8064.5, 8065.5, 8067.0, 8068.5, 8069.5, 8070.5, 8071.5, 8073.0, 8074.5, 8075.5, 8076.5, 8078.0, 8079.5, 8080.5, 8081.5, 8083.0, 8084.5, 8085.0, 8086.5, 8087.5, 8088.5, 8089.5,8093.5, 8097.5, 8098.5, 8100.5, 8101.5, 8102.5, 8103.0, 8104.5, 8105.0, 8106.5, 8107.0, 8108.5, 8109.5, 8110.5, 8111.5, 8112.5, 8114.0, 8115.5, 8116.5, 8117.5, 8118.5, 8120.0, 8121.5, 8122.5, 8123.5, 8126.5, 8129.5, 8130.5, 8132.0, 8133.5, 8135.0, 8137.5, 8138.5, 8139.5, 8140.5, 8141.5, 8143.0, 8144.5, 8150.0, 8156.5, 8157.5, 8158.5, 8159.5, 8160.5, 8161.5, 8163.5, 8165.5, 8166.5, 8168.0, 8169.5, 8171.5, 8173.5, 8174.5, 8175.5, 8176.5, 8178.5, 8180.5, 8181.5, 8183.5, 8184.5, 8185.0, 8186.5, 8187.0, 8188.5, 8189.0, 8190.5, 8191.5, 8192.5, 8194.0, 8195.5, 8197.0, 8198.5, 8200.0, 8201.5, 8202.5, 8203.5, 8204.5, 8205.5, 8206.5, 8207.5, 8208.5, 8209.5, 8210.5, 8211.5, 8212.5, 8214.0, 8215.5, 8216.0, 8217.5, 8218.0, 8219.5, 8220.5, 8224.5, 8228.5, 8229.5, 8230.5, 8233.0, 8235.5, 8236.5, 8237.0, 8238.5, 8239.5, 8241.0, 8242.5, 8243.0, 8244.5, 8245.5, 8246.5, 8247.5, 8248.5, 8249.0, 8250.5, 8251.5, 8252.0, 8253.5, 8254.5, 8255.5, 8256.5, 8257.5, 8261.5, 8265.5, 8266.5, 8267.5, 8269.0, 8270.5, 8271.5, 8272.0, 8273.5, 8274.5, 8275.5, 8276.5, 8277.5, 8278.5, 8280.0, 8281.5, 8282.5, 8284.5, 8285.0, 8286.5, 8287.5, 8288.5, 8289.5, 8290.5, 8291.0, 8292.5, 8293.5, 8294.5, 8296.0, 8297.5, 8298.5, 8299.5, 8302.5, 8305.5, 8306.5}} );
//correction_task.AddEventAxis( {"runId", {0,7300,7550,7650,7760,7930,8100,9000}} );
  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f1.AddCut( "fhcalModId", [&f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f2.AddCut( "fhcalModId", [&f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f3.AddCut( "fhcalModId", [&f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  // VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  // Tneg.SetHarmonicArray( {1, 2} );
  // Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  // Tneg.AddCut( "trCharge", [](double charge){
  //   return charge < 0.0;
  //   }, "charge" );
  // Tneg.AddCut( "trEta", [](double eta){
  //   return 1.5 < eta && eta < 4.0;
  //   }, "eta cut" );
  // Tneg.AddCut( "trPt", [](double pT){
  //   return pT > 0.2;
  //   }, "pT cut" );
  // Tneg.AddCut( "trFhcalX", [](double pos){
  //   return pos < -40.0 || pos > 170;
  //   }, "cut on x-pos in fhcal plane" );
  // Tneg.AddCut( "trFhcalY", [](double pos){
  //   return pos < -100.0 || pos > 100;
  //   }, "cut on y-pos in fhcal plane" );
  // correction_task.AddVector(Tneg);

  // VectorConfig Tpos( "Tpos", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  // Tpos.SetHarmonicArray( {1, 2, 3} );
  // Tpos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  // Tpos.AddCut( "trCharge", [](double charge){
  //   return charge >= 0.0;
  //   }, "charge" );
  // Tpos.AddCut( "trEta", [](double eta){
  //   return 2.0 < eta && eta < 3.0;
  // }, "eta cut" );
  // Tpos.AddCut( "trPt", [](double pT){
  //   return pT > 0.2;
  // }, "pT cut" );
  // Tpos.AddCut( "trFhcalX", [](double pos){
  //   return pos < -40.0 || pos > 170;
  //   }, "cut on x-pos in fhcal plane" );
  // Tpos.AddCut( "trFhcalY", [](double pos){
  //   return pos < -100.0 || pos > 100;
  //   }, "cut on y-pos in fhcal plane" );
  // correction_task.AddVector(Tpos);

  std::vector<Qn::AxisD> proton_axes{
        { "trProtonY", 4, -0.2, 1.4 },
        { "trPt", 5, 0.0, 2.0 },
  };
  
  VectorConfig proton( "proton", "trPhi", "trWeight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2, 3} );
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING,CORRECTION::TWIST_RESCALING  } );
  proton.SetCorrectionAxes( proton_axes );
//  proton.SetAlignmentReference( "F2" );
  proton.AddCut( "trNsigmaProton", [](double n_sigma){
    return n_sigma < 3;
  }, "proton cut" );
  proton.AddCut( "trFhcalX", [](double pos){
    return pos < -30.0 || pos > 160;
  }, "cut on x-pos in fhcal plane" );
  proton.AddCut( "trFhcalY", [](double pos){
    return pos < -60.0 || pos > 60;
  }, "cut on y-pos in fhcal plane" );
  proton.AddCut( "trStsNhits", [](double nhits){
    return nhits > 2.5;
  }, "cut on fake tracks" );
  proton.AddCut( "trDcaR", [](double dca){
    return dca < 5.0;
  }, "DCA cut" );
  proton.AddCut( "trStsChi2", [](double chi2){
    return chi2 < 5.0;
  }, "cut on chi2 in sts" );
  proton.AddHisto2D({{"trProtonY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}});
  proton.AddHisto2D({{"runId",{6666.5, 6667.5, 6668.5, 6669.5, 6670.5, 6671.5, 6672.5, 6673.5, 6674.5, 6675.5, 6676.5, 6677.5, 6678.5, 6679.5, 6680.5, 6682.0, 6683.5, 6684.5, 6686.0, 6688.0, 6690.5, 6691.5, 6692.5, 6693.0, 6694.5, 6695.5, 6697.0, 6698.5, 6715.5, 6731.5, 6732.5, 6733.5, 6735.5, 6737.5, 6738.5, 6739.5, 6742.5, 6748.5, 6752.5, 6756.5, 6760.5, 6763.0, 6765.5, 6766.5, 6767.5, 6768.5, 6770.0, 6771.5, 6772.5, 6773.5, 6776.5, 6779.5, 6781.0, 6782.5, 6784.0, 6785.5, 6787.0, 6790.5, 6794.5, 6796.0, 6797.5, 6798.0, 6799.5, 6801.5, 6803.5, 6815.5, 6816.5, 6817.5, 6818.5, 6819.5, 6820.5, 6821.5, 6822.5, 6879.5, 6881.0, 6882.5, 6883.5, 6885.0, 6886.5, 6888.0, 6890.0, 6895.5, 6900.5, 6901.5, 6902.5, 6903.5, 6904.5, 6905.5, 6906.5, 6907.5, 6908.5, 6909.5, 6910.5, 6913.0, 6915.5, 6917.0, 6918.5, 6919.5, 6920.5, 6922.0, 6923.5, 6925.0, 6926.5, 6927.5, 6928.5, 6929.5, 6930.5, 6931.5, 6932.5, 6933.5, 6934.5, 6935.5, 6936.5, 6938.0, 6939.5, 6954.0, 6967.5, 6968.5, 6970.0, 6971.5, 6972.5, 6973.5, 6974.5, 6975.5, 6976.5, 6977.5, 6978.5, 6979.5, 6980.5, 6981.5, 6982.5, 6983.5, 6987.0, 6990.5, 6991.5, 6992.5, 6993.5, 6994.5, 6996.0, 6997.5, 6998.5, 6999.5, 7001.0, 7002.5, 7003.5, 7004.5, 7005.5, 7007.0, 7008.5, 7009.5, 7010.5, 7011.5, 7012.5, 7021.0, 7030.5, 7031.5, 7032.5, 7033.5, 7034.5, 7036.0, 7037.5, 7039.0, 7040.5, 7041.5, 7042.5, 7043.5, 7045.0, 7046.5, 7047.5, 7048.5, 7049.5, 7050.5, 7051.5, 7052.5, 7053.5, 7054.5, 7055.5, 7056.5, 7065.5, 7075.5, 7076.5, 7077.5, 7078.5, 7080.0, 7081.5, 7082.5, 7083.5, 7085.0, 7086.5, 7088.0, 7091.5, 7092.5, 7093.5, 7095.0, 7096.5, 7097.5, 7098.5, 7101.5, 7102.5, 7103.5, 7104.5, 7114.5, 7125.5, 7126.5, 7127.5, 7128.5, 7129.5, 7130.5, 7131.5, 7132.5, 7134.0, 7135.5, 7136.5, 7137.5, 7141.5, 7146.5, 7147.5, 7148.0, 7149.5, 7150.5, 7152.5, 7154.5, 7155.5, 7156.5, 7158.0, 7159.5, 7160.5, 7161.5, 7162.5, 7163.5, 7164.5, 7165.5, 7166.5, 7167.5, 7170.5, 7173.5, 7174.5, 7175.5, 7176.5, 7177.5, 7178.5, 7179.5, 7180.5, 7181.5, 7183.0, 7184.5, 7185.0, 7186.5, 7187.5, 7189.0, 7191.5, 7192.5, 7193.5, 7194.5, 7195.5, 7197.5, 7200.5, 7201.0, 7202.5, 7204.0, 7205.5, 7206.5, 7207.5, 7208.5, 7210.0, 7211.5, 7212.5, 7213.5, 7214.5, 7215.5, 7216.5, 7217.5, 7218.5, 7219.5, 7221.5, 7223.5, 7224.0, 7240.0, 7255.5, 7256.5, 7258.5, 7259.5, 7261.5, 7262.0, 7263.5, 7264.5, 7265.5, 7266.5, 7267.5, 7268.5, 7269.5, 7270.5, 7271.5, 7273.0, 7274.5, 7275.0, 7276.5, 7277.5, 7278.5, 7280.0, 7281.5, 7282.5, 7284.5, 7285.0, 7286.5, 7287.0, 7288.5, 7289.5, 7290.5, 7301.5, 7311.5, 7312.5, 7316.5, 7320.5, 7321.5, 7322.5, 7323.5, 7324.0, 7325.5, 7326.5, 7327.5, 7328.5, 7329.0, 7336.5, 7337.5, 7339.5, 7342.5, 7343.5, 7344.5, 7345.5, 7346.5, 7347.0, 7348.5, 7349.5, 7350.5, 7351.5, 7352.5, 7353.5, 7354.5, 7355.5, 7356.5, 7357.5, 7358.5, 7360.0, 7361.5, 7362.0, 7363.5, 7364.5, 7365.5, 7366.0, 7367.5, 7368.0, 7369.5, 7371.5, 7374.5, 7375.0, 7376.5, 7377.5, 7378.5, 7379.5, 7380.5, 7381.5, 7382.5, 7384.0, 7386.5, 7387.5, 7388.5, 7389.5, 7390.5, 7391.5, 7392.5, 7393.5, 7394.0, 7395.5, 7396.5, 7397.5, 7398.5, 7399.5, 7400.5, 7401.5, 7402.5, 7404.0, 7405.5, 7406.5, 7407.0, 7408.5, 7409.5, 7410.5, 7411.5, 7412.5, 7413.5, 7414.5, 7416.0, 7417.5, 7418.5, 7420.0, 7421.5, 7422.5, 7424.0, 7425.5, 7426.0, 7427.5, 7428.5, 7430.0, 7431.5, 7432.5, 7433.5, 7434.5, 7436.0, 7437.5, 7438.0, 7439.5, 7440.5, 7441.5, 7443.0, 7444.5, 7445.5, 7446.5, 7448.0, 7449.5, 7450.0, 7451.5, 7452.5, 7453.5, 7454.5, 7455.5, 7456.5, 7457.5, 7459.0, 7460.5, 7465.0, 7469.5, 7470.5, 7471.5, 7472.5, 7473.5, 7475.5, 7477.5, 7479.0, 7480.5, 7481.5,7482.5, 7483.5, 7485.5, 7487.5, 7488.5, 7489.5, 7490.5, 7491.5, 7492.5, 7494.0, 7495.5, 7496.0, 7497.5, 7499.0, 7500.5, 7501.5, 7507.5, 7513.5, 7514.5, 7516.0, 7517.5, 7518.0, 7519.5, 7520.5, 7524.5, 7528.5, 7529.5, 7530.5, 7531.5, 7532.5, 7533.5, 7535.5, 7537.5, 7538.5, 7540.5, 7542.5, 7543.5, 7544.0, 7545.5, 7546.5, 7548.0, 7549.5, 7550.5, 7551.5, 7552.5, 7553.5, 7558.5, 7564.5, 7565.5, 7566.5, 7568.0, 7569.5, 7571.0, 7572.5, 7573.5, 7574.5, 7576.0, 7577.5, 7578.0, 7579.5, 7580.5, 7581.5, 7582.5, 7584.5, 7585.5, 7586.5, 7587.5, 7588.5, 7590.5, 7591.5, 7594.0, 7596.5, 7597.5, 7598.0, 7599.5, 7600.5, 7602.0, 7604.5, 7605.5, 7606.5, 7607.5, 7608.5, 7609.5, 7610.0, 7611.5, 7612.5, 7617.5, 7622.5, 7623.5, 7624.0, 7625.5, 7626.5, 7627.5, 7628.5, 7629.0, 7630.5, 7631.5, 7632.0, 7633.5, 7634.5, 7635.5, 7636.5, 7637.0, 7638.5, 7639.5, 7640.5, 7641.5, 7642.0, 7643.5, 7644.5, 7645.5, 7646.5, 7648.0, 7649.5, 7652.0, 7655.5, 7656.5, 7658.0, 7659.5, 7661.0, 7662.5, 7663.5, 7664.5, 7665.5, 7667.0, 7668.5, 7669.5, 7670.5, 7672.0, 7673.5, 7674.5, 7675.5, 7676.5, 7677.5, 7678.5, 7680.0, 7681.5, 7683.0, 7684.5, 7686.0, 7687.5, 7688.5, 7689.5, 7691.0, 7692.5, 7693.5, 7695.0, 7696.5, 7697.0, 7698.5, 7699.0, 7700.5, 7701.5, 7702.5, 7703.5, 7704.5, 7706.0, 7707.5, 7710.5, 7711.0, 7712.5, 7713.5, 7714.5, 7715.5, 7716.5, 7717.5, 7718.5, 7719.5, 7721.5, 7722.0, 7723.5, 7724.5, 7725.5, 7726.5, 7727.5, 7728.5, 7729.5, 7730.5, 7731.0, 7732.5, 7733.5, 7734.5, 7736.0, 7737.5, 7744.0, 7751.5, 7752.5, 7754.0, 7755.5, 7756.5, 7758.5, 7761.5, 7762.5, 7763.5, 7764.5, 7765.0, 7766.5, 7767.5, 7768.5, 7770.0, 7771.5, 7773.5, 7775.5, 7776.5, 7777.0, 7778.5, 7779.5, 7780.5, 7781.5, 7782.0, 7783.5, 7784.5, 7785.5, 7787.0, 7788.5, 7789.5, 7790.5, 7792.5, 7794.5, 7795.5, 7796.5, 7797.5, 7798.5, 7799.5, 7801.5, 7802.5, 7803.5, 7808.5, 7814.5, 7815.0, 7816.5, 7817.5, 7819.5, 7820.0, 7821.5, 7822.5, 7824.5, 7825.5, 7826.5, 7828.5, 7829.5, 7830.5, 7831.5, 7833.0, 7834.5, 7835.5, 7838.5, 7842.5, 7844.0, 7845.5, 7846.5, 7848.0, 7849.5, 7854.5, 7861.5, 7865.5, 7866.5, 7868.5, 7869.5, 7871.0, 7872.5, 7873.5, 7874.5, 7875.0, 7876.5, 7877.5, 7878.5, 7879.0, 7880.5, 7881.0, 7882.5, 7883.5, 7884.5, 7885.5, 7886.5, 7888.0, 7889.5, 7890.5, 7891.5, 7892.5, 7893.5, 7895.0, 7896.5, 7897.5, 7898.5, 7899.5, 7900.5, 7902.0, 7903.5, 7904.5, 7905.5, 7906.5, 7908.0, 7909.5, 7910.5, 7911.5, 7912.5, 7913.5, 7922.5, 7931.5, 7932.5, 7933.5, 7934.0, 7935.5, 7936.0, 7937.5, 7938.5, 7940.0, 7941.5, 7942.5, 7943.0, 7944.5, 7946.0, 7948.5, 7949.5, 7951.0, 7952.5, 7953.0, 7954.5, 7955.5, 7956.5, 7957.5, 7958.5, 7959.0, 7960.5, 7961.5, 7962.5, 7964.0, 7965.5, 7966.5, 7967.5, 7971.0, 7975.5, 7976.0, 7977.5, 7978.5, 7980.0, 7981.5, 7982.5, 7984.0, 7986.5, 7987.0, 7988.5, 7989.5, 7990.5, 7991.5, 7993.5, 7995.5, 7996.5, 7997.5, 7998.5, 7999.5, 8000.5, 8001.5, 8003.0, 8004.5, 8005.5, 8006.5, 8007.5, 8008.5, 8010.5, 8013.5, 8014.5, 8015.5, 8017.0, 8018.5, 8019.0, 8020.5, 8021.5, 8022.5, 8024.5, 8026.5, 8027.5, 8028.5, 8029.5, 8030.5, 8031.5, 8032.5, 8035.5, 8038.5, 8039.5, 8040.5, 8041.5, 8043.0, 8044.5, 8045.5, 8046.5, 8047.5, 8049.0, 8050.5, 8051.5, 8052.5, 8054.0, 8055.5, 8056.5, 8057.5, 8058.5, 8060.0, 8061.5, 8062.0, 8063.5, 8064.5, 8065.5, 8067.0, 8068.5, 8069.5, 8070.5, 8071.5, 8073.0, 8074.5, 8075.5, 8076.5, 8078.0, 8079.5, 8080.5, 8081.5, 8083.0, 8084.5, 8085.0, 8086.5, 8087.5, 8088.5, 8089.5,8093.5, 8097.5, 8098.5, 8100.5, 8101.5, 8102.5, 8103.0, 8104.5, 8105.0, 8106.5, 8107.0, 8108.5, 8109.5, 8110.5, 8111.5, 8112.5, 8114.0, 8115.5, 8116.5, 8117.5, 8118.5, 8120.0, 8121.5, 8122.5, 8123.5, 8126.5, 8129.5, 8130.5, 8132.0, 8133.5, 8135.0, 8137.5, 8138.5, 8139.5, 8140.5, 8141.5, 8143.0, 8144.5, 8150.0, 8156.5, 8157.5, 8158.5, 8159.5, 8160.5, 8161.5, 8163.5, 8165.5, 8166.5, 8168.0, 8169.5, 8171.5, 8173.5, 8174.5, 8175.5, 8176.5, 8178.5, 8180.5, 8181.5, 8183.5, 8184.5, 8185.0, 8186.5, 8187.0, 8188.5, 8189.0, 8190.5, 8191.5, 8192.5, 8194.0, 8195.5, 8197.0, 8198.5, 8200.0, 8201.5, 8202.5, 8203.5, 8204.5, 8205.5, 8206.5, 8207.5, 8208.5, 8209.5, 8210.5, 8211.5, 8212.5, 8214.0, 8215.5, 8216.0, 8217.5, 8218.0, 8219.5, 8220.5, 8224.5, 8228.5, 8229.5, 8230.5, 8233.0, 8235.5, 8236.5, 8237.0, 8238.5, 8239.5, 8241.0, 8242.5, 8243.0, 8244.5, 8245.5, 8246.5, 8247.5, 8248.5, 8249.0, 8250.5, 8251.5, 8252.0, 8253.5, 8254.5, 8255.5, 8256.5, 8257.5, 8261.5, 8265.5, 8266.5, 8267.5, 8269.0, 8270.5, 8271.5, 8272.0, 8273.5, 8274.5, 8275.5, 8276.5, 8277.5, 8278.5, 8280.0, 8281.5, 8282.5, 8284.5, 8285.0, 8286.5, 8287.5, 8288.5, 8289.5, 8290.5, 8291.0, 8292.5, 8293.5, 8294.5, 8296.0, 8297.5, 8298.5, 8299.5, 8302.5, 8305.5, 8306.5}},{"trPhi",100,-8,8}});
  correction_task.AddVector(proton);

  std::cout << "Initialized" << std::endl;

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}
