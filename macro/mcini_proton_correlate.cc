#include "QnDataFrame.hpp"

vector <vector<string>> u1_psi_rp=
{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u1_rescaled_Q1_rescaled=
{
  {"tru_proton_PLAIN", "S1_PLAIN"},
  {"tru_proton_PLAIN", "S2_PLAIN"},
  {"tru_proton_PLAIN", "S3_PLAIN"},
};

vector <vector<string>> Q1_rescaled_Q1_rescaled =
{
  {"S1_PLAIN", "psi_rp_PLAIN"},
  {"S2_PLAIN", "psi_rp_PLAIN"},
  {"S3_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q1Q1_rescaled=
{
  {"tru_proton_PLAIN", "S1_PLAIN", "S1_PLAIN"},
  {"tru_proton_PLAIN", "S2_PLAIN", "S2_PLAIN"},
  {"tru_proton_PLAIN", "S3_PLAIN", "S3_PLAIN"},
  
  {"tru_proton_PLAIN", "S1_PLAIN", "S2_PLAIN"},
  {"tru_proton_PLAIN", "S1_PLAIN", "S3_PLAIN"},
  {"tru_proton_PLAIN", "S2_PLAIN", "S3_PLAIN"},
};

vector <vector<string>> u3Q1Q1Q1_rescaled=
{
  {"tru_proton_PLAIN", "S1_PLAIN", "S1_PLAIN", "S1_PLAIN"},
  {"tru_proton_PLAIN", "S2_PLAIN", "S2_PLAIN", "S2_PLAIN"},
  {"tru_proton_PLAIN", "S3_PLAIN", "S3_PLAIN", "S3_PLAIN"},
  
  {"tru_proton_PLAIN", "S1_PLAIN", "S2_PLAIN", "S3_PLAIN"},
};

namespace P13 {
  inline auto x3x3() {
    return [](const Qn::QVector &a, const Qn::QVector &b) {
      return pow( a.x(1)*b.x(1), 3 );
    };
  }
  inline auto y3y3() {
    return [](const Qn::QVector &a, const Qn::QVector &b) {
      return pow( a.y(1)*b.y(1), 3 );
    };
  }
}

namespace P4{

  inline auto xxxx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.x(h_b) * c.x(h_c) * d.x(h_d);
    };
  }

  inline auto xxyy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.x(h_b) * c.y(h_c) * d.y(h_d);
    };
  }

  inline auto xyxy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.y(h_b) * c.x(h_c) * d.y(h_d);
    };
  }


  inline auto xyyx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.y(h_b) * c.y(h_c) * d.x(h_d);
    };
  }

  inline auto yxxy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.x(h_b) * c.x(h_c) * d.y(h_d);
    };
  }

  inline auto yxyx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.x(h_b) * c.y(h_c) * d.x(h_d);
    };
  }

  inline auto yyxx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.y(h_b) * c.x(h_c) * d.x(h_d);
    };
  }

  inline auto yyyy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.y(h_b) * c.y(h_c) * d.y(h_d);
    };
  }
}

void mcini_proton_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 4, 0, 40});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wUnity3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return 1; };

  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wSumWu4part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};

  for ( auto &corr: u1_psi_rp )
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x3", P2::xx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y3", P2::yy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y3", P2::xy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x3", P2::yx(3, 3), wSumWu, wy, qn, qn);
  }

  for ( auto &corr: u1_rescaled_Q1_rescaled )
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1", P2::xx(3, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1", P2::yy(3, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1", P2::xy(3, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1", P2::yx(3, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr: Q1_rescaled_Q1_rescaled)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x3", P13::x3x3(), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y3", P13::y3y3(), wUnity, wn, qn, qn);
  }

  for ( auto &corr: u2Q1Q1_rescaled )
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, qn, qn);
  }

  for ( auto &corr: u3Q1Q1Q1_rescaled )
  {
    std::array<std::string, 4> qn{corr.at(0), corr.at(1), corr.at(2), corr.at(3)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2)+"."+corr.at(3);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1x1x1", P4::xxxx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1y1y1", P4::yyyy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1x1x1", P4::xxxx(1, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1y1y1", P4::yyyy(1, 1, 1, 1), wSumWu4part, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
}
