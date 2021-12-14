#ifndef MINIMAC4_RECOMBINATION_HPP
#define MINIMAC4_RECOMBINATION_HPP

#include "variant.hpp"

#include <shrinkwrap/istream.hpp>

#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

class recombination
{
public:
  //static constexpr float recom_min = 1e-5f;
//private:
  struct map_file_line
  {
   std::string chrom;
   std::size_t pos = 0;
   float map_value = 0.;
  };
public:
  static bool parse_map_file(const std::string& map_file_path, std::vector<target_variant>& sites, float recom_min);
  static double haldane(double cm) { return (1. - std::exp(-cm/50.))/2.;}
  static double cm_to_switch_prob(double cm) { return 1. - std::exp(-cm/100.);}
  static double haldane_inverse(double recom_prob) { return 50. * std::log(1. / (1. - 2. * recom_prob)); }
  static double switch_prob_to_cm(double recom_prob) { return 100. * std::log(1. / (1. - recom_prob)); }
private:
  static bool read_entry(std::istream& ifs, map_file_line& entry, bool new_format);
};

class genetic_map_file
{
public:
  struct record
  {
    std::string chrom;
    std::size_t pos = 0;
    double map_value = 0.;
  };
private:
  shrinkwrap::istream ifs_;
  std::string target_chrom_;
  record prev_rec_;
  record cur_rec_;
  bool good_;
  bool new_format_;
public:
  genetic_map_file(const std::string& map_file_path, const std::string& chrom);

  bool good() const { return good_; }
  operator bool() const { return good_; }

  // This must be called in order. The pos argument must be >= to pos passed in previous call
  float interpolate_centimorgan(std::size_t variant_pos);
private:
  bool read_record(record& rec);
};

#endif // MINIMAC4_RECOMBINATION_HPP