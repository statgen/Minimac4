#ifndef MINIMAC4_PROG_ARGS_HPP
#define MINIMAC4_PROG_ARGS_HPP

#include "getopt_wrapper.hpp"

#include <savvy/reader.hpp>

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <fstream>
#include <iterator>

class prog_args : public getopt_wrapper
{
private:
  std::string ref_path_;
  std::string tar_path_;
  std::string map_path_;
  std::string out_path_ = "/dev/stdout";
  std::string prefix_; // deprecated
  std::string emp_out_path_;
  std::string sites_out_path_;
  savvy::file::format out_format_ = savvy::file::format::sav;
  std::uint8_t out_compression_ = 6;
  std::vector<std::string> fmt_fields_ = {"HDS"};
  std::unordered_set<std::string> sample_ids_;
  savvy::genomic_region reg_ = {""};
  std::size_t temp_buffer_ = 200;
  std::int64_t chunk_size_ = 20000000;
  std::int64_t overlap_ = 3000000;
  std::int16_t threads_ = 1;
  float min_r2_ = -1.f;
  float min_ratio_ = 1e-5f;
  float prob_threshold_ = 0.01f;
  float prob_threshold_s1_ = -1.f;
  float diff_threshold_ = 0.01f;
  float min_recom_ = 1e-5f;
  float error_param_ = 0.01f;
  bool all_typed_sites_ = false;
  bool update_m3vcf_ = false;
  bool pass_only_ = false;
  bool meta_ = false; // deprecated
  bool help_ = false;
  bool version_ = false;

public:
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  const std::string& ref_path() const { return ref_path_; }
  const std::string& tar_path() const { return tar_path_; }
  const std::string& map_path() const { return map_path_; }
  const std::string& out_path() const { return out_path_; }
  const std::string& emp_out_path() const { return emp_out_path_; }
  const std::string& sites_out_path() const { return sites_out_path_; }
  savvy::file::format out_format() const { return out_format_; }
  std::uint8_t out_compression() const { return out_compression_; }
  const std::vector<std::string>& fmt_fields() const { return fmt_fields_; }
  const std::unordered_set<std::string>& sample_ids() const { return sample_ids_; }
  const savvy::genomic_region& region() const { return reg_; }
  std::int64_t chunk_size() const { return chunk_size_; }
  std::int64_t overlap() const { return overlap_; }
  std::int16_t threads() const { return threads_; }
  std::size_t temp_buffer() const { return temp_buffer_ ; }
  float min_r2() const { return min_r2_; }
  float min_ratio() const { return min_ratio_; }
  float prob_threshold() const { return prob_threshold_; }
  float prob_threshold_s1() const { return prob_threshold_s1_; }
  float diff_threshold() const { return diff_threshold_; }
  float min_recom() const { return min_recom_; }
  float error_param() const { return error_param_; }
  bool all_typed_sites() const { return all_typed_sites_; }
  bool update_m3vcf() const { return update_m3vcf_; }
  bool pass_only() const { return pass_only_; }

  prog_args() :
    getopt_wrapper(
      "Usage: minimac4 [opts ...] <reference.msav> <target.{sav,bcf,vcf.gz}>",
      {
        {"all-typed-sites", no_argument, 0, 'a', "Include in the output sites that exist only in target VCF"},
        {"temp-buffer", required_argument, 0, 'b', "Number of samples to impute before writing to temporary files (default: 200)"},
        {"chunk", required_argument, 0, 'c', "Maximum chunk length in base pairs to impute at once (default: 20000000"},
        {"empirical-output", required_argument, 0, 'e', "Output path for empirical dosages"},
        {"help", no_argument, 0, 'h', "Print usage"},
        {"format", required_argument, 0, 'f', "Comma-separated list of format fields to generate (GT, HDS, DS, GP, or SD; default: HDS)"},
        {"map", required_argument, 0, 'm', "Genetic map file"},
        {"output", required_argument, 0, 'o', "Output path (default: /dev/stdout)"},
        {"output-format", required_argument, 0, 'O', "Output file format (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: bcf)"},
        //{"pass-only", no_argument, 0, 'p', "Only imports variants with FILTER column set to PASS"},
        {"region", required_argument, 0, 'r', "Genomic region to impute"},
        {"sites", required_argument, 0, 's', "Output path for sites-only file"},
        {"threads", required_argument, 0, 't', "Number of threads (default: 1)"},
        {"version", no_argument, 0, 'v', "Print version"},
        {"overlap", required_argument, 0, 'w', "Size (in base pairs) of overlap before and after impute region to use as input to HMM (default: 3000000)"},
        {"min-r2", required_argument, 0, '\x02', "Minimum estimated r-square for output variants"},
        {"min-ratio", required_argument, 0, '\x02', "Minimum ratio of number of target sites to reference sites (default: 0.00001)"},
        {"match-error", required_argument, 0, '\x02', "Error parameter for HMM match probabilities (default: 0.01)"},
        {"min-recom", required_argument, 0, '\x02', "Minimum recombination probability (default: 0.00001)"},
        {"prob-threshold", required_argument, 0, '\x02', "Probability threshold used for template selection"},
        {"prob-threshold-s1", required_argument, 0, '\x02', ""}, // "Probability threshold used for template selection in original state space"},
        {"diff-threshold", required_argument, 0, '\x02', "Probability diff threshold used in template selection"},
        {"sample-ids", required_argument, 0, '\x02', "Comma-separated list of sample IDs to subset from reference panel"},
        {"sample-ids-file", required_argument, 0, '\x02', "Text file containing sample IDs to subset from reference panel (one ID per line)"},
        {"update-m3vcf", no_argument, 0, '\x01', nullptr},
        // vvvv deprecated vvvv //
        {"allTypedSites", no_argument, 0, '\x01', nullptr},
        {"rsid", no_argument, 0, '\x01', nullptr},
        //{"passOnly", no_argument, 0, '\x01', nullptr},
        {"meta", no_argument, 0, '\x01', nullptr},
        {"haps", required_argument, 0, '\x02', nullptr},
        {"refHaps", required_argument, 0, '\x02', nullptr},
        {"prefix", required_argument, 0, '\x02', nullptr},
        {"mapFile", required_argument, 0, '\x02', nullptr},
        {"chr", required_argument, 0, '\x02', nullptr},
        {"start", required_argument, 0, '\x02', nullptr},
        {"end", required_argument, 0, '\x02', nullptr},
        {"window", required_argument, 0, '\x02', nullptr},
        {"ChunkOverlapMb", required_argument, 0, '\x02', nullptr},
        {"ChunkLengthMb", required_argument, 0, '\x02', nullptr},
        {"cpus", required_argument, 0, '\x02', nullptr},
        {"minRatio", no_argument, 0, '\x02', nullptr}
      })
  {
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, short_opt_string_.c_str(), long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'a':
        all_typed_sites_ = true;
        break;
      case 'b':
        temp_buffer_ = std::size_t(std::atoll(optarg ? optarg : ""));
        break;
      case 'c':
        chunk_size_ = std::atoll(optarg ? optarg : "");
        break;
      case 'e':
        emp_out_path_ = optarg ? optarg : "";
        break;
      case 'h':
        help_ = true;
        return true;
      case 'f':
        {
          fmt_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
          std::unordered_set<std::string> allowed = {"GT", "GP", "DS", "HDS", "SD"};
          for (auto it = fmt_fields_.begin(); it != fmt_fields_.end(); ++it)
          {
            if (allowed.find(*it) == allowed.end())
              return std::cerr << "Error: Invalid --format option (" << *it << ")\n", false;
          }
          break;
        }
      case 'm':
        map_path_ = optarg ? optarg : "";
        break;
      case 'o':
        out_path_ = optarg ? optarg : "";
        break;
      case 'O':
        {
          using fmt = savvy::file::format;
          std::string ot = optarg ? optarg : "";
          if (ot == "vcf")
          {
            out_format_ = fmt::vcf;
            out_compression_ = 0;
          }
          else if (ot == "vcf.gz")
          {
            out_format_ = fmt::vcf;
          }
          else if (ot == "bcf")
          {
            out_format_ = fmt::bcf;
          }
          else if (ot == "ubcf")
          {
            out_format_ = fmt::bcf;
            out_compression_ = 0;
          }
          else if (ot == "sav")
          {
            out_format_ = fmt::sav;
          }
          else if (ot == "usav")
          {
            out_format_ = fmt::sav;
            out_compression_ = 0;
          }
          else
          {
            std::cerr << "Invalid --output-format: " << ot << std::endl;
            return false;
          }
          break;
        }
      case 'p':
        pass_only_ = true;
        break;
      case 'r':
        reg_ = string_to_region(optarg ? optarg : "");
        break;
      case 's':
        sites_out_path_ = optarg ? optarg : "";
        break;
      case 't':
        threads_ = atoi(optarg ? optarg : "");
        break;
      case 'v':
        version_ = true;
        return true;
      case 'w':
        overlap_ = std::atoll(optarg ? optarg : "");
        break;
      case '\x01':
        if (std::string(long_options_[long_index].name) == "update-m3vcf")
        {
          update_m3vcf_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "allTypedSites")
        {
          std::cerr << "Warning: --allTypedSites is deprecated in favor of --all-typed-sites\n";
          all_typed_sites_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "rsid")
        {
          std::cerr << "Warning: --rsid is deprecated (on by default)\n";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "passOnly")
        {
          std::cerr << "Warning: --passOnly is deprecated in favor of --pass-only\n";
          pass_only_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "meta")
        {
          std::cerr << "Warning: --meta is deprecated in favor of --empirical-output\n";
          meta_ = true;
          break;
        }
        // else pass through to default
      case '\x02':
        {
          std::string long_opt_str = std::string(long_options_[long_index].name);
          if (long_opt_str == "min-r2")
          {
            min_r2_ = std::atof(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "min-ratio")
          {
            min_ratio_ = std::min(1., std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "match-error")
          {
            error_param_ = std::min(0.5, std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "min-recom")
          {
            min_recom_ = std::min(0.5, std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "prob-threshold")
          {
            prob_threshold_ = std::min(1., std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "prob-threshold-s1")
          {
            prob_threshold_s1_ = std::min(1., std::atof(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "diff-threshold")
          {
            diff_threshold_ = std::max(0., std::atof(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "sample-ids")
          {
            auto tmp_ids = split_string_to_vector(optarg ? optarg : "", ',');
            sample_ids_.insert(tmp_ids.begin(), tmp_ids.end());
            break;
          }
          else if (long_opt_str == "sample-ids-file")
          {
            std::ifstream ifs(optarg ? optarg : "");
            sample_ids_.insert(std::istream_iterator<std::string>(ifs), std::istream_iterator<std::string>());
            break;
          }
          else if (long_opt_str == "haps")
          {
            std::cerr << "Warning: --haps is deprecated\n";
            tar_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "refHaps")
          {
            std::cerr << "Warning: --refHaps is deprecated\n";
            ref_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "chr")
          {
            std::cerr << "Warning: --chr is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(optarg ? optarg : "", reg_.from(), reg_.to());
            break;
          }
          else if (long_opt_str == "start")
          {
            std::cerr << "Warning: --start is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(reg_.chromosome(), std::atoll(optarg ? optarg : ""), reg_.to());
            break;
          }
          else if (long_opt_str == "end")
          {
            std::cerr << "Warning: --end is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(reg_.chromosome(), reg_.from(), std::atoll(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "prefix")
          {
            std::cerr << "Warning: --prefix is deprecated in favor of --output, --empirical-output, and --sites\n";
            prefix_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "mapFile")
          {
            std::cerr << "Warning: --mapFile is deprecated in favor of --map\n";
            map_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "window")
          {
            std::cerr << "Warning: --window is deprecated in favor of --overlap\n";
            overlap_ = std::atoll(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "ChunkLengthMb")
          {
            std::cerr << "Warning: --ChunkLengthMb is deprecated in favor of --chunk\n";
            chunk_size_ = std::atoll(optarg ? optarg : "") * 1000000;
            break;
          }
          else if (long_opt_str == "ChunkOverlapMb")
          {
            std::cerr << "Warning: --ChunkOverlapMb is deprecated in favor of --overlap\n";
            overlap_ = std::atoll(optarg ? optarg : "") * 1000000;
            break;
          }
          else if (long_opt_str == "cpus")
          {
            std::cerr << "Warning: --cpus is deprecated in favor of --threads\n";
            threads_ = atoi(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "minRatio")
          {
            std::cerr << "Warning: --minRatio is deprecated in favor of --min-ratio\n";
            min_ratio_ = std::atof(optarg ? optarg : "");
            break;
          }
          // else pass through to default
        }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 2)
    {
      ref_path_ = argv[optind];
      tar_path_ = argv[optind + 1];
    }
    else if (update_m3vcf_ && remaining_arg_count == 1)
    {
      ref_path_ = argv[optind];
    }
    else if (remaining_arg_count < 2)
    {
      if (ref_path_.empty() || tar_path_.empty())
      {
        std::cerr << "Too few arguments\n";
        return false;
      }
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (!prefix_.empty())
    {
      std::string suffix = "sav";
      if (out_format_ == savvy::file::format::bcf)
        suffix = "bcf";
      else if (out_format_ == savvy::file::format::vcf)
      {
        suffix = "vcf";
        if (out_compression_)
          suffix += ".gz";
      }

      out_path_ = prefix_ + ".dose." + suffix;
      sites_out_path_ = prefix_ + ".sites." + suffix;
      if (meta_)
        emp_out_path_ = prefix_ + ".empiricalDose." + suffix;
    }

    if (!emp_out_path_.empty() && std::find(fmt_fields_.begin(), fmt_fields_.end(), "HDS") == fmt_fields_.end())
      fmt_fields_.emplace_back("HDS");

    return true;
  }
private:
  savvy::genomic_region string_to_region(const std::string& s)
  {
    const std::size_t colon_pos = s.find(':');
    if (colon_pos == std::string::npos)
    {
      return savvy::genomic_region(s);
    }
    else
    {
      std::string chr = s.substr(0, colon_pos);
      const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
      if (hyphen_pos == std::string::npos)
      {
        std::string slocus = s.substr(colon_pos + 1);
        std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
        return savvy::genomic_region(chr, ilocus, ilocus);
      }
      else
      {
        std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
        std::string send = s.substr(hyphen_pos + 1);
        if (send.empty())
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
        }
        else
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
        }
      }
    }
  }

  std::vector<std::string> split_string_to_vector(const char* in, char delim)
  {
    std::vector<std::string> ret;
    const char* d = nullptr;
    std::string token;
    const char* s = in;
    const char*const e = in + strlen(in);
    while ((d = std::find(s, e,  delim)) != e)
    {
      ret.emplace_back(std::string(s, d));
      s = d ? d + 1 : d;
    }
    ret.emplace_back(std::string(s,d));
    return ret;
  }
};

#endif // MINIMAC4_PROG_ARGS_HPP