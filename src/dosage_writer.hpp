#ifndef MINIMAC4_DOSAGE_WRITER_HPP
#define MINIMAC4_DOSAGE_WRITER_HPP

#include "variant.hpp"
#include "unique_haplotype.hpp"
#include "hidden_markov_model.hpp" // TODO: change imputed data class

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <memory>

class dosage_writer
{
private:
  struct accuracy_statistics
  {
    double er2_sum = 0.;
    std::size_t n_var = 0;
  };

  savvy::writer out_file_;
  std::unique_ptr<savvy::writer> emp_out_file_;
  std::unique_ptr<savvy::writer> sites_out_file_;
  savvy::file::format file_format_;
  std::vector<std::string> fmt_fields_;
  std::unordered_set<std::string> fmt_field_set_;
  std::vector<accuracy_statistics> accuracy_stats_;
  std::size_t n_samples_ = 0;
  float min_r2_ = -1.f;
  bool is_temp_file_;

  // buffers
  savvy::compressed_vector<std::int8_t> sparse_gt_;
  std::vector<float> dense_float_vec_;
  std::vector<float> dense_zero_vec_;

  struct variant_update_ctx
  {
    savvy::compressed_vector<std::int8_t> sparse_gt;
    savvy::compressed_vector<float> sparse_dosages;
    std::vector<float> gp_vec;
  };
public:
  dosage_writer(const std::string& file_path, const std::string& emp_file_path, const std::string& sites_file_path,
    savvy::file::format file_format,
    std::uint8_t out_compression,
    const std::vector<std::string>& sample_ids,
    const std::vector<std::string>& fmt_fields,
    const std::string& chromosome,
    float min_r2, bool is_temp);

  bool merge_temp_files(std::list<savvy::reader>& temp_files, std::list<savvy::reader>& temp_emp_files);
  bool merge_temp_files(std::list<std::string>& temp_file_paths, std::list<std::string>& temp_emp_file_paths);

  bool write_dosages(const full_dosages_results& hmm_results, const std::vector<target_variant>& tar_variants, const std::vector<target_variant>& tar_only_variants, std::pair<std::size_t, std::size_t> observed_range, const reduced_haplotypes& full_reference_data, const savvy::region& impute_region);

  void print_mean_er2(std::ostream& os) const;
private:
  static std::vector<std::pair<std::string, std::string>> gen_headers(const std::vector<std::string>& fmt_fields, const std::string& chromosome, bool is_temp);
  static std::vector<std::pair<std::string, std::string>> gen_emp_headers(const std::string& chromosome);

  static bool sites_match(const target_variant& t, const reference_site_info& r);
  bool has_good_r2(savvy::site_info& site);

  static float calc_er2(double s_x, double s_xx, double s_y, double s_yy, double s_xy, std::size_t n);
  static float calc_r2(double s_x, double s_xx, std::size_t n);

  void set_info_fields(savvy::variant& out_var, const savvy::compressed_vector<float>& sparse_dosages, const std::vector<float>& loo_dosages, const std::vector<std::int8_t>& observed);
  void set_format_fields(savvy::variant& out_var, savvy::compressed_vector<float>& sparse_dosages);


  struct plus_ignore_missing
  {
    float operator()(const float& l, const float& r)
    {
      if (std::isnan(r))
        return l;
      if (std::isnan(l))
        return r;
      return l + r;
    }

    std::int32_t operator()(const std::int32_t& l, const std::int32_t& r)
    {
      if (r < 0)
        return l;
      if (l < 0)
        return r;
      return l + r;
    }
  };
};

#endif // MINIMAC4_DOSAGE_WRITER_HPP
