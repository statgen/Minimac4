#ifndef MINIMAC4_DOSAGE_WRITER_HPP
#define MINIMAC4_DOSAGE_WRITER_HPP

#include "variant.hpp"
#include "unique_haplotype.hpp"
#include "hidden_markov_model.hpp" // TODO: change imputed data class

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

class dosage_writer
{
private:
  savvy::writer out_file_;
  savvy::file::format file_format_;
  std::vector<std::string> fmt_fields_;
  std::unordered_set<std::string> fmt_field_set_;
  std::size_t n_samples_ = 0;
  const std::int16_t bin_scalar_ = 100; //256;

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
  dosage_writer(const std::string& file_path, savvy::file::format file_format, std::uint8_t out_compression, const std::vector<std::string>& sample_ids, const std::vector<std::string>& fmt_fields, const std::string& chromosome);

  bool merge_temp_files(const std::vector<std::string>& temp_files_paths);

  bool write_dosages(const full_dosages_results& hmm_results, const std::vector<target_variant>& tar_variants, std::pair<std::size_t, std::size_t> observed_range, const reduced_haplotypes& full_reference_data);

  bool piecewise_constant( const best_templates_results& hmm_results, const std::vector<target_variant>& typed_variants, const reduced_haplotypes& reduced_reference, const std::string& reference_path, const savvy::genomic_region& reg, const std::string& output_path);
private:
  static std::vector<std::pair<std::string, std::string>> gen_headers(const std::vector<std::string>& fmt_fields, const std::string& chromosome);

  static bool sites_match(const target_variant& t, const reference_site_info& r);

  void prepare_output_variant(savvy::variant& out_var, variant_update_ctx& ctx, const std::vector<float>& dosages, const std::vector<float>& loo_dosages, const std::vector<std::int8_t>& observed);

  void generate_extra_format_fields(savvy::variant& out_var, savvy::compressed_vector<float>& sparse_dosages);
};

#endif // MINIMAC4_DOSAGE_WRITER_HPP
