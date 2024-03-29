#ifndef MINIMAC4_INPUT_PREP_HPP
#define MINIMAC4_INPUT_PREP_HPP

#include "unique_haplotype.hpp"

#include <savvy/reader.hpp>

bool stat_tar_panel(const std::string& tar_file_path, std::vector<std::string>& sample_ids);

bool stat_ref_panel(const std::string& ref_file_path, std::string& chrom, std::uint64_t& end_pos);

bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids);

bool load_reference_haplotypes(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data,
  genetic_map_file* map_file,
  float min_recom,
  float default_match_error);

bool load_reference_haplotypes_old_recom_approach(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data,
  genetic_map_file* map_file);

std::vector<target_variant> separate_target_only_variants(std::vector<target_variant>& target_sites);

bool load_variant_hmm_params(std::vector<target_variant>& tar_variants, reduced_haplotypes& typed_only_reference_data, float default_error_param, float recom_min, const std::string& map_file_path);

std::vector<std::vector<std::vector<std::size_t>>> generate_reverse_maps(const reduced_haplotypes& typed_only_reference_data);

bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path, const std::string& map_file_path = "");
bool compress_reference_panel(const std::string& input_path, const std::string& output_path,
  std::size_t min_block_size = 10,
  std::size_t max_block_size = 0xFFFF, // max s1r block size minus 1 partition record
  std::size_t slope_unit = 10, 
  const std::string& map_file_path = "");


#endif // MINIMAC4_INPUT_PREP_HPP
