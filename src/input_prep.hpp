#ifndef MINIMAC4_INPUT_PREP_HPP
#define MINIMAC4_INPUT_PREP_HPP

#include "unique_haplotype.hpp"

#include <savvy/reader.hpp>

bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids);

bool load_reference_haplotypes(const std::string& file_path, const savvy::genomic_region& extended_reg, const savvy::genomic_region& impute_reg, std::vector<target_variant>& target_sites, reduced_haplotypes& typed_only_reference_data, reduced_haplotypes& full_reference_data);

std::vector<target_variant> separate_target_only_variants(std::vector<target_variant>& target_sites);

std::vector<std::vector<std::vector<std::size_t>>> generate_reverse_maps(const reduced_haplotypes& typed_only_reference_data);

bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path);


#endif // MINIMAC4_INPUT_PREP_HPP