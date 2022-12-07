#include "input_prep.hpp"
#include "recombination.hpp"

#include <algorithm>
#include <sys/stat.h>

bool stat_tar_panel(const std::string& tar_file_path, std::vector<std::string>& sample_ids)
{
  savvy::reader temp_rdr(tar_file_path);
  if (!temp_rdr)
    return std::cerr << "Error: could not open target file\n", false;

  sample_ids = temp_rdr.samples();

  return true;
}

bool stat_ref_panel(const std::string& ref_file_path, std::string& chrom, std::uint64_t& end_pos)
{
  std::string separate_s1r_path = ref_file_path + ".s1r";
  struct stat st;
  std::vector<savvy::s1r::index_statistics> s1r_stats = savvy::s1r::stat_index(stat(separate_s1r_path.c_str(), &st) == 0 ? separate_s1r_path : ref_file_path);
  if (s1r_stats.size())
  {
    if (chrom.size())
    {
      for (auto it = s1r_stats.begin(); it != s1r_stats.end(); ++it)
      {
        if (it->contig == chrom)
        {
          end_pos = std::min(end_pos, std::uint64_t(it->max_position));
          return true;
        }
      }

      std::cerr << "Error: reference file does not contain chromosome " << chrom << "\n";
      return false;
    }
    else if (s1r_stats.size() == 1)
    {
      chrom = s1r_stats.front().contig;
      end_pos = std::min(end_pos, std::uint64_t(s1r_stats.front().max_position));
      return true;
    }

    std::cerr << "Error: reference file contains multiple chromosomes so --region is required\n";
    return false;
  }
  else if (stat((ref_file_path + ".csi").c_str(), &st) == 0 || stat((ref_file_path + ".tbi").c_str(), &st) == 0)
  {
    savvy::reader stat_rdr(ref_file_path);
    if (chrom.empty())
    {
      savvy::variant var;
      stat_rdr >> var;
      chrom = var.chromosome();
    }

    if (chrom.size())
    {
      for (auto it = stat_rdr.headers().begin(); it != stat_rdr.headers().end(); ++it)
      {
        if (it->first == "contig" && chrom == savvy::parse_header_sub_field(it->second, "ID"))
        {
          std::string length_str = savvy::parse_header_sub_field(it->second, "length");
          if (length_str.size())
          {
            end_pos = std::min(end_pos, std::uint64_t(std::atoll(length_str.c_str())));
            return true;
          }
          break;
        }
      }
      std::cerr << "Error: could not parse chromosome length from VCF/BCF header so --region is required" << std::endl;
      return false;
    }

    std::cerr << "Error: could not determine chromosome from reference" << std::endl;
    return false;
  }

  std::cerr << "Error: could not load reference file index (reference must be an indexed MVCF)\n";
  std::cerr << "Notice: M3VCF files must be updated to an MVCF encoded file. This can be done by running `minimac4 --update-m3vcf input.m3vcf.gz > output.msav`\n";
  return false;
}

void init_ploidies(std::vector<std::uint8_t>& ploidies, const std::vector<std::int8_t>& gt_vec)
{
  std::uint8_t max_ploidy = gt_vec.size() / ploidies.size();
  std::fill(ploidies.begin(), ploidies.end(), max_ploidy);
  for (std::size_t i = 0; i < gt_vec.size(); ++i)
  {
    if (savvy::typed_value::is_end_of_vector(gt_vec[i]))
      --ploidies[i / max_ploidy];
  }
}

std::int64_t check_ploidies(const std::vector<std::uint8_t>& ploidies, const std::vector<std::int8_t>& gt_vec)
{
  std::uint8_t max_ploidy = gt_vec.size() / ploidies.size();
  for (std::size_t i = 0; i < ploidies.size(); ++i)
  {
    std::uint8_t p = max_ploidy;
    for (std::size_t j = 0; j < max_ploidy; ++j)
    {
      if (savvy::typed_value::is_end_of_vector(gt_vec[i * max_ploidy + j]))
        --p;
    }

    if (p != ploidies[i])
      return std::int64_t(i);
  }

  return -1;
}

bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids)
{
  savvy::reader input(file_path);
  if (!input)
    return std::cerr << "Error: cannot open target file\n", false;

  sample_ids = input.samples();
  input.reset_bounds(reg);
  if (!input)
    return std::cerr << "Error: cannot query region (" << reg.chromosome() << ":" << reg.from() << "-" << reg.to() << ") from target file. Target file must be indexed.\n", false;

  const auto nan = std::numeric_limits<float>::quiet_NaN();
  std::vector<std::uint8_t> ploidies(sample_ids.size());
  savvy::variant var;
  std::vector<std::int8_t> tmp_geno;
  while (input >> var)
  {
    var.get_format("GT", tmp_geno);

    std::int64_t ploidy_res = -1;
    if (ploidies[0] == 0)
      init_ploidies(ploidies, tmp_geno);
    else
      ploidy_res = check_ploidies(ploidies, tmp_geno);

    if (ploidy_res >= 0)
    {
      std::cerr << "Error: Sample " << sample_ids[ploidy_res] << " changes ploidy at " << var.chrom() << ":" << var.pos() << "\n";
      if (var.chrom() == "X" || var.chrom() == "chrX")
        std::cerr << "Notice: PAR and non-PAR regions on chromosome X should be imputed separately\n";
      return false;
    }

    for (std::size_t i = 0; i < var.alts().size(); ++i)
    {
      std::size_t allele_idx = i + 1;
      target_sites.push_back({var.chromosome(), var.position(), var.id(), var.ref(), var.alts()[i], true, false, nan, nan, nan, {}});
      if (var.alts().size() == 1)
        tmp_geno.swap(target_sites.back().gt);
      else
      {
        target_sites.back().gt.resize(tmp_geno.size());
        std::int8_t* dest_ptr = target_sites.back().gt.data();
        for (std::size_t j = 0; j < tmp_geno.size(); ++j)
          dest_ptr[j] = std::int8_t(tmp_geno[j] == allele_idx);
      }
    }
  }

  return !input.bad();
}

bool load_reference_haplotypes(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data)
{
  savvy::reader input(file_path);

  if (input)
  {
    if (!input.reset_bounds(extended_reg, savvy::bounding_point::any))
      return std::cerr << "Error: reference file must be indexed MVCF\n", false;

    bool is_m3vcf_v3 = false;
    for (auto it = input.headers().begin(); !is_m3vcf_v3 && it != input.headers().end(); ++it)
    {
      if (it->first == "subfileformat" && (it->second == "M3VCFv3.0" || it->second == "MVCFv3.0"))
        is_m3vcf_v3 = true;
    }

    if (!is_m3vcf_v3)
      return std::cerr << "Error: reference file must be an MVCF\n", false;

    if (subset_ids.size() && input.subset_samples(subset_ids).empty())
      return std::cerr << "Error: no reference samples overlap subset IDs\n", false;

    savvy::variant var;
    if (!input.read(var))
      return std::cerr << "Notice: no variant records in reference query region (" << extended_reg.chromosome() << ":" << extended_reg.from() << "-" << extended_reg.to() << ")\n", input.bad() ? false : true;

    double start_cm = 0.;
    std::vector<std::int8_t> tmp_geno;
    unique_haplotype_block block;
    auto tar_it = target_sites.begin();
    int res;
    while ((res = block.deserialize(input, var)) > 0)
    {
      if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
        break;

      block.fill_cm_from_recom(start_cm);

      for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
      {
        while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];

            typed_only_reference_data.compress_variant({ref_it->chrom, ref_it->pos, ref_it->id, ref_it->ref, ref_it->alt, ref_it->err, ref_it->recom, ref_it->cm}, tmp_geno);

            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;

            if (it != tar_it)
              std::swap(*it, *tar_it);

            ++tar_it;
            break;
          }
        }
      }

      block.trim(impute_reg.from(), impute_reg.to());
      if (!block.variants().empty())
        full_reference_data.append_block(block);
    }

    if (res < 0)
      return false;


    return !input.bad();
  }
  else
  {
    shrinkwrap::gz::istream input_file(file_path);
    std::string line;

    std::uint8_t m3vcf_version = 0;
    const std::string m3vcf_version_line = "##fileformat=M3VCF";
    while (std::getline(input_file, line))
    {
      if (line.substr(0, m3vcf_version_line.size()) == m3vcf_version_line)
      {
        if (line == "##fileformat=M3VCFv2.0")
          m3vcf_version = 2;
        else
          m3vcf_version = 1;
        break;
      }

      if (line.size() < 2 || line[1] != '#')
      {
        std::cerr << "Error: invalid reference file" << std::endl;
        return false;
      }
    }

    std::size_t n_samples = 0;
    while (std::getline(input_file, line))
    {
      if (line.size() < 2 || line[1] != '#')
      {
        n_samples = std::count(line.begin(), line.end(), '\t') - 8;
        break;
      }
    }

    double start_cm = 0.;
    auto tar_it = target_sites.begin();
    unique_haplotype_block block;
    std::vector<std::int8_t> tmp_geno;
    while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples))
    {
      if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
        break;

      block.fill_cm_from_recom(start_cm);

      for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
      {
        while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];

            typed_only_reference_data.compress_variant({ref_it->chrom, ref_it->pos, ref_it->id, ref_it->ref, ref_it->alt, ref_it->err, ref_it->recom, ref_it->cm}, tmp_geno);

            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;

            if (it != tar_it)
              std::swap(*it, *tar_it);

            ++tar_it;
            break;
          }
        }
      }

      block.trim(impute_reg.from(), impute_reg.to());
      if (!block.variants().empty())
        full_reference_data.append_block(block);
    }
  }

  return !input.bad();
}

// Old approach to setting recombination probs.
// saving this to compare with new appraoch.
bool load_reference_haplotypes_old_recom_approach(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data,
  genetic_map_file* map_file)
{
  savvy::reader input(file_path);

  if (input)
  {
    if (!input.reset_bounds(extended_reg, savvy::bounding_point::any))
      return std::cerr << "Error: reference file must be indexed MVCF\n", false;

    bool is_m3vcf_v3 = false;
    for (auto it = input.headers().begin(); !is_m3vcf_v3 && it != input.headers().end(); ++it)
    {
      if (it->first == "subfileformat" && (it->second == "M3VCFv3.0" || it->second == "MVCFv3.0"))
        is_m3vcf_v3 = true;
    }

    if (!is_m3vcf_v3)
      return std::cerr << "Error: reference file must be an MVCF\n", false;

    if (subset_ids.size() && input.subset_samples(subset_ids).empty())
      return std::cerr << "Error: no reference samples overlap subset IDs\n", false;

    savvy::variant var;
    if (!input.read(var))
      return std::cerr << "Notice: no variant records in reference query region (" << extended_reg.chromosome() << ":" << extended_reg.from() << "-" << extended_reg.to() << ")\n", input.bad() ? false : true;

    std::size_t ref_read_cnt = 0;
    float prev_recom = 0.f;
    double prev_cm = 0.;
    double no_recom_prob = 1.;
    auto recom_it = target_sites.end();
    std::vector<std::int8_t> tmp_geno;
    unique_haplotype_block block;
    auto tar_it = target_sites.begin();
    int res;
    while ((res = block.deserialize(input, var)) > 0)
    {
      if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
        break;

      for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
      {
        while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            if (recom_it != target_sites.end())
            {
              if (recom_it->pos != ref_it->pos)
              {
                recom_it->recom = std::max(1. - no_recom_prob, 1e-5);
                no_recom_prob = 1.;
              }
            }

            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];

            typed_only_reference_data.compress_variant({ref_it->chrom, ref_it->pos, ref_it->id, ref_it->ref, ref_it->alt, ref_it->err, ref_it->recom, ref_it->cm}, tmp_geno);

            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;

            if (it != tar_it)
              std::swap(*it, *tar_it);

            no_recom_prob = 1.;
            recom_it = tar_it;
            ++tar_it;
            break;
          }
        }

        if (map_file)
        {
          double cm = map_file->interpolate_centimorgan(ref_it->pos);
          prev_recom = recombination::cm_to_switch_prob(cm - prev_cm);
          prev_cm = cm;
        }
        if (ref_read_cnt)
          no_recom_prob *= 1. - prev_recom; //std::max(prev_recom, 1e-5f); // TODO: cm -> recom
        prev_recom = ref_it->cm; // TODO: cm -> recom
        ++ref_read_cnt;
      }

      block.trim(impute_reg.from(), impute_reg.to());
      if (!block.variants().empty())
        full_reference_data.append_block(block);
    }

    assert(recom_it != target_sites.end() && recom_it + 1 == target_sites.end());
    recom_it->recom = 0.f;

    if (res < 0)
      return false;


    return !input.bad();
  }
  else
  {
    shrinkwrap::gz::istream input_file(file_path);
    std::string line;

    std::uint8_t m3vcf_version = 0;
    const std::string m3vcf_version_line = "##fileformat=M3VCF";
    while (std::getline(input_file, line))
    {
      if (line.substr(0, m3vcf_version_line.size()) == m3vcf_version_line)
      {
        if (line == "##fileformat=M3VCFv2.0")
          m3vcf_version = 2;
        else
          m3vcf_version = 1;
        break;
      }

      if (line.size() < 2 || line[1] != '#')
      {
        std::cerr << "Error: invalid reference file" << std::endl;
        return false;
      }
    }

    std::size_t n_samples = 0;
    while (std::getline(input_file, line))
    {
      if (line.size() < 2 || line[1] != '#')
      {
        n_samples = std::count(line.begin(), line.end(), '\t') - 8;
        break;
      }
    }

    auto tar_it = target_sites.begin();
    unique_haplotype_block block;
    std::vector<std::int8_t> tmp_geno;
    while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples))
    {
      if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
        break;

      for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
      {
        while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];

            typed_only_reference_data.compress_variant({ref_it->chrom, ref_it->pos, ref_it->id, ref_it->ref, ref_it->alt, ref_it->err, ref_it->recom, ref_it->cm}, tmp_geno);

            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;

            if (it != tar_it)
              std::swap(*it, *tar_it);

            ++tar_it;
            break;
          }
        }
      }

      block.trim(impute_reg.from(), impute_reg.to());
      if (!block.variants().empty())
        full_reference_data.append_block(block);
    }
  }

  return !input.bad();
}

std::vector<target_variant> separate_target_only_variants(std::vector<target_variant>& target_sites)
{
  std::vector<target_variant> target_only_sites;
  std::size_t shift_idx = 0;
  for (std::size_t i = 0; i < target_sites.size(); ++i)
  {
    if (!target_sites[i].in_ref)
    {
      target_only_sites.emplace_back();
      std::swap(target_only_sites.back(), target_sites[i]);
    }
    else
    {
      if (i != shift_idx)
      {
        std::swap(target_sites[i], target_sites[shift_idx]);
      }
      ++shift_idx;
    }
  }

  target_sites.resize(shift_idx);
  return target_only_sites;
}

bool load_variant_hmm_params(std::vector<target_variant>& tar_variants, reduced_haplotypes& typed_only_reference_data, float default_error_param, float recom_min, const std::string& map_file_path)
{
  assert(tar_variants.size() == typed_only_reference_data.variant_size());
  auto tar_var_it = tar_variants.begin();
  auto tar_var_end = tar_variants.cend();
  if (tar_var_it == tar_var_end)
    return false;

  if (!map_file_path.empty())
  {
    genetic_map_file map_file(map_file_path, tar_var_it->chrom);
    if (!map_file)
      return std::cerr << "Error: could not open genetic map file\n", false;

    typed_only_reference_data.fill_cm(map_file);
  }

  const auto ref_var_end = typed_only_reference_data.end();
  for (auto ref_var_it = typed_only_reference_data.begin(); ref_var_it != ref_var_end; )
  {
    tar_var_it->err = std::isnan(ref_var_it->err) ? default_error_param : ref_var_it->err;
//    ++tar_var_it; ++ref_var_it; continue;

    auto prev_ref_var_it = ref_var_it++;
    if (ref_var_it == ref_var_end)
      (tar_var_it++)->recom = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
    else
    {
      if (ref_var_it->pos == prev_ref_var_it->pos)
        (tar_var_it++)->recom = 0.f;
      else
      {
        float delta = (ref_var_it->cm - prev_ref_var_it->cm);
        //float recom = std::max(0.01f * delta, recom_min);
        float recom = delta >= 0. ? std::max<float>(recombination::cm_to_switch_prob(delta), recom_min) : recom_min;
        (tar_var_it++)->recom = recom;
      }
    }
  }

  assert(tar_var_it == tar_var_end);

  return true;
}

std::vector<std::vector<std::vector<std::size_t>>> generate_reverse_maps(const reduced_haplotypes& typed_only_reference_data)
{
  std::vector<std::vector<std::vector<std::size_t>>> reverse_maps;

  reverse_maps.reserve(typed_only_reference_data.blocks().size());
  for (auto it = typed_only_reference_data.blocks().begin(); it != typed_only_reference_data.blocks().end(); ++it)
  {
    reverse_maps.emplace_back();
    auto& map = reverse_maps.back();
    for (std::size_t i = 0; i < it->cardinalities().size(); ++i)
    {
      map.emplace_back();
      map.back().reserve(it->cardinalities()[i]);
    }

    for (std::size_t i = 0; i < it->unique_map().size(); ++i)
    { assert(it->unique_map()[i] < map.size());
      map[it->unique_map()[i]].push_back(i);
    }
  }

  return reverse_maps;
}

bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path, const std::string& map_file_path)
{
  std::vector<std::pair<std::string, std::string>> headers;
  std::vector<std::string> ids;

  shrinkwrap::gz::istream input_file(input_path);
  std::string line;

  bool phasing_header_present = false;
  bool contig_header_present = false;
  std::uint8_t m3vcf_version = 0;
  const std::string m3vcf_version_line = "##fileformat=M3VCF";
  const std::string vcf_version_line = "##fileformat=VCF";
  while (std::getline(input_file, line))
  {
    std::size_t equal_pos = line.find('=');
    if (equal_pos != std::string::npos)
    {
      std::string key = line.substr(0, equal_pos);
      std::string val = line.substr(equal_pos + 1);
      key.erase(0, key.find_first_not_of('#'));

      if (key == "fileformat")
      {
        if (val.substr(0, 5) == "M3VCF")
        {
          if (val == "M3VCFv2.0")
            m3vcf_version = 2;
          else
            m3vcf_version = 1;
        }
      }
      else if (key != "INFO" && key != "FORMAT")
      {
        if (!phasing_header_present && key == "phasing")
          phasing_header_present = true;
        else if (!contig_header_present && key == "contig")
          contig_header_present = true;

        headers.emplace_back(std::move(key), std::move(val));
      }
    }
    else
    {

      break;
    }
  }

  if (line.size() < 1 || line[0] != '#')
  {
    std::cerr << "Error: invalid reference file" << std::endl;
    return false;
  }

  headers.insert(headers.begin(), {"subfileformat","MVCFv3.0"});
  headers.insert(headers.begin(), {"fileformat","VCFv4.2"});
  if (!phasing_header_present)
    headers.emplace_back("phasing","full");
  headers.emplace_back("INFO", "<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">");
  headers.emplace_back("INFO", "<ID=AN,Number=1,Type=Float,Description=\"Total number of alleles in called genotypes\">");
  headers.emplace_back("INFO","<ID=REPS,Number=1,Type=Integer,Description=\"Number of distinct haplotypes in block\">");
  headers.emplace_back("INFO","<ID=VARIANTS,Number=1,Type=Integer,Description=\"Number of variants in block\">");
  headers.emplace_back("INFO","<ID=ERR,Number=1,Type=Float,Description=\"Error parameter for HMM\">");
  headers.emplace_back("INFO","<ID=RECOM,Number=1,Type=Float,Description=\"Recombination probability\">");
  if (!map_file_path.empty())
    headers.emplace_back("INFO","<ID=CM,Number=1,Type=Float,Description=\"Centimorgan\">");
  headers.emplace_back("INFO","<ID=END,Number=1,Type=Integer,Description=\"End position of record\">");
  headers.emplace_back("INFO","<ID=UHA,Number=.,Type=Integer,Description=\"Unique haplotype alleles\">");
  headers.emplace_back("FORMAT","<ID=UHM,Number=.,Type=Integer,Description=\"Unique haplotype mapping\">");


  //headers.emplace_back("ALT","<ID=DUP,Description=\"Duplication\">");


  std::size_t tab_cnt = 0;
  std::size_t last_pos = 0;
  std::size_t tab_pos = 0;
  while ((tab_pos = line.find('\t', tab_pos)) != std::string::npos)
  {
    if (tab_cnt >= 9)
    {
      ids.emplace_back(line.substr(last_pos, tab_pos - last_pos));
    }
    last_pos = ++tab_pos;
    ++tab_cnt;
  }

  ids.emplace_back(line.substr(last_pos, tab_pos - last_pos));
  std::size_t n_samples = ids.size();

  unique_haplotype_block block;
  block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples);
  if (!contig_header_present && block.variants().size())
    headers.emplace_back("contig","<ID=" + block.variants()[0].chrom + ">");

  std::string last_3;
  if (output_path.size() >= 3)
    last_3 = output_path.substr(output_path.size() - 3);
  savvy::writer output_file(output_path, last_3 == "bcf" ? savvy::file::format::bcf : savvy::file::format::sav, headers, ids, 6);


  std::size_t block_cnt = 0;

  std::unique_ptr<genetic_map_file> map_file;
  if (!map_file_path.empty() && !block.variants().empty())
  {
    map_file.reset(new genetic_map_file(map_file_path, block.variants()[0].chrom));
    if (!map_file->good())
      return std::cerr << "Error: could not open map file\n", false;
  }

  std::vector<std::int8_t> tmp_geno;
  do
  {
    if (block.variants().empty())
      break;

    if (map_file)
      block.fill_cm(*map_file);

    if (!block.serialize(output_file))
      return false;
    ++block_cnt;
  } while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples));

  return !input_file.bad() && output_file.good();
}

bool compress_reference_panel(const std::string& input_path, const std::string& output_path, const std::string& map_file_path)
{
  savvy::reader input_file(input_path);
  savvy::variant var;
  std::vector<std::int8_t> gts;

  std::vector<std::pair<std::string, std::string>> headers = {
    {"subfileformat","MVCFv3.0"},
    {"fileformat","VCFv4.2"},
    {"phasing","full"},
    {"INFO", "<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">"},
    {"INFO", "<ID=AN,Number=1,Type=Float,Description=\"Total number of alleles in called genotypes\">"},
    {"INFO","<ID=REPS,Number=1,Type=Integer,Description=\"Number of distinct haplotypes in block\">"},
    {"INFO","<ID=VARIANTS,Number=1,Type=Integer,Description=\"Number of variants in block\">"},
    {"INFO","<ID=ERR,Number=1,Type=Float,Description=\"Error parameter for HMM\">"},
    {"INFO","<ID=RECOM,Number=1,Type=Float,Description=\"Recombination probability\">"},
    /*{"INFO","<ID=CM,Number=1,Type=Float,Description=\"Centimorgan\">"},*/
    {"INFO","<ID=END,Number=1,Type=Integer,Description=\"End position of record\">"},
    {"INFO","<ID=UHA,Number=.,Type=Integer,Description=\"Unique haplotype alleles\">"},
    {"FORMAT","<ID=UHM,Number=.,Type=Integer,Description=\"Unique haplotype mapping\">"},
  };

  std::unordered_set<std::string> header_contigs;
  for (auto it = input_file.headers().begin(); it != input_file.headers().end(); ++it)
  {
    if (it->first == "phasing")
    {
      if (it->second == "none" || it->second == "partial")
        return std::cerr << "Error: phased genotypes are required (phasing header status is '" << it->second << "')\n", false;
      continue;
    }

    if (it->first == "contig")
      header_contigs.insert(savvy::parse_header_sub_field(it->second, "ID"));

    if (it->first != "fileformat" && it->first != "subfileformat" && it->first != "INFO" && it->first != "FORMAT")
      headers.emplace_back(it->first, it->second);
  }

  if (!input_file.read(var) || !var.get_format("GT", gts))
    return std::cerr << "Error: could not read GT from first variant record in input file\n", false;

  if (header_contigs.find(var.chrom()) == header_contigs.end())
    headers.emplace_back("contig", "<ID=" + var.chrom() + ">");

  float flt_nan = std::numeric_limits<float>::quiet_NaN();
  unique_haplotype_block block;
  block.compress_variant(reference_site_info(var.chrom(), var.pos(), var.id(), var.ref(), var.alts().size() ? var.alts()[0] : "", flt_nan, flt_nan, std::numeric_limits<double>::quiet_NaN()), gts);
  std::size_t variant_cnt = 1;
  std::size_t block_cnt = 0;

  std::string last_3;
  if (output_path.size() >= 3)
    last_3 = output_path.substr(output_path.size() - 3);
  savvy::writer output_file(output_path, last_3 == "bcf" ? savvy::file::format::bcf : savvy::file::format::sav, headers, input_file.samples(), 6);

  auto comp_ratio = [](const unique_haplotype_block& b)
  {
    return float(b.expanded_haplotype_size() + b.unique_haplotype_size() * b.variant_size()) / float(b.expanded_haplotype_size() * b.variant_size());
  };

  const std::size_t min_block_size = 10;
  const std::size_t max_block_size = 0xFFFF; // max s1r block size minus 1 partition record

  bool flush_block = false;
  while (input_file >> var)
  {
    float old_cr = comp_ratio(block);

    if (!var.get_format("GT", gts))
      return std::cerr << "Error: could not read GT from variant record\n", false;

    bool ret = block.compress_variant(
      reference_site_info(var.chrom(), var.pos(), var.id(), var.ref(), var.alts().size() ? var.alts()[0] : "", flt_nan, flt_nan, std::numeric_limits<double>::quiet_NaN()),
      gts);

    if (!ret)
      return std::cerr << "Error: compressing variant failed\n", false;

    ++variant_cnt;

    std::size_t cnt = block.variant_size();
    if (cnt >= min_block_size)
    {
      float new_cr = comp_ratio(block);
      if (new_cr > old_cr || cnt >= max_block_size)
      {
        // CM INFO field would need to be a double in order to have enough precision
        // if (map_file)
        //  block.fill_cm(*map_file);

        if (!block.serialize(output_file))
          return std::cerr << "Error: serializing block failed\n", false;
        block.clear();
        ++block_cnt;
      }
    }
  }

  // CM INFO field would need to be a double in order to have enough precision
  // if (map_file)
  //  block.fill_cm(*map_file);

  if (block.variants().size())
  {
    if (!block.serialize(output_file))
      return std::cerr << "Error: serializing final block failed\n", false;
    ++block_cnt;
  }

  return !input_file.bad() && output_file.good();
}
