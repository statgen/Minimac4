#include "input_prep.hpp"
#include "recombination.hpp"


bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids)
{
  savvy::reader input(file_path);
  sample_ids = input.samples();
  input.reset_bounds(reg);
  savvy::variant var;
  std::vector<std::int8_t> tmp_geno;
  while (input >> var)
  {
    var.get_format("GT", tmp_geno);
    for (std::size_t i = 0; i < var.alts().size(); ++i)
    {
      std::size_t allele_idx = i + 1;
      target_sites.push_back({var.chromosome(), var.position(), var.ref(), var.alts()[i], true, false, std::numeric_limits<float>::quiet_NaN(), 0.00999, recombination::recom_min, 0, {}});
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

bool load_reference_haplotypes(const std::string& file_path, const savvy::genomic_region& extended_reg, const savvy::genomic_region& impute_reg, std::vector<target_variant>& target_sites, reduced_haplotypes& typed_only_reference_data, reduced_haplotypes* full_reference_data)
{
  savvy::reader input(file_path);

  if (input)
  {
    input.reset_bounds(extended_reg, savvy::bounding_point::any);
    savvy::variant var;
    std::vector<std::int8_t> tmp_geno;

    if (full_reference_data)
    {
      bool is_m3vcf_v3 = false;
      for (auto it = input.headers().begin(); !is_m3vcf_v3 && it != input.headers().end(); ++it)
      {
        if (it->first == "subfileformat" && (it->second == "M3VCFv3.0" || it->second == "MVCFv3.0"))
          is_m3vcf_v3 = true;
      }

      if (!is_m3vcf_v3)
        return std::cerr << "Error: reference file must be an M3VCF\n", false;

      savvy::variant var;
      if (!input.read(var))
        return std::cerr << "Error: no variant records in reference query region\n", false;

      auto tar_it = target_sites.begin();
      unique_haplotype_block block;
      std::size_t ref_cnt = 0;
      int res;
      while ((res = block.deserialize(input, var)) > 0)
      {
        if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
          break;

        for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
        {
          while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          {
            if (ref_cnt)
            {
              tar_it->ref_cnt = ref_cnt;
              ref_cnt = 0;
            }
            ++tar_it;
          }

          for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
          {
            if (it->ref == ref_it->ref && it->alt == ref_it->alt)
            {
              tmp_geno.resize(block.unique_map().size());
              for (std::size_t i = 0; i < tmp_geno.size(); ++i)
                tmp_geno[i] = ref_it->gt[block.unique_map()[i]];
              typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
              it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
              it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
              it->in_ref = true;
              if (it != tar_it)
                std::swap(*it, *tar_it);
              if (ref_cnt)
              {
                tar_it->ref_cnt = ref_cnt;
                ref_cnt = 0;
              }
              ++tar_it;
              break;
            }
          }

          if (ref_it->pos >= extended_reg.from())
            ++ref_cnt;
        }

        if (full_reference_data)
        {
          // TODO: remove redundant variants
          block.trim(impute_reg.from(), impute_reg.to());
          if (!block.variants().empty())
            full_reference_data->append_block(block);
        }
      }

      if (res < 0)
        return false;
    }
    else
    {
      auto tar_it = target_sites.begin();
      while (input >> var)
      {
        while (tar_it != target_sites.end() && tar_it->pos < var.position())
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == var.position(); ++it)
        {
          if (it->ref == var.ref() && it->alt == (var.alts().size() ? var.alts()[0] : ""))
          {
            var.get_format("GT", tmp_geno);
            typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
            // freq.push_back(std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size());
            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size(); // TODO; remove
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;
            if (it != tar_it)
              std::swap(*it, *tar_it);
            ++tar_it;
            break;
          }
        }
      }
    }

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

    std::size_t ref_cnt = 0;

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
        {
          if (ref_cnt)
          {
            tar_it->ref_cnt = ref_cnt;
            ref_cnt = 0;
          }
          ++tar_it;
        }

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];
            typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;
            if (it != tar_it)
              std::swap(*it, *tar_it);
            if (ref_cnt)
            {
              tar_it->ref_cnt = ref_cnt;
              ref_cnt = 0;
            }
            ++tar_it;
            break;
          }
        }

        if (ref_it->pos >= extended_reg.from())
          ++ref_cnt;
      }

      if (full_reference_data)
      {
        // TODO: remove redundant variants
        block.trim(impute_reg.from(), impute_reg.to());
        if (!block.variants().empty())
          full_reference_data->append_block(block);
      }
    }

    if (ref_cnt && tar_it != target_sites.end())
    {
      tar_it->ref_cnt = ref_cnt;
      ref_cnt = 0;
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

bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path)
{
  std::vector<std::pair<std::string, std::string>> headers;
  std::vector<std::string> ids;

  shrinkwrap::gz::istream input_file(input_path);
  std::string line;

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
      else
      {
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

  headers.insert(headers.begin(), {"subfileformat","M3VCFv3.0"});
  headers.insert(headers.begin(), {"fileformat","VCFv4.2"});
  headers.emplace_back("INFO","<ID=REPS,Number=1,Type=Integer,Description=\"Number of distinct haplotypes in block\">");
  headers.emplace_back("INFO","<ID=VARIANTS,Number=1,Type=Integer,Description=\"Number of variants in block\">");
  headers.emplace_back("INFO","<ID=END,Number=1,Type=Integer,Description=\"End position of record\">");
  headers.emplace_back("INFO","<ID=UHA,Number=.,Type=Integer,Description=\"Unique haplotype alleles\">");
  headers.emplace_back("FORMAT","<ID=UHM,Number=.,Type=Integer,Description=\"Unique haplotype mapping\">");

  //headers.emplace_back("ALT","<ID=DUP,Description=\"Duplication\">");


  std::size_t tab_cnt = 0;
  std::size_t last_pos = 0;
  std::size_t tab_pos;
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

  savvy::writer output_file(output_path, savvy::file::format::bcf, headers, ids, 6);


  std::size_t block_cnt = 0;
  unique_haplotype_block block;
  std::vector<std::int8_t> tmp_geno;
  while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples))
  {
    if (block.variants().empty())
      break;

    block.serialize(output_file);
    ++block_cnt;
  }

  return !input_file.bad() && output_file.good();
}