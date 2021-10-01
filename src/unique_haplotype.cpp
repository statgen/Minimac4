
#include "unique_haplotype.hpp"

#include <iostream>
#include <numeric>
#include <cassert>
#include <cstring>
#include <algorithm>

bool unique_haplotype_block::compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles)
{
  if (alleles.empty())
    return false;

  variants_.emplace_back(site_info.chrom, site_info.pos, site_info.ref, site_info.alt, 0, std::vector<std::int8_t>());

  if (variants_.size() == 1)
  {
    unique_map_.resize(alleles.size());
    variants_[0].gt.push_back(alleles[0]);
    cardinalities_.push_back(1);
    unique_map_[0] = 0;

    for (std::size_t i = 1; i < alleles.size(); ++i)
    {
      std::size_t j = 0;
      for (; j < variants_[0].gt.size(); ++j)
      {
        if (variants_[0].gt[j] == alleles[i])
          break;
      }

      if (j == variants_[0].gt.size())
      {
        variants_[0].gt.push_back(alleles[i]);
        cardinalities_.push_back(0);
      }
      unique_map_[i] = j;
      ++cardinalities_[j];
    }
  }
  else
  {
    if (alleles.size() != unique_map_.size())
      return false;

    std::size_t original_hap_cnt = variants_[0].gt.size();
    variants_.back().gt.resize(original_hap_cnt, std::int8_t(-1));

    for (std::size_t i = 0; i < alleles.size(); ++i)
    {
      std::size_t original_hap_idx = unique_map_[i];
      if (variants_.back().gt[original_hap_idx] != alleles[i])
      {
        if (variants_.back().gt[original_hap_idx] == -1)
        {
          // this is the first haplotype mapped to this column, so set allele for new variant
          variants_.back().gt[original_hap_idx] = alleles[i];
        }
        else
        {
          // this haplotype no longer matches its column
          --cardinalities_[original_hap_idx];
          std::size_t j = original_hap_cnt;
          std::size_t j_end = variants_[0].gt.size();
          const std::size_t k_end = variants_.size() - 1;

          for ( ; j < j_end; ++j)
          {
            // check if this haplotype matches any newly created column
            std::size_t k = 0;
            for ( ; k < k_end; ++k)
            {
              if (variants_[k].gt[j] != variants_[k].gt[original_hap_idx])
                break;
            }

            if (k == k_end && variants_[k].gt[j] == alleles[i])
              break;
          }

          if (j == j_end)
          {
            // does not match a newly created column, so insert new column
            for (std::size_t k = 0; k < k_end; ++k)
              variants_[k].gt.push_back(variants_[k].gt[original_hap_idx]);
            variants_[k_end].gt.push_back(alleles[i]);
            cardinalities_.push_back(0);
          }

          unique_map_[i] = j;
          ++cardinalities_[j];
        }
      }
    }
  }

  assert(cardinalities_.size() == variants_.back().gt.size());
  for (std::size_t i = 0; i < cardinalities_.size(); ++i)
  {
    if (variants_.back().gt[i])
      variants_.back().ac += cardinalities_[i];
  }

  assert(std::accumulate(cardinalities_.begin(), cardinalities_.end(), std::size_t(0)) == unique_map_.size());

  return true;
}

void unique_haplotype_block::clear()
{
  variants_.clear();
  unique_map_.clear();
  cardinalities_.clear();
}

void unique_haplotype_block::trim(std::size_t min_pos, std::size_t max_pos)
{
  if (variants_.size())
  {
    if (variants_.front().pos > max_pos || variants_.back().pos < min_pos)
    {
      clear();
    }
    else
    {
      std::size_t before_size = variants_.size();
      auto it = --variants_.end();
      for ( ; it != variants_.begin() && it->pos > max_pos; --it) { }
      variants_.erase(it + 1, variants_.end());

      it = variants_.begin();
      for ( ; it != variants_.end() && it->pos < min_pos; ++it) { }
      variants_.erase(variants_.begin(), it);
      std::size_t after_size = variant_size();
      if (before_size != after_size)
      {
        auto a = 0;
      }
    }
  }
}

void unique_haplotype_block::pop_variant()
{
  variants_.pop_back();
}

bool unique_haplotype_block::deserialize(std::istream& is, std::uint8_t m3vcf_version, std::size_t n_haplotypes)
{
  clear();
  unique_map_.reserve(n_haplotypes);

  auto split_string_to_vector = [](const char* in, char delim)
  {
    std::vector<std::string> ret;
    const char* d = nullptr;
    std::string token;
    const char* s = in;
    const char*const e = in + std::strlen(in);
    while ((d = std::find(s, e,  delim)) != e)
    {
      ret.emplace_back(std::string(s, d));
      s = d ? d + 1 : d;
    }
    ret.emplace_back(std::string(s,d));
    return ret;
  };

  std::string line;
  std::getline(is, line);

  std::size_t n_variants = 0;
  std::size_t n_reps = 0;
  std::size_t col_idx = 0;
  std::size_t last_tab_pos = 0;
  std::size_t tab_pos = 0;
  while ((tab_pos = line.find('\t', tab_pos)) != std::string::npos)
  {
    if (col_idx == 7) //INFO
    {
      std::vector<std::string> info_fields = split_string_to_vector(line.substr(last_tab_pos, tab_pos - last_tab_pos).c_str(), ';');
      for (auto it = info_fields.begin(); it != info_fields.end(); ++it)
      {
        if (it->compare(0, 9, "VARIANTS=") == 0)
          n_variants = std::atoll(it->c_str() + 9);
        else if (it->compare(0, 5, "REPS=") == 0)
          n_reps = std::atoll(it->c_str() + 5);
      }
    }
    else if (col_idx >= 9)
    {
      char* p = nullptr;
      unique_map_.push_back(std::strtoll(line.c_str() + last_tab_pos, &p, 10));
      if (m3vcf_version == 2)
      {
        if (*p != '|')
        {
          clear();
          std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
          return false;
        }
        unique_map_.push_back(std::strtoll(++p, &p, 10));
      }
    }

    last_tab_pos = ++tab_pos;
    ++col_idx;
  }

  {
    char* p = nullptr;
    unique_map_.push_back(std::strtoll(line.c_str() + last_tab_pos, &p, 10));
    if (m3vcf_version == 2)
    {
      if (*p != '|')
      {
        clear();
        std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
        return false;
      }
      unique_map_.push_back(std::strtoll(++p, &p, 10));
    }
  }

  if (unique_map_.size() != n_haplotypes)
  {
    clear();
    std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
    return false;
  }

  cardinalities_.resize(n_reps);
  for (auto it = unique_map_.begin(); it != unique_map_.end(); ++it)
    ++cardinalities_[*it];

  variants_.resize(n_variants);
  for (std::size_t i = 0; i < variants_.size(); ++i)
  {
    if (!std::getline(is, line))
    {
      clear();
      std::cerr << "Error: truncated m3vcf v" << m3vcf_version << " file\n";
      return false;
    }

    std::vector<std::string> cols = split_string_to_vector(line.c_str(), '\t');
    if (cols.size() != 9)
    {
      clear();
      std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
      return false;
    }

    variants_[i].chrom = cols[0];
    variants_[i].pos = std::atoll(cols[1].c_str());
    variants_[i].ref = cols[3];
    variants_[i].alt = cols[4];

    if (m3vcf_version == 2)
    {
      variants_[i].gt.resize(n_reps);
      const char* p = cols[8].data();
      char* p_end;
      std::size_t prev_offset = 0;
      do
      {
        std::size_t offset = std::strtoll(p, &p_end, 10);
        if (offset >= n_reps)
        {
          clear();
          std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
          return false;
        }
        std::size_t uniq_idx = prev_offset + offset;
        variants_[i].gt[uniq_idx] = 1;
        variants_[i].ac += cardinalities_[uniq_idx];
        prev_offset = uniq_idx;
        p = p_end + 1;
      } while (*p_end);
    }
    else
    {
      variants_[i].gt.reserve(n_reps);
      const char* p = cols[8].data();
      char* p_end;
      do
      {
        std::uint8_t val = std::strtoll(p, &p_end, 10);
        variants_[i].gt.push_back(val);
        p = p_end + 1;
      } while (*p_end);
      variants_[i].ac = std::inner_product(variants_[i].gt.begin(), variants_[i].gt.end(), cardinalities_.begin(), 0ull);

      if (variants_[i].gt.size() != n_reps)
      {
        clear();
        std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
        return false;
      }
    }
  }

  return is.good();
}

reduced_haplotypes::reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size)
{
  blocks_.emplace_back();
  block_offsets_.emplace_back(0);
  min_block_size_ = std::max(std::size_t(1), min_block_size);
  max_block_size_ = std::max(std::size_t(1), max_block_size);
}

bool reduced_haplotypes::compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles)
{
  auto comp_ratio = [](const unique_haplotype_block& b)
    {
      return float(b.expanded_haplotype_size() + b.unique_haplotype_size() * b.variant_size()) / float(b.expanded_haplotype_size() * b.variant_size());
    };

  float old_cr = comp_ratio(blocks_.back());
  bool ret = blocks_.back().compress_variant(site_info, alleles);
  if (ret)
    ++variant_count_;

  std::size_t cnt = blocks_.back().variant_size();
  if (cnt >= min_block_size_)
  {
    float new_cr = comp_ratio(blocks_.back());
    if (new_cr > old_cr)
    {
      blocks_.emplace_back();
      block_offsets_.push_back(block_offsets_.back() + cnt);
      //std::cerr << "compression old/new: " << old_cr << " / " << new_cr << std::endl;
      return ret;
    }
  }

  //std::cerr << "compression old: " << old_cr << std::endl;

  return ret;
}

bool reduced_haplotypes::append_block(const unique_haplotype_block& block)
{
  assert(block.variants().size());
  if (!blocks_.empty())
  {
    const auto& last_var_prev_block = blocks_.back().variants().back();
    const auto& first_var_new_block = block.variants().front();
    if (last_var_prev_block.pos == first_var_new_block.pos && last_var_prev_block.alt == first_var_new_block.alt)
    {
      blocks_.back().pop_variant();
      --variant_count_;
    }
  }

  block_offsets_.push_back(variant_count_);
  variant_count_ += block.variants().size();
  blocks_.push_back(block);
}

float reduced_haplotypes::compression_ratio() const
{
  float num = 0.f, denom = 0.f;
  for (auto it = blocks_.begin(); it != blocks_.end(); ++it)
  {
    num += it->expanded_haplotype_size() + it->unique_haplotype_size() * it->variant_size();
    denom += it->expanded_haplotype_size() * it->variant_size();
  }

  return num / denom;
}