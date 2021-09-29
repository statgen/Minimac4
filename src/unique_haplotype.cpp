
#include "unique_haplotype.hpp"

#include <iostream>
#include <numeric>
#include <cassert>

bool unique_haplotype_block::compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles)
{
  if (alleles.empty())
    return false;

  sites_.emplace_back(site_info.chrom, site_info.pos, site_info.ref, site_info.alt, 0, std::vector<std::int8_t>());

  if (sites_.size() == 1)
  {
    unique_map_.resize(alleles.size());
    sites_[0].gt.push_back(alleles[0]);
    cardinalities_.push_back(1);
    unique_map_[0] = 0;

    for (std::size_t i = 1; i < alleles.size(); ++i)
    {
      std::size_t j = 0;
      for (; j < sites_[0].gt.size(); ++j)
      {
        if (sites_[0].gt[j] == alleles[i])
          break;
      }

      if (j == sites_[0].gt.size())
      {
        sites_[0].gt.push_back(alleles[i]);
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

    std::size_t original_hap_cnt = sites_[0].gt.size();
    sites_.back().gt.resize(original_hap_cnt, std::int8_t(-1));

    for (std::size_t i = 0; i < alleles.size(); ++i)
    {
      std::size_t original_hap_idx = unique_map_[i];
      if (sites_.back().gt[original_hap_idx] != alleles[i])
      {
        if (sites_.back().gt[original_hap_idx] == -1)
        {
          // this is the first haplotype mapped to this column, so set allele for new variant
          sites_.back().gt[original_hap_idx] = alleles[i];
        }
        else
        {
          // this haplotype no longer matches its column
          --cardinalities_[original_hap_idx];
          std::size_t j = original_hap_cnt;
          std::size_t j_end = sites_[0].gt.size();
          const std::size_t k_end = sites_.size() - 1;

          for ( ; j < j_end; ++j)
          {
            // check if this haplotype matches any newly created column
            std::size_t k = 0;
            for ( ; k < k_end; ++k)
            {
              if (sites_[k].gt[j] != sites_[k].gt[original_hap_idx])
                break;
            }

            if (k == k_end && sites_[k].gt[j] == alleles[i])
              break;
          }

          if (j == j_end)
          {
            // does not match a newly created column, so insert new column
            for (std::size_t k = 0; k < k_end; ++k)
              sites_[k].gt.push_back(sites_[k].gt[original_hap_idx]);
            sites_[k_end].gt.push_back(alleles[i]);
            cardinalities_.push_back(0);
          }

          unique_map_[i] = j;
          ++cardinalities_[j];
        }
      }
    }
  }

  assert(cardinalities_.size() == sites_.back().gt.size());
  for (std::size_t i = 0; i < cardinalities_.size(); ++i)
  {
    if (sites_.back().gt[i])
      sites_.back().ac += cardinalities_[i];
  }

  assert(std::accumulate(cardinalities_.begin(), cardinalities_.end(), std::size_t(0)) == unique_map_.size());

  return true;
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