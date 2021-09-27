
#include "unique_haplotype.hpp"

#include <iostream>
#include <numeric>
#include <cassert>

bool unique_haplotype_block::append_variant(const reference_site& site, const std::vector<std::int8_t>& allele_vec)
{
  if (allele_vec.empty())
    return false;

  if (sites_.size() == 0)
  {
    unique_map_.resize(allele_vec.size());
    unique_hap_matrix_.emplace_back(1, allele_vec[0]);
    cardinalities_.push_back(1);
    unique_map_[0] = 0;

    for (std::size_t i = 1; i < allele_vec.size(); ++i)
    {
      std::size_t j = 0;
      for (; j < unique_hap_matrix_[0].size(); ++j)
      {
        if (unique_hap_matrix_[0][j] == allele_vec[i])
          break;
      }

      if (j == unique_hap_matrix_[0].size())
      {
        unique_hap_matrix_[0].emplace_back(allele_vec[i]);
        cardinalities_.push_back(0);
      }
      unique_map_[i] = j;
      ++cardinalities_[j];
    }
  }
  else
  {
    if (allele_vec.size() != unique_map_.size())
      return false;

    std::size_t original_hap_cnt = unique_hap_matrix_[0].size();
    unique_hap_matrix_.emplace_back(original_hap_cnt, std::int8_t(-1));

    for (std::size_t i = 0; i < allele_vec.size(); ++i)
    {
      std::size_t original_hap_idx = unique_map_[i];
      if (unique_hap_matrix_.back()[original_hap_idx] != allele_vec[i])
      {
        if (unique_hap_matrix_.back()[original_hap_idx] == -1)
        {
          // this is the first haplotype mapped to this column, so set allele for new variant
          unique_hap_matrix_.back()[original_hap_idx] = allele_vec[i];
        }
        else
        {
          // this haplotype no longer matches its column
          --cardinalities_[original_hap_idx];
          std::size_t j = original_hap_cnt;
          std::size_t j_end = unique_hap_matrix_[0].size();
          const std::size_t k_end = unique_hap_matrix_.size() - 1;

          for ( ; j < j_end; ++j)
          {
            // check if this haplotype matches any newly created column
            std::size_t k = 0;
            for ( ; k < k_end; ++k)
            {
              if (unique_hap_matrix_[k][j] != unique_hap_matrix_[k][original_hap_idx])
                break;
            }

            if (k == k_end && unique_hap_matrix_[k][j] == allele_vec[i])
              break;
          }

          if (j == j_end)
          {
            // does not match a newly created column, so insert new column
            for (std::size_t k = 0; k < k_end; ++k)
              unique_hap_matrix_[k].push_back(unique_hap_matrix_[k][original_hap_idx]);
            unique_hap_matrix_[k_end].push_back(allele_vec[i]);
            cardinalities_.push_back(0);
          }

          unique_map_[i] = j;
          ++cardinalities_[j];
        }
      }
    }
  }

  sites_.emplace_back(site);

  assert(std::accumulate(cardinalities_.begin(), cardinalities_.end(), std::size_t(0)) == unique_map_.size());

  return true;
}

reduced_haplotypes::reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size)
{
  blocks_.emplace_back();
  min_block_size_ = std::max(std::size_t(1), min_block_size);
  max_block_size_ = std::max(std::size_t(1), max_block_size);
}

bool reduced_haplotypes::append_variant(const reference_site& site, const std::vector<std::int8_t>& alleles)
{
  auto comp_ratio = [](const unique_haplotype_block& b)
    {
      return float(b.expanded_haplotype_size() + b.unique_haplotype_size() * b.variant_size()) / float(b.expanded_haplotype_size() * b.variant_size());
    };

  float old_cr = comp_ratio(blocks_.back());
  bool ret = blocks_.back().append_variant(site, alleles);
  if (ret)
    ++variant_count_;

  std::size_t cnt = blocks_.back().variant_size();
  if (cnt >= min_block_size_)
  {
    float new_cr = comp_ratio(blocks_.back());
    if (new_cr > old_cr)
    {
      blocks_.emplace_back();
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