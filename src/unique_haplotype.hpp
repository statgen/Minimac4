#ifndef MINIMAC4_UNIQUE_HAPLOTYPE_HPP
#define MINIMAC4_UNIQUE_HAPLOTYPE_HPP

#include "variant.hpp"

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <limits>
#include <cassert>

class unique_haplotype_block
{
private:
  std::vector<std::int64_t> unique_map_;
  std::vector<std::size_t> cardinalities_;
  std::vector<reference_variant> variants_;
public:
  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles);

  const std::vector<reference_variant>& variants() const { return variants_; }
  const std::vector<std::int64_t>& unique_map() const { return unique_map_; }
  std::size_t expanded_haplotype_size() const { return unique_map_.size(); }
  std::size_t unique_haplotype_size() const { return variants_.empty() ? 0 : variants_[0].gt.size(); }
  std::size_t variant_size() const { return variants_.size(); }
  const std::vector<std::size_t>& cardinalities() const { return cardinalities_; }

  void clear();
  void trim(std::size_t min_pos, std::size_t max_pos);
  void pop_variant();
  bool deserialize(std::istream& is, int m3vcf_version, std::size_t n_haplotypes);
  int deserialize(savvy::reader& input_file, savvy::variant& var); // return <0 is error, 0 is EOF, >0 good
  bool serialize(savvy::writer& output_file);
};

class reduced_haplotypes
{
private:
  std::vector<std::size_t> block_offsets_;
  std::deque<unique_haplotype_block> blocks_;
  std::size_t variant_count_ = 0;
  std::size_t min_block_size_ = 1;
  std::size_t max_block_size_ = std::numeric_limits<std::size_t>::max();
public:
  class iterator
  {
  private:
    const reduced_haplotypes* parent_;
    std::size_t block_idx_;
    std::size_t variant_idx_;
  public:
    iterator(const reduced_haplotypes& parent, std::size_t block_idx, std::size_t variant_idx) :
      parent_(&parent),
      block_idx_(block_idx),
      variant_idx_(variant_idx)
    {
    }

    iterator& operator++()
    {
      ++variant_idx_;
      assert(block_idx_ < parent_->blocks_.size());
      if (variant_idx_ >= parent_->blocks_[block_idx_].variant_size())
      {
        ++block_idx_;
        variant_idx_ = 0;
      }
      return *this;
    }

    iterator operator++(int)
    {
      iterator ret(*this);
      ++(*this);
      return ret;
    }

    iterator& operator--()
    {
      if (variant_idx_ == 0)
      {
        --block_idx_;
        if (block_idx_ < parent_->blocks_.size())
        {
          (*this) = iterator(*parent_, block_idx_, parent_->blocks_[block_idx_].variant_size());
          --(*this);
        }
      }
      else
      {
        --variant_idx_;
      }

      return *this;
    }

    iterator operator--(int)
    {
      iterator ret(*this);
      --(*this);
      return ret;
    }

    const reference_variant& operator*() const
    {
      return parent_->blocks_[block_idx_].variants()[variant_idx_];
    }

    const reference_variant* operator->() const
    {
      return &(parent_->blocks_[block_idx_].variants()[variant_idx_]);
    }

    std::size_t block_idx() const { return block_idx_; }
    std::size_t block_local_idx() const { return variant_idx_; }
    std::size_t global_idx() const { return parent_->block_offsets_[block_idx_] + variant_idx_; }
    const std::vector<std::int64_t>& unique_map() const { return parent_->blocks_[block_idx_].unique_map(); }
    const std::vector<std::size_t>& cardinalities() const { return parent_->blocks_[block_idx_].cardinalities(); }
  };

  iterator begin() const { return iterator(*this, 0, 0); }
  iterator end() const { return iterator(*this, blocks_.size(), 0); }

  reduced_haplotypes() {}
  reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size);
  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles, bool flush_block = false);
  void append_block(const unique_haplotype_block& block);
  float compression_ratio() const;
  const std::deque<unique_haplotype_block>& blocks() const { return blocks_; }
  std::size_t variant_size() const { return variant_count_; }
};

inline bool operator==(const reduced_haplotypes::iterator& lhs, const reduced_haplotypes::iterator& rhs)
{
  return lhs.block_idx() == rhs.block_idx() && lhs.block_local_idx() == rhs.block_local_idx();
}

inline bool operator!=(const reduced_haplotypes::iterator& lhs, const reduced_haplotypes::iterator& rhs)
{
  return !(lhs == rhs);
}



#endif // MINIMAC4_UNIQUE_HAPLOTYPE_HPP