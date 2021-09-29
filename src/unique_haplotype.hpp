#ifndef MINIMAC4_UNIQUE_HAPLOTYPE_HPP
#define MINIMAC4_UNIQUE_HAPLOTYPE_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <limits>
#include <cassert>

struct reference_site_info
{
  std::string chrom;
  std::uint32_t pos;
  std::string ref;
  std::string alt;

  reference_site_info(std::string _chrom,
    std::uint32_t _pos,
    std::string _ref,
    std::string _alt)
  :
    chrom(std::move(_chrom)),
    pos(_pos),
    ref(std::move(_ref)),
    alt(std::move(_alt))
  {
  }
};

struct reference_variant : public reference_site_info
{
  reference_variant(const std::string& _chrom,
    std::uint32_t _pos,
    const std::string& _ref,
    const std::string& _alt,
    std::size_t _ac,
    std::vector<std::int8_t> _gt)
    :
    reference_site_info(_chrom, _pos, _ref, _alt),
    ac(_ac),
    gt(std::move(_gt))
  {
  }
  std::size_t ac;
  std::vector<std::int8_t> gt;
};

class unique_haplotype_block
{
private:
  std::vector<std::size_t> unique_map_;
  std::vector<std::size_t> cardinalities_;
  std::vector<std::vector<std::size_t>> reverse_map_;
  std::vector<reference_variant> sites_;
  //std::vector<std::vector<std::int8_t>> unique_hap_matrix_;
public:
//  class unique_matrix_wrapper
//  {
//  public:
//    unique_matrix_wrapper(const unique_haplotype_block& parent) : parent_(parent) {}
//    const std::vector<std::int8_t>& operator[](std::size_t row_idx) const { return parent_.unique_hap_matrix_[row_idx]; }
//    std::size_t n_rows() const { return parent_.unique_hap_matrix_.size(); }
//    std::size_t n_columns() const { return parent_.unique_hap_matrix_.size() ? parent_.unique_hap_matrix_[0].size() : 0; }
//  private:
//    const unique_haplotype_block& parent_;
//  };

  unique_haplotype_block& operator=(const unique_haplotype_block&) = default;
  unique_haplotype_block& operator=(unique_haplotype_block&&) = default;

  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles);

  //const std::vector<std::vector<std::int8_t>>& compressed_matrix() const { return unique_hap_matrix_; }
  const std::vector<reference_variant>& variants() const { return sites_; }
  const std::vector<std::size_t>& unique_map() const { return unique_map_; }
  std::size_t expanded_haplotype_size() const { return unique_map_.size(); }
  std::size_t unique_haplotype_size() const { return sites_.empty() ? 0 : sites_[0].gt.size(); }
  std::size_t variant_size() const { return sites_.size(); }
  const std::vector<std::size_t>& cardinalities() const { return cardinalities_; }
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
        ++block_idx_;
      return *this;
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
    std::size_t global_idx() const { return parent_->block_offsets_[block_idx_] + parent_->blocks_[block_idx_].variant_size(); }
    const std::vector<std::size_t>& unique_map() const { parent_->blocks_[block_idx_].unique_map(); }
    const std::vector<std::size_t>& cardinalities() const { return parent_->blocks_[block_idx_].cardinalities(); }
  };

  iterator begin() const { return iterator(*this, 0, 0); }
  iterator end() const { return iterator(*this, blocks_.size() - 1, blocks_.size() ? blocks_.back().variant_size() : 0); }

  reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size);
  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles);
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