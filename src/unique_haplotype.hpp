#ifndef MINIMAC4_UNIQUE_HAPLOTYPE_HPP
#define MINIMAC4_UNIQUE_HAPLOTYPE_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <limits>

struct reference_site
{
  std::string chrom;
  std::uint32_t pos;
  std::string ref;
  std::string alt;
};

class unique_haplotype_block
{
private:
  std::vector<std::size_t> unique_map_;
  std::vector<std::size_t> cardinalities_;
  std::vector<std::vector<std::size_t>> reverse_map_;
  std::vector<reference_site> sites_;
  std::vector<std::vector<std::int8_t>> unique_hap_matrix_;
public:
  class unique_matrix_wrapper
  {
  public:
    unique_matrix_wrapper(const unique_haplotype_block& parent) : parent_(parent) {}
    const std::vector<std::int8_t>& operator[](std::size_t row_idx) const { return parent_.unique_hap_matrix_[row_idx]; }
    std::size_t n_rows() const { return parent_.unique_hap_matrix_.size(); }
    std::size_t n_columns() const { return parent_.unique_hap_matrix_.size() ? parent_.unique_hap_matrix_[0].size() : 0; }
  private:
    const unique_haplotype_block& parent_;
  };

  unique_haplotype_block& operator=(const unique_haplotype_block&) = default;
  unique_haplotype_block& operator=(unique_haplotype_block&&) = default;

  bool append_variant(const reference_site& site_info, const std::vector<std::int8_t>& alleles);

  const std::vector<std::vector<std::int8_t>>& compressed_matrix() const { return unique_hap_matrix_; }
  const std::vector<std::size_t>& unique_map() const { return unique_map_; }
  std::size_t expanded_haplotype_size() const { return unique_map_.size(); }
  std::size_t unique_haplotype_size() const { return unique_hap_matrix_.empty() ? 0 : unique_hap_matrix_[0].size(); }
  std::size_t variant_size() const { return unique_hap_matrix_.size(); }
  const std::vector<std::size_t>& cardinalities() const { return cardinalities_; }
};

class reduced_haplotypes
{
private:
  std::deque<unique_haplotype_block> blocks_;
  std::size_t variant_count_ = 0;
  std::size_t min_block_size_ = 1;
  std::size_t max_block_size_ = std::numeric_limits<std::size_t>::max();
public:
  reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size);
  bool append_variant(const reference_site& site_info, const std::vector<std::int8_t>& alleles);
  float compression_ratio() const;
  const std::deque<unique_haplotype_block>& blocks() const { return blocks_; }
  std::size_t variant_size() const { return variant_count_; }
};

#endif // MINIMAC4_UNIQUE_HAPLOTYPE_HPP