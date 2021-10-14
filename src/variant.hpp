#ifndef MINIMAC4_VARIANT_HPP
#define MINIMAC4_VARIANT_HPP

#include <string>
#include <vector>
#include <cstdint>

struct target_variant
{
  std::string chrom;
  std::uint32_t pos;
  std::string ref;
  std::string alt;
  bool in_tar; // site exists in target file
  bool in_ref; // site exists in reference file
  float af;
  float err;
  float recom;
  std::size_t ref_cnt;
  std::vector<std::int8_t> gt;
};

struct reference_site_info
{
  std::string chrom;
  std::uint32_t pos = 0;
  std::string ref;
  std::string alt;
  //float recom = 0.f;

  reference_site_info() {}

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
  reference_variant() {}

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
  std::size_t ac = 0;
  std::vector<std::int8_t> gt;
};

struct sparse_ref_variant : public reference_site_info
{
  std::size_t ac;
  std::vector<std::uint32_t> alt_allele_offsets;
  sparse_ref_variant(const std::string& _chrom,
    std::uint32_t _pos,
    const std::string& _ref,
    const std::string& _alt,
    std::size_t _ac,
    const std::size_t* off_it, const std::size_t* off_it_end)
    :
    reference_site_info(_chrom, _pos, _ref, _alt),
    ac(_ac),
    alt_allele_offsets(off_it, off_it_end)
  {
  }
};

#endif // MINIMAC4_VARIANT_HPP