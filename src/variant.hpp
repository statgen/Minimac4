#ifndef MINIMAC4_VARIANT_HPP
#define MINIMAC4_VARIANT_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

struct target_variant
{
  std::string chrom;
  std::uint32_t pos;
  std::string id;
  std::string ref;
  std::string alt;
  bool in_tar; // site exists in target file
  bool in_ref; // site exists in reference file
  float af;
  float err;
  float recom;
  std::vector<std::int8_t> gt;
};

struct reference_site_info
{
  std::string chrom;
  std::uint32_t pos = 0;
  std::string id;
  std::string ref;
  std::string alt;
  float err = std::numeric_limits<float>::quiet_NaN();
  float recom = std::numeric_limits<float>::quiet_NaN();
  double cm = std::numeric_limits<double>::quiet_NaN();
  //float recom = std::numeric_limits<float>::quiet_NaN();

  reference_site_info() {}

  reference_site_info(std::string _chrom,
    std::uint32_t _pos,
    std::string _id,
    std::string _ref,
    std::string _alt,
    float _err,
    float _recom,
    double _cm)
    :
    chrom(std::move(_chrom)),
    pos(_pos),
    id(_id),
    ref(std::move(_ref)),
    alt(std::move(_alt)),
    err(_err),
    recom(_recom),
    cm(_cm)
  {
  }
};

struct reference_variant : public reference_site_info
{
  reference_variant() {}

  reference_variant(const std::string& _chrom,
    std::uint32_t _pos,
    const std::string& _id,
    const std::string& _ref,
    const std::string& _alt,
    float _err,
    float _recom,
    double _cm,
    std::size_t _ac,
    std::vector<std::int8_t> _gt)
    :
    reference_site_info(_chrom, _pos, _id, _ref, _alt, _err, _recom, _cm),
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
    const std::string& _id,
    const std::string& _ref,
    const std::string& _alt,
    float _err,
    float _recom,
    double _cm,
    std::size_t _ac,
    const std::size_t* off_it, const std::size_t* off_it_end)
    :
    reference_site_info(_chrom, _pos, _id, _ref, _alt, _err, _recom, _cm),
    ac(_ac),
    alt_allele_offsets(off_it, off_it_end)
  {
  }
};

#endif // MINIMAC4_VARIANT_HPP