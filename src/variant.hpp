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
  std::vector<std::int8_t> gt;
};

#endif // MINIMAC4_VARIANT_HPP