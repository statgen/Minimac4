#ifndef MINIMAC4_RECOMBINATION_HPP
#define MINIMAC4_RECOMBINATION_HPP

#include "variant.hpp"

#include <shrinkwrap/gz.hpp>

#include <string>
#include <vector>
#include <cstdint>

class recombination
{
public:
 static constexpr float recom_min = 1e-5f;
private:
 struct map_file_line
 {
   std::string chrom;
   std::size_t pos = 0;
   float map_value = 0.;
 };
public:
 static bool parse_map_file(const std::string& map_file_path, std::vector<target_variant>& sites);
private:
 static bool read_entry(std::istream& ifs, map_file_line& entry, bool new_format);
};

#endif // MINIMAC4_RECOMBINATION_HPP