#ifndef MINIMAC4_RECOMBINATION_HPP
#define MINIMAC4_RECOMBINATION_HPP

#include "variant.hpp"

#include <shrinkwrap/gz.hpp>

#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>

class recombination
{
public:
  //static constexpr float recom_min = 1e-5f;
private:
  struct map_file_line
  {
   std::string chrom;
   std::size_t pos = 0;
   float map_value = 0.;
  };
public:
  template <typename VarIter>
  static bool parse_map_file(const std::string& map_file_path, VarIter var_beg, VarIter var_end, float recom_min);
private:
  static bool read_entry(std::istream& ifs, map_file_line& entry, bool new_format);
};

template <typename VarIter>
bool recombination::parse_map_file(const std::string& map_file_path, VarIter var_beg, VarIter var_end, float recom_min)
{
  if (var_beg == var_end)
    return false;

  std::string target_chrom = var_beg->chrom;

  bool new_format = false;
  shrinkwrap::gz::istream ifs(map_file_path);
  while (ifs.peek() == '#') // skip header line
  {
    std::string line;
    std::getline(ifs, line);
    if (std::count(line.begin(), line.end(), '\t') != 2)
      return false;
    new_format = true;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // First align genetic map to target sites
  //std::cerr << "Scanning map file for target chromosome (" << target_chrom << ") ..." << std::endl;
  map_file_line last_entry;
  do
  {
    if (!read_entry(ifs, last_entry, new_format))
      return false;
  } while (last_entry.chrom != target_chrom);

  double basepair_cm = 0.;

  //std::cerr << "First entry: " << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << std::endl;
  auto var_it = var_beg;
  while(var_it != var_end && var_it->pos < last_entry.pos)
  {
    basepair_cm = last_entry.map_value / last_entry.pos;
    var_it->recom = var_it->pos * basepair_cm;
    ++var_it;
  }

  //std::cerr << "Site index less than first entry: " << site_index << std::endl;

  map_file_line entry;
  while (read_entry(ifs, entry, new_format) && entry.chrom == target_chrom)
  {
    assert(entry.pos != last_entry.pos); //TODO: handle gracefully
    assert(entry.pos - last_entry.pos < entry.pos); //TODO: handle gracefully

    basepair_cm = (entry.map_value - last_entry.map_value) / double(entry.pos - last_entry.pos);
    while (var_it != var_end && var_it->pos < entry.pos)
    {
      assert(var_it->pos - last_entry.pos < var_it->pos); //TODO: handle gracefully

      var_it->recom = last_entry.map_value + double(var_it->pos - last_entry.pos) * basepair_cm;
      ++var_it;
      //      if (site_index % 1000 == 0)
      //        std::cerr << "Map file progress: " << site_index << "\t" << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << "\t" << entry.chrom << ":" << entry.pos << " -> " << entry.map_value << std::endl;
    }

    last_entry = entry;
  }

  //std::cerr << "Last entry: " << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << std::endl;
  //std::cerr << "Site index greater than last entry: " << site_index << std::endl;
  while (var_it != var_end)
  {
    var_it->recom = last_entry.map_value + double(var_it->pos - last_entry.pos) * basepair_cm; // Should we assume that basepair_cm is zero?
    ++var_it;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // convert from aligned genetic map to switch probabilities
  for (var_it = var_beg; var_it != var_end; )
  {
    auto prev_var_it = var_it++;
    if (var_it == var_end)
      prev_var_it->recom = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
    else
    {
      if (var_it->pos == prev_var_it->pos)
        prev_var_it->recom = 0.f;
      else
      {
        float delta = (var_it->recom - prev_var_it->recom);
        prev_var_it->recom = std::max(0.01f * delta, recom_min);
        //sites[i].recom = std::min(0.01f, std::max(delta, recom_min));
        //sites[i].recom = (1. - std::exp(-delta/50.))/2.;
      }
    }
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  return true;
}

#endif // MINIMAC4_RECOMBINATION_HPP