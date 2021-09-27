#include "recombination.hpp"

#include <cassert>
#include <algorithm>
#include <cmath>

constexpr float recombination::recom_min;

bool recombination::parse_map_file(std::vector<float>& recom_probs, const std::string& map_file_path, const std::vector<target_variant>& sites)
{
  if (sites.empty())
    return false;

  recom_probs.resize(sites.size());
  std::cerr << "Num sites: " << recom_probs.size() << std::endl;

  std::string target_chrom = sites[0].chrom;

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
  std::cerr << "Scanning map file for target chromosome (" << target_chrom << ") ..." << std::endl;
  map_file_line last_entry;
  do
  {
    if (!read_entry(ifs, last_entry, new_format))
      return false;
  } while (last_entry.chrom != target_chrom);

  double basepair_cm = 0.;

  std::cerr << "First entry: " << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << std::endl;
  std::size_t site_index = 0;
  while(site_index < sites.size() && sites[site_index].pos < last_entry.pos)
  {
    basepair_cm = last_entry.map_value / last_entry.pos;
    recom_probs[site_index] = sites[site_index].pos * basepair_cm;
    ++site_index;
  }

  std::cerr << "Site index less than first entry: " << site_index << std::endl;

  map_file_line entry;
  while (read_entry(ifs, entry, new_format) && entry.chrom == target_chrom)
  {
    assert(entry.pos != last_entry.pos); //TODO: handle gracefully
    assert(entry.pos - last_entry.pos < entry.pos); //TODO: handle gracefully

    basepair_cm = (entry.map_value - last_entry.map_value) / double(entry.pos - last_entry.pos);
    while (site_index < sites.size() && sites[site_index].pos < entry.pos)
    {
      assert(sites[site_index].pos - last_entry.pos < sites[site_index].pos); //TODO: handle gracefully

      recom_probs[site_index] = last_entry.map_value + double(sites[site_index].pos - last_entry.pos) * basepair_cm;
      ++site_index;
      if (site_index % 1000 == 0)
        std::cerr << "Map file progress: " << site_index << "\t" << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << "\t" << entry.chrom << ":" << entry.pos << " -> " << entry.map_value << std::endl;
    }

    last_entry = entry;
  }

  std::cerr << "Last entry: " << last_entry.chrom << ":" << last_entry.pos << " -> " << last_entry.map_value << std::endl;
  std::cerr << "Site index greater than last entry: " << site_index << std::endl;
  while (site_index < sites.size())
  {
    recom_probs[site_index] = last_entry.map_value + double(sites[site_index].pos - last_entry.pos) * basepair_cm; // Should we assume that basepair_cm is zero?
    ++site_index;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // convert from aligned genetic map to switch probabilities
  for (std::size_t i = 0; i < recom_probs.size(); ++i)
  {
    if (i + 1 == recom_probs.size())
      recom_probs[i] = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
    else
    {
      float delta = (recom_probs[i + 1] - recom_probs[i]);
      recom_probs[i] = std::max(0.01f * delta, recom_min);
      //recom_probs[i] = std::min(0.01f, std::max(/*0.01f * */delta, recom_min));
      //recom_probs[i] = (1. - std::exp(-delta/50.))/2.;
    }
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  return true;
}

bool recombination::read_entry(std::istream& ifs, map_file_line& entry, bool new_format)
{
  std::string discard;
  if (new_format)
    ifs >> entry.chrom >> entry.pos >> entry.map_value;
  else
    ifs >> entry.chrom >> discard >> entry.map_value >> entry.pos;

  if (!ifs || entry.chrom.empty())
    return false;
  return true;
}