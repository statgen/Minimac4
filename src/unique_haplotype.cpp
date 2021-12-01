
#include "unique_haplotype.hpp"

#include <iostream>
#include <numeric>
#include <cassert>
#include <cstring>
#include <algorithm>

bool unique_haplotype_block::compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles)
{
  if (alleles.empty())
    return false;

  variants_.emplace_back(site_info.chrom, site_info.pos, site_info.ref, site_info.alt, site_info.err, site_info.cm, 0, std::vector<std::int8_t>());

  if (variants_.size() == 1)
  {
    unique_map_.resize(alleles.size());
    variants_[0].gt.push_back(alleles[0]);
    cardinalities_.push_back(1);
    unique_map_[0] = 0;

    for (std::size_t i = 1; i < alleles.size(); ++i)
    {
      std::size_t j = 0;
      for (; j < variants_[0].gt.size(); ++j)
      {
        if (variants_[0].gt[j] == alleles[i])
          break;
      }

      if (j == variants_[0].gt.size())
      {
        variants_[0].gt.push_back(alleles[i]);
        cardinalities_.push_back(0);
      }
      unique_map_[i] = j;
      ++cardinalities_[j];
    }
  }
  else
  {
    if (alleles.size() != unique_map_.size())
      return false;

    std::size_t original_hap_cnt = variants_[0].gt.size();
    variants_.back().gt.resize(original_hap_cnt, std::int8_t(-1));

    for (std::size_t i = 0; i < alleles.size(); ++i)
    {
      std::size_t original_hap_idx = unique_map_[i];
      if (variants_.back().gt[original_hap_idx] != alleles[i])
      {
        if (variants_.back().gt[original_hap_idx] == -1)
        {
          // this is the first haplotype mapped to this column, so set allele for new variant
          variants_.back().gt[original_hap_idx] = alleles[i];
        }
        else
        {
          // this haplotype no longer matches its column
          --cardinalities_[original_hap_idx];
          std::size_t j = original_hap_cnt;
          std::size_t j_end = variants_[0].gt.size();
          const std::size_t k_end = variants_.size() - 1;

          for ( ; j < j_end; ++j)
          {
            // check if this haplotype matches any newly created column
            std::size_t k = 0;
            for ( ; k < k_end; ++k)
            {
              if (variants_[k].gt[j] != variants_[k].gt[original_hap_idx])
                break;
            }

            if (k == k_end && variants_[k].gt[j] == alleles[i])
              break;
          }

          if (j == j_end)
          {
            // does not match a newly created column, so insert new column
            for (std::size_t k = 0; k < k_end; ++k)
              variants_[k].gt.push_back(variants_[k].gt[original_hap_idx]);
            variants_[k_end].gt.push_back(alleles[i]);
            cardinalities_.push_back(0);
          }

          unique_map_[i] = j;
          ++cardinalities_[j];
        }
      }
    }
  }

  assert(cardinalities_.size() == variants_.back().gt.size());
  for (std::size_t i = 0; i < cardinalities_.size(); ++i)
  {
    if (variants_.back().gt[i])
      variants_.back().ac += cardinalities_[i];
  }

  assert(std::accumulate(cardinalities_.begin(), cardinalities_.end(), std::size_t(0)) == unique_map_.size());

  return true;
}

void unique_haplotype_block::clear()
{
  variants_.clear();
  unique_map_.clear();
  cardinalities_.clear();
}

void unique_haplotype_block::trim(std::size_t min_pos, std::size_t max_pos)
{
  if (variants_.size())
  {
    if (variants_.front().pos > max_pos || variants_.back().pos < min_pos)
    {
      clear();
    }
    else
    {
      std::size_t before_size = variants_.size();
      auto it = --variants_.end();
      for ( ; it != variants_.begin() && it->pos > max_pos; --it) { }
      variants_.erase(it + 1, variants_.end());

      it = variants_.begin();
      for ( ; it != variants_.end() && it->pos < min_pos; ++it) { }
      variants_.erase(variants_.begin(), it);
      std::size_t after_size = variant_size();
      if (before_size != after_size)
      {
        auto a = 0;
      }
    }
  }
}

void unique_haplotype_block::pop_variant()
{
  variants_.pop_back();
}

void unique_haplotype_block::fill_cm(genetic_map_file& map_file)
{
  for (auto it = variants_.begin(); it != variants_.end(); ++it)
    it->cm = map_file.interpolate_centimorgan(it->pos);
}

bool unique_haplotype_block::serialize(savvy::writer& output_file)
{
  savvy::variant var;
  if (!variants_.empty())
  {
    var = savvy::site_info(variants_.front().chrom,
      variants_.front().pos,
      variants_.front().ref, {"<BLOCK>"});


//    std::int32_t last_end_val;
//    if (variants_.back().get_info("END", last_end_val))
//      var.set_info("END", last_end_val);
//    else
    var.set_info("END", std::int32_t(variants_.back().pos + std::max(variants_.back().ref.size(), variants_.back().alt.size()) - 1));
    var.set_info("VARIANTS", std::int32_t(variants_.size()));
    var.set_info("REPS", std::int32_t(cardinalities_.size()));

    var.set_format("UHM", unique_map_);

    output_file << var;
    output_file.set_block_size(0); // Using set_block_size as a workaround to align zstd blocks with m3vcf blocks.

    for (auto it = variants_.begin(); it != variants_.end() && output_file; ++it)
    {
      var = savvy::variant(it->chrom, it->pos, it->ref, {it->alt});
      var.set_info("AC", std::int32_t(it->ac));
      var.set_info("AN", std::int32_t(unique_map_.size()));
      if (!std::isnan(it->err))
        var.set_info("ERR", it->err);
      if (!std::isnan(it->cm))
        var.set_info("CM", it->cm);
      var.set_info("UHA", it->gt);

      output_file << var;
    }
    output_file.set_block_size(1 + variants_.size());
    return output_file.good();
  }

  return std::cerr << "Error: Cannot write empty block\n", false;
}

int unique_haplotype_block::deserialize(savvy::reader& input_file, savvy::variant& var)
{
  clear();

  if (!input_file.good() || !var.get_format("UHM", unique_map_))
    return input_file.bad() ? -1 : 0;

  std::int64_t n_variants = 0;
  var.get_info("VARIANTS", n_variants);
  std::int64_t n_reps = 0;
  var.get_info("REPS", n_reps);

  cardinalities_.resize(n_reps);
  for (auto it = unique_map_.begin(); it != unique_map_.end(); ++it)
    ++cardinalities_[*it];

  variants_.reserve(n_variants);
  while (input_file >> var)
  {
    if (!var.alts().empty() && var.alts()[0] == "<BLOCK>")
      break;

    variants_.emplace_back();

    variants_.back().chrom = var.chrom();
    variants_.back().pos = var.position();
    variants_.back().ref = var.ref();
    variants_.back().alt = var.alts().size() ? var.alts()[0] : "";
    var.get_info("ERR", variants_.back().err);
    var.get_info("CM", variants_.back().cm);
    var.get_info("UHA", variants_.back().gt);
    if (variants_.back().gt.size() != cardinalities_.size())
      return -1;
    variants_.back().ac = std::inner_product(variants_.back().gt.begin(), variants_.back().gt.end(), cardinalities_.begin(), 0ull);
  }


  if (input_file.bad())
    return -1;
  else
  {
    assert(variants_.size());
    return variants_.size() + 1;
  }
}

bool unique_haplotype_block::deserialize(std::istream& is, int m3vcf_version, std::size_t n_haplotypes)
{
  clear();
  if (is.peek() == EOF)
    return is.get(), false;


  unique_map_.reserve(n_haplotypes);

  auto split_string_to_vector = [](const char* in, char delim)
  {
    std::vector<std::string> ret;
    const char* d = nullptr;
    std::string token;
    const char* s = in;
    const char*const e = in + std::strlen(in);
    while ((d = std::find(s, e,  delim)) != e)
    {
      ret.emplace_back(std::string(s, d));
      s = d ? d + 1 : d;
    }
    ret.emplace_back(std::string(s,d));
    return ret;
  };

  std::string line;
  std::getline(is, line);

  std::size_t n_variants = 0;
  std::size_t n_reps = 0;
  std::size_t col_idx = 0;
  std::size_t last_tab_pos = 0;
  std::size_t tab_pos = 0;
  while ((tab_pos = line.find('\t', tab_pos)) != std::string::npos)
  {
    if (col_idx == 7) //INFO
    {
      std::vector<std::string> info_fields = split_string_to_vector(line.substr(last_tab_pos, tab_pos - last_tab_pos).c_str(), ';');
      for (auto it = info_fields.begin(); it != info_fields.end(); ++it)
      {
        if (it->compare(0, 9, "VARIANTS=") == 0)
          n_variants = std::atoll(it->c_str() + 9);
        else if (it->compare(0, 5, "REPS=") == 0)
          n_reps = std::atoll(it->c_str() + 5);
      }
    }
    else if (col_idx >= 9)
    {
      char* p = nullptr;
      unique_map_.push_back(std::strtoll(line.c_str() + last_tab_pos, &p, 10));
      if (m3vcf_version == 2)
      {
        if (*p != '|')
        {
          clear();
          std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
          return false;
        }
        unique_map_.push_back(std::strtoll(++p, &p, 10));
      }
    }

    last_tab_pos = ++tab_pos;
    ++col_idx;
  }

  {
    char* p = nullptr;
    unique_map_.push_back(std::strtoll(line.c_str() + last_tab_pos, &p, 10));
    if (m3vcf_version == 2)
    {
      if (*p != '|')
      {
        clear();
        is.setstate(is.rdstate() | std::ios::badbit);
        std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
        return false;
      }
      unique_map_.push_back(std::strtoll(++p, &p, 10));
    }
  }

  if (unique_map_.size() != n_haplotypes)
  {
    clear();
    is.setstate(is.rdstate() | std::ios::badbit);
    std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
    return false;
  }

  cardinalities_.resize(n_reps);
  for (auto it = unique_map_.begin(); it != unique_map_.end(); ++it)
    ++cardinalities_[*it];

  variants_.resize(n_variants);
  for (std::size_t i = 0; i < variants_.size(); ++i)
  {
    if (!std::getline(is, line))
    {
      clear();
      is.setstate(is.rdstate() | std::ios::badbit);
      std::cerr << "Error: truncated m3vcf v" << m3vcf_version << " file\n";
      return false;
    }

    std::vector<std::string> cols = split_string_to_vector(line.c_str(), '\t');
    if (cols.size() != 9)
    {
      clear();
      is.setstate(is.rdstate() | std::ios::badbit);
      std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
      return false;
    }

    variants_[i].chrom = cols[0];
    variants_[i].pos = std::atoll(cols[1].c_str());
    variants_[i].ref = cols[3];
    variants_[i].alt = cols[4];

    std::vector<std::string> info_fields = split_string_to_vector(cols[7].c_str(), ';');
    for (auto it = info_fields.begin(); it != info_fields.end(); ++it)
    {
      if (it->compare(0, 4, "ERR=") == 0 || it->compare(0, 4, "Err=") == 0)
        variants_[i].err = std::atof(it->c_str() + 4);
//      else if (it->compare(0, 6, "RECOM=") == 0 || it->compare(0, 6, "Recom=") == 0)
//        variants_[i].recom = std::atof(it->c_str() + 6);
    }

    if (m3vcf_version == 2)
    {
      variants_[i].gt.resize(n_reps);
      const char* p = cols[8].data();
      char* p_end;
      std::size_t prev_offset = 0;
      do
      {
        std::size_t offset = std::strtoll(p, &p_end, 10);
        if (offset >= n_reps)
        {
          clear();
          is.setstate(is.rdstate() | std::ios::badbit);
          std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
          return false;
        }
        std::size_t uniq_idx = prev_offset + offset;
        variants_[i].gt[uniq_idx] = 1;
        variants_[i].ac += cardinalities_[uniq_idx];
        prev_offset = uniq_idx;
        p = p_end + 1;
      } while (*p_end);
    }
    else
    {
      variants_[i].gt.reserve(n_reps);
      for (auto c = cols[8].begin(); c != cols[8].end(); ++c)
        variants_[i].gt.emplace_back((*c) - '0');

      variants_[i].ac = std::inner_product(variants_[i].gt.begin(), variants_[i].gt.end(), cardinalities_.begin(), 0ull);

      if (variants_[i].gt.size() != n_reps)
      {
        clear();
        is.setstate(is.rdstate() | std::ios::badbit);
        std::cerr << "Error: invalid m3vcf v" << m3vcf_version << " file\n";
        return false;
      }
    }
  }

  if (is.good())
    return true;
  is.setstate(is.rdstate() | std::ios::badbit);
  return false;
}

reduced_haplotypes::reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size)
{
  blocks_.emplace_back();
  block_offsets_.emplace_back(0);
  min_block_size_ = std::max(std::size_t(1), min_block_size);
  max_block_size_ = std::max(std::size_t(1), max_block_size);
}

bool reduced_haplotypes::compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles, bool flush_block)
{
  auto comp_ratio = [](const unique_haplotype_block& b)
    {
      return float(b.expanded_haplotype_size() + b.unique_haplotype_size() * b.variant_size()) / float(b.expanded_haplotype_size() * b.variant_size());
    };

  if (flush_block && variant_count_)
  {
    float cr = comp_ratio(blocks_.back());
    block_offsets_.push_back(block_offsets_.back() + blocks_.back().variant_size());
    blocks_.emplace_back();
    return blocks_.back().compress_variant(site_info, alleles);
  }

  float old_cr = comp_ratio(blocks_.back());
  bool ret = blocks_.back().compress_variant(site_info, alleles);
  if (ret)
    ++variant_count_;

  std::size_t cnt = blocks_.back().variant_size();
  if (cnt >= min_block_size_)
  {
    float new_cr = comp_ratio(blocks_.back());
    if (new_cr > old_cr)
    {
      blocks_.emplace_back();
      block_offsets_.push_back(block_offsets_.back() + cnt);
      //std::cerr << "compression old/new: " << old_cr << " / " << new_cr << std::endl;
      return ret;
    }
  }

  //std::cerr << "compression old: " << old_cr << std::endl;

  return ret;
}

void reduced_haplotypes::append_block(const unique_haplotype_block& block)
{
  assert(block.variants().size());
  if (!blocks_.empty())
  {
    const auto& last_var_prev_block = blocks_.back().variants().back();
    const auto& first_var_new_block = block.variants().front();
    if (last_var_prev_block.pos == first_var_new_block.pos && last_var_prev_block.alt == first_var_new_block.alt)
    {
      blocks_.back().pop_variant();
      --variant_count_;
      if (blocks_.back().variant_size() == 0)
      {
        blocks_.pop_back();
        block_offsets_.pop_back();
      }
    }
  }

  block_offsets_.push_back(variant_count_);
  variant_count_ += block.variants().size();
  blocks_.push_back(block);
  assert(blocks_.size() == block_offsets_.size());
}

void reduced_haplotypes::fill_cm(genetic_map_file& map_file)
{
  for (auto it = blocks_.begin(); it != blocks_.end(); ++it)
    it->fill_cm(map_file);
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