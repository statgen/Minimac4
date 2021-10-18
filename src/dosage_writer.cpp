#include "dosage_writer.hpp"

dosage_writer::dosage_writer(const std::string& file_path, const std::string& emp_file_path, const std::string& sites_file_path, savvy::file::format file_format, std::uint8_t out_compression, const std::vector<std::string>& sample_ids, const std::vector<std::string>& fmt_fields, const std::string& chromosome, bool is_temp) :
  out_file_(file_path, file_format, gen_headers(fmt_fields, chromosome, is_temp), sample_ids, out_compression),
  emp_out_file_(emp_file_path.empty() ? nullptr : new savvy::writer(emp_file_path, file_format, gen_emp_headers(chromosome), sample_ids, out_compression)),
  sites_out_file_(sites_file_path.empty() ? nullptr : new savvy::writer(sites_file_path, file_format, gen_headers({}, chromosome, false), {}, out_compression)),
  file_format_(file_format),
  fmt_fields_(fmt_fields),
  fmt_field_set_(fmt_fields.begin(), fmt_fields.end()),
  n_samples_(sample_ids.size()),
  is_temp_file_(is_temp)
{

}

std::vector<std::pair<std::string, std::string>> dosage_writer::gen_headers(const std::vector<std::string>& fmt_fields, const std::string& chromosome, bool is_temp)
{
  std::time_t t = std::time(nullptr);
  char datestr[11];
  std::string filedate(datestr, std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));
  assert(filedate.size());

  std::vector<std::pair<std::string, std::string>> headers;

  if (is_temp)
  {
    headers = {
      {"fileformat","VCFv4.2"},
      {"filedate", filedate},
      {"source","Minimac v" + std::string(VERSION)},
      {"phasing","full"},
      {"contig","<ID=" + std::string(chromosome) + ">"},
      {"INFO","<ID=S_X,Number=1,Type=Float,Description=\"Sum of dosages\">"},
      {"INFO","<ID=S_XX,Number=1,Type=Float,Description=\"Sum of squared dosages\">"},
      {"INFO","<ID=S_CS,Number=1,Type=Float,Description=\"Sum of call scores\">"},
      {"INFO","<ID=LOO_S_X,Number=1,Type=Float,Description=\"Sum of LOO dosages\">"},
      {"INFO","<ID=LOO_S_XX,Number=1,Type=Float,Description=\"Sum of squared LOO dosages\">"},
      {"INFO","<ID=LOO_S_Y,Number=1,Type=Float,Description=\"Sum of observed genotypes\">"},
      {"INFO","<ID=LOO_S_YY,Number=1,Type=Float,Description=\"Sum of squared observed genotypes\">"},
      {"INFO","<ID=LOO_S_XY,Number=1,Type=Float,Description=\"Dot product of LOO dosages and observed genotypes\">"},
      {"INFO","<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">"},
      {"INFO","<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">"},
      {"INFO","<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Marker was genotyped but NOT imputed\">"}};
  }
  else
  {
    headers = {
      {"fileformat", "VCFv4.2"},
      {"filedate", filedate},
      {"source", "Minimac v" + std::string(VERSION)},
      {"phasing", "full"},
      {"contig", "<ID=" + std::string(chromosome) + ">"},
      {"INFO", "<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">"},
      {"INFO", "<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">"},
      {"INFO", "<ID=AVGCS,Number=1,Type=Float,Description=\"Average Call Score\">"},
      {"INFO", "<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">"},
      {"INFO", "<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">"},
      {"INFO", "<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">"},
      {"INFO", "<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">"},
      {"INFO", "<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Marker was genotyped but NOT imputed\">"}};
  }

  headers.reserve(headers.size() + 5);
  for (const auto& f : fmt_fields)
  {
    if (f == "GT")
      headers.emplace_back("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    else if (f == "DS")
      headers.emplace_back("FORMAT", "<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">");
    else if (f == "HDS")
      headers.emplace_back("FORMAT", "<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">");
    else if (f == "GP")
      headers.emplace_back("FORMAT", "<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">");
    else if (f == "SD")
      headers.emplace_back("FORMAT", "<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">");
  }

  // TODO: command string

  return headers;
}

std::vector<std::pair<std::string, std::string>> dosage_writer::gen_emp_headers(const std::string& chromosome)
{
  std::time_t t = std::time(nullptr);
  char datestr[11];
  std::string filedate(datestr, std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));
  assert(filedate.size());

  std::vector<std::pair<std::string, std::string>> headers = {
    {"fileformat", "VCFv4.2"},
    {"filedate", filedate},
    {"source", "Minimac v" + std::string(VERSION)},
    {"phasing", "full"},
    {"contig", "<ID=" + std::string(chromosome) + ">"},
    {"INFO", "<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">"},
    {"INFO", "<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">"},
    {"INFO", "<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Marker was genotyped but NOT imputed\">"},
    {"FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotyped alleles from Array\">"},
    {"FORMAT","<ID=LDS,Number=2,Type=Float,Description=\"Leave-one-out Imputed Dosage : Estimated Haploid Alternate Allele Dosage assuming site was NOT genotyped\">"}};

  // TODO: command string

  return headers;
}

bool dosage_writer::sites_match(const target_variant& t, const reference_site_info& r)
{
  return t.pos == r.pos && t.alt == r.alt && t.ref == r.ref;
}

bool dosage_writer::merge_temp_files(std::list<std::string>& temp_file_paths, std::list<std::string>& temp_emp_file_paths)
{
  std::list<savvy::reader> temp_files(temp_file_paths.begin(), temp_file_paths.end());
  std::list<savvy::reader> temp_emp_files(temp_emp_file_paths.begin(), temp_emp_file_paths.end());
  std::for_each(temp_file_paths.begin(), temp_file_paths.end(), [](std::string&s) { std::remove(s.c_str()); });
  std::for_each(temp_emp_file_paths.begin(), temp_emp_file_paths.end(), [](std::string&s) { std::remove(s.c_str()); });
  return merge_temp_files(temp_files, temp_emp_files);
}

bool dosage_writer::merge_temp_files(std::list<savvy::reader>& temp_files, std::list<savvy::reader>& temp_emp_files)
{
  if (temp_files.empty())
    return false;
#if 0 // For when temp file readers are created before fully written.
  for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
    it->reset_bounds(savvy::slice_bounds(0));

  for (auto it = temp_emp_files.begin(); it != temp_emp_files.end(); ++it)
    it->reset_bounds(savvy::slice_bounds(0));
#endif
  savvy::variant out_var;
  savvy::variant out_var_emp;
  savvy::compressed_vector<float> pasted_hds;
  savvy::compressed_vector<float> partial_hds;

  std::vector<float> pasted_lds;
  std::vector<float> partial_lds;
  std::vector<std::int8_t> pasted_gt;
  std::vector<std::int8_t> partial_gt;

  int good_count = temp_files.size();
  while (good_count == temp_files.size())
  {
    float s_cs{}, s_x{}, s_xx{}, loo_s_x{}, loo_s_xx{}, loo_s_y{}, loo_s_yy{}, loo_s_xy{};
    bool is_typed = false;

    pasted_hds.clear();

    good_count = 0;
    for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
    {
      good_count += (int)it->read(out_var).good();
      out_var.get_format("HDS", partial_hds);
      std::size_t old_size = pasted_hds.size();
      pasted_hds.resize(old_size + partial_hds.size());
      for (auto jt = partial_hds.begin(); jt != partial_hds.end(); ++jt)
        pasted_hds[old_size + jt.offset()] = *jt;

      float tmp;
      if (out_var.get_info("S_X", tmp)) s_x += tmp;
      if (out_var.get_info("S_XX", tmp)) s_xx += tmp;
      if (out_var.get_info("S_CS", tmp)) s_cs += tmp;
      if (out_var.get_info("LOO_S_X", tmp))
      {
        loo_s_x += tmp;
        if (out_var.get_info("LOO_S_XX", tmp)) loo_s_xx += tmp;
        if (out_var.get_info("LOO_S_Y", tmp)) loo_s_y += tmp;
        // if (out_var.get_info("LOO_S_YY", tmp)) loo_s_yy += tmp;
        if (out_var.get_info("LOO_S_XY", tmp)) loo_s_xy += tmp;
        is_typed = true;
      }
    }

    loo_s_yy = loo_s_y;

    if (good_count)
    {
      if (good_count < temp_files.size())
        return std::cerr << "Error: record mismatch in temp files" << std::endl, false;

      out_var.remove_info("S_X");
      out_var.remove_info("S_XX");
      out_var.remove_info("S_CS");

      std::size_t n = pasted_hds.size();
//      float s_x = std::accumulate(pasted_hds.begin(), pasted_hds.end(), 0.f);
//      float s_xx = std::inner_product(pasted_hds.begin(), pasted_hds.end(), pasted_hds.begin(), 0.f);

      float af = s_x / n;
      out_var.set_info("AF", af);
      out_var.set_info("MAF", af > 0.5f ? 1.f - af : af);
      out_var.set_info("AVGCS", s_cs / n);

      set_r2_info_field(out_var, s_x, s_xx, n);

      if (is_typed)
      {
        assert(std::find_if(out_var.info_fields().begin(), out_var.info_fields().end(), [](const std::pair<std::string, savvy::typed_value>& v) { return v.first == "TYPED"; }) != out_var.info_fields().end());
        out_var.remove_info("LOO_S_X");
        out_var.remove_info("LOO_S_XX");
        out_var.remove_info("LOO_S_Y");
        out_var.remove_info("LOO_S_YY");
        out_var.remove_info("LOO_S_XY");

        set_er2_info_field(out_var, loo_s_x, loo_s_xx, loo_s_y, loo_s_yy, loo_s_xy, n);

        if (emp_out_file_)
        {
          pasted_lds.clear();
          pasted_gt.clear();
          pasted_lds.reserve(n);
          pasted_gt.reserve(n);

          int good_count_emp = 0;
          for (auto it = temp_emp_files.begin(); it != temp_emp_files.end(); ++it)
          {
            good_count_emp += (int)it->read(out_var_emp).good();

            out_var_emp.get_format("LDS", partial_lds);
            for (auto jt = partial_lds.begin(); jt != partial_lds.end(); ++jt)
              pasted_lds.push_back(*jt);

            out_var_emp.get_format("GT", partial_gt);
            for (auto jt = partial_gt.begin(); jt != partial_gt.end(); ++jt)
              pasted_gt.push_back(*jt);
          }

          if (good_count_emp < temp_emp_files.size())
            return std::cerr << "Error: record mismatch in empirical temp files" << std::endl, false;

          out_var_emp.set_format("GT", pasted_gt);
          out_var_emp.set_format("LDS", pasted_lds);
          emp_out_file_->write(out_var_emp);
        }
      }

      if (sites_out_file_)
      {
        savvy::variant site_var;
        dynamic_cast<savvy::site_info&>(site_var) = out_var;
        sites_out_file_->write(site_var);
      }

      set_format_fields(out_var, pasted_hds);
      out_file_ << out_var;
    }
  }

  int bad_count = 0;
  for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
    bad_count += (int)it->bad();

  if (bad_count || !out_file_.good())
    return std::cerr << "Error: I/O failed while merging" << std::endl, false;

  if (emp_out_file_)
  {
    bad_count = 0;
    for (auto it = temp_emp_files.begin(); it != temp_emp_files.end(); ++it)
      bad_count += (int)it->bad();

    if (bad_count || !emp_out_file_->good())
      return std::cerr << "Error: I/O failed while merging empirical" << std::endl, false;
  }

  return true;
}

bool dosage_writer::write_dosages(const full_dosages_results& hmm_results, const std::vector<target_variant>& tar_variants, std::pair<std::size_t, std::size_t> observed_range, const reduced_haplotypes& full_reference_data)
{
  assert(hmm_results.dimensions()[0] == full_reference_data.variant_size());

  variant_update_ctx update_context;
  savvy::variant out_var;
  savvy::variant out_var_emp;

  assert(!tar_variants.empty());
//  std::vector<float> dosages(tar_variants[0].gt.size());
//  std::vector<float> loo_dosages(tar_variants[0].gt.size());
  savvy::compressed_vector<float> sparse_dosages;
  const std::vector<std::int8_t> empty_gt_vec;
  const std::vector<float> empty_loo_vec;

  assert(!is_temp_file_ || observed_range.first != 0 || observed_range.second != tar_variants[0].gt.size());

  std::size_t tar_idx = 0;
  std::size_t i = 0;

  for (auto ref_it = full_reference_data.begin(); ref_it != full_reference_data.end(); ++ref_it,++i)
  {
    assert(i == ref_it.global_idx());
    while (tar_idx < tar_variants.size() && tar_variants[tar_idx].pos < ref_it->pos)
      ++tar_idx;

    bool ref_matches_tar = tar_idx < tar_variants.size() && sites_match(tar_variants[tar_idx], *ref_it);
    //      for (std::size_t j = 0; j < dosages.size(); ++j)
    //        dosages[j] = float(std::int16_t(hmm_results.dosage(i, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
    //
    //      if (ref_matches_tar)
    //      {
    //        for (std::size_t j = 0; j < loo_dosages.size(); ++j)
    //          loo_dosages[j] = float(std::int16_t(hmm_results.loo_dosage(tar_idx, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
    //      }



    out_var = savvy::site_info(ref_it->chrom, ref_it->pos, ref_it->ref, {ref_it->alt}, ""/*ref_var.id()*/);
    sparse_dosages.assign(hmm_results.dosages_[i].begin(), hmm_results.dosages_[i].end());
    if (ref_matches_tar)
    {
      std::vector<std::int8_t> observed(tar_variants[tar_idx].gt.begin() + observed_range.first, tar_variants[tar_idx].gt.begin() + observed_range.second);
      set_info_fields(out_var, sparse_dosages, hmm_results.loo_dosages_[tar_idx], observed);
      if (emp_out_file_)
      {
        out_var_emp = savvy::site_info(ref_it->chrom, ref_it->pos, ref_it->ref, {ref_it->alt}, ""/*ref_var.id()*/);
        out_var_emp.set_info("TYPED", std::vector<std::int8_t>());
        out_var_emp.set_format("GT", observed);
        out_var_emp.set_format("LDS", hmm_results.loo_dosages_[tar_idx]);
        emp_out_file_->write(out_var_emp);
      }
    }
    else
    {
      set_info_fields(out_var, sparse_dosages, {}, {});
    }

    if (sites_out_file_)
    {
      savvy::variant site_var;
      dynamic_cast<savvy::site_info&>(site_var) = out_var;
      sites_out_file_->write(site_var);
    }

    set_format_fields(out_var, sparse_dosages);
    out_file_ << out_var;
  }

  return out_file_.good();
}

void dosage_writer::set_r2_info_field(savvy::variant& out_var, float s_x, float s_xx, std::size_t n)
{
  float af = s_x / float(n);
  float r2 = savvy::typed_value::missing_value<float>();
  if (af > 0.f && af < 1.f)
    r2 = ((s_xx - s_x * s_x / n) / n) / (af * (1.f - af));

  out_var.set_info("R2", r2);
}

void dosage_writer::set_er2_info_field(savvy::variant& out_var, float s_x, float s_xx, float s_y, float s_yy, float s_xy, std::size_t n)
{
  //                         n * Sum xy - Sum x * Sum y
  //  r = -------------------------------------------------------------------
  //      Sqrt(n * Sum xx - Sum x * Sum x) * Sqrt(n * Sum yy - Sum y * Sum y)
  float emp_r = (n * s_xy - s_x * s_y) / (std::sqrt(n * s_xx - s_x * s_x) * std::sqrt(n * s_yy - s_y * s_y));
  out_var.set_info("ER2", std::isnan(emp_r) ? savvy::typed_value::missing_value<float>() : emp_r * emp_r);
}

void dosage_writer::set_info_fields(savvy::variant& out_var, const savvy::compressed_vector<float>& sparse_dosages, const std::vector<float>& loo_dosages, const std::vector<std::int8_t>& observed)
{
  std::size_t n = sparse_dosages.size();

  float s_x = std::accumulate(sparse_dosages.begin(), sparse_dosages.end(), 0.f);
  float s_xx = std::inner_product(sparse_dosages.begin(), sparse_dosages.end(), sparse_dosages.begin(), 0.f);
  float af = s_x / n;
  float s_cs(sparse_dosages.size() - sparse_dosages.non_zero_size());
  for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
    s_cs += *it > 0.5f ? *it : 1.f - *it;

  if (is_temp_file_)
  {
    out_var.set_info("S_X", s_x);
    out_var.set_info("S_XX", s_xx);
    out_var.set_info("S_CS", s_cs);
  }
  else
  {
    out_var.set_info("AF", af);
    out_var.set_info("MAF", af > 0.5f ? 1.f - af : af);
    out_var.set_info("AVGCS", s_cs / n);
    set_r2_info_field(out_var, s_x, s_xx, n);
  }

  if (observed.size())
  {
    assert(observed.size() == loo_dosages.size());
    //sparse_loo_dosages.assign(loo_dosages.begin(), loo_dosages.end());

    s_x = std::accumulate(loo_dosages.begin(), loo_dosages.end(), 0.f);
    s_xx = std::inner_product(loo_dosages.begin(), loo_dosages.end(), loo_dosages.begin(), 0.f);
    float s_y = std::accumulate(observed.begin(), observed.end(), 0.f);
    // since observed can only be 0 or 1, s_yy is the same as s_y
    float s_yy = s_y; //std::inner_product(sparse_gt.begin(), sparse_gt.end(), sparse_gt.begin(), 0.f); // TODO: allow for missing oberserved genotypes.
    float s_xy = 0.f;
    s_xy = std::inner_product(loo_dosages.begin(), loo_dosages.end(), observed.begin(), s_xy);
//    for (auto it = ctx.sparse_gt.begin(); it != ctx.sparse_gt.end(); ++it)
//      s_xy += *it * loo_dosages[it.offset()];

    if (is_temp_file_)
    {
      out_var.set_info("LOO_S_X", s_x);
      out_var.set_info("LOO_S_XX", s_xx);
      out_var.set_info("LOO_S_Y", s_y);
      out_var.set_info("LOO_S_YY", s_yy);
      out_var.set_info("LOO_S_XY", s_xy);
    }
    else
    {
      set_er2_info_field(out_var, s_x, s_xx, s_y, s_yy, s_xy, n);
    }
  }

  out_var.set_info(observed.size() ? "TYPED" : "IMPUTED", std::vector<std::int8_t>());
}

void dosage_writer::set_format_fields(savvy::variant& out_var, savvy::compressed_vector<float>& sparse_dosages)
{
  std::size_t stride = sparse_dosages.size() / n_samples_;

  if (fmt_field_set_.find("GT") != fmt_field_set_.end())
  {
    out_var.set_format("HDS", {});

    sparse_gt_.assign(sparse_dosages.value_data(), sparse_dosages.value_data() + sparse_dosages.non_zero_size(), sparse_dosages.index_data(), sparse_dosages.size(), [](float v)
      {
        if (savvy::typed_value::is_end_of_vector(v))
          return savvy::typed_value::end_of_vector_value<std::int8_t>();
        return std::int8_t(v < 0.5f ? 0 : 1);
      });
    out_var.set_format("GT", sparse_gt_);
  }

  if (fmt_field_set_.find("HDS") != fmt_field_set_.end())
  {
    out_var.set_format("HDS", sparse_dosages);
  }
  else
  {
    out_var.set_format("HDS", {});
  }

  if (fmt_field_set_.find("GP") != fmt_field_set_.end() || fmt_field_set_.find("SD") != fmt_field_set_.end())
  {
    // set dense dosage vector
    dense_zero_vec_.resize(sparse_dosages.size());
    for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
      dense_zero_vec_[it.offset()] = *it;

    std::vector<float>& dense_hds = dense_zero_vec_;

    if (fmt_field_set_.find("GP") != fmt_field_set_.end())
    {
      if (stride == 1)
      {
        // All samples are haploid
        dense_float_vec_.resize(n_samples_ * 2);
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          std::size_t dest_idx = i * 2;
          dense_float_vec_[dest_idx] = 1.f - dense_hds[i];
          dense_float_vec_[dest_idx + 1] = dense_hds[i];
        }
      }
      else if (stride == 2)
      {
        dense_float_vec_.resize(n_samples_ * 3);
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          std::size_t src_idx = i * 2;
          std::size_t dest_idx = i * 3;
          float x = dense_hds[src_idx];
          float y = dense_hds[src_idx + 1];
          if (savvy::typed_value::is_end_of_vector(y))
          {
            // haploid
            dense_float_vec_[dest_idx] = 1.f - x;
            dense_float_vec_[dest_idx + 1] = x;
            dense_float_vec_[dest_idx + 2] = y;
          }
          else
          {
            // diploid
            dense_float_vec_[dest_idx] = (1.f - x) * (1.f - y);
            dense_float_vec_[dest_idx + 1] = x * (1.f - y) + y * (1.f - x);
            dense_float_vec_[dest_idx + 2] = x * y;
          }
        }
      }

      out_var.set_format("GP", dense_float_vec_);
    }

    if (fmt_field_set_.find("SD") != fmt_field_set_.end())
    {
      dense_float_vec_.resize(n_samples_);
      if (stride == 1)
      {
        // All samples are haploid
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          dense_float_vec_[i] = dense_hds[i] * (1.f - dense_hds[i]);
        }

        out_var.set_format("SD", dense_float_vec_);
      }
      else if (stride == 2)
      {
        for (std::size_t i = 0; i < dense_hds.size(); i += 2)
        {
          float x = dense_hds[i];
          float y = dense_hds[i + 1];
          if (savvy::typed_value::is_end_of_vector(y)) // haploid
            dense_float_vec_[i / 2] = x * (1.f - x);
          else // diploid
            dense_float_vec_[i / 2] = x * (1.f - x) + y * (1.f - y);
        }

        out_var.set_format("SD", dense_float_vec_);
      }
      else
      {
        // TODO: suppress error excessive error messages
        std::cerr << "Error: only haploid and diploid samples are supported when generating SD\n";
      }
    }

    // unset dense dosage vector
    for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
      dense_zero_vec_[it.offset()] = 0.f;
  }

  if (fmt_field_set_.find("DS") != fmt_field_set_.end())
  {
    savvy::stride_reduce(sparse_dosages, sparse_dosages.size() / n_samples_);
    out_var.set_format("DS", sparse_dosages);
  }
}
