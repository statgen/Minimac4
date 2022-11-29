#include "hidden_markov_model.hpp"
#include "variant.hpp"

#include <cassert>
#include <numeric>
#include <unordered_map>

constexpr float hidden_markov_model::jump_fix;
constexpr float hidden_markov_model::jump_threshold;

hidden_markov_model::hidden_markov_model(float s3_prob_threshold, float s1_prob_threshold, float diff_threshold, float background_error) :
  prob_threshold_(s3_prob_threshold),
  s1_prob_threshold_(s1_prob_threshold),
  diff_threshold_(diff_threshold),
  background_error_(background_error)
{
}

void hidden_markov_model::initialize_likelihoods(std::vector<float>& probs, std::vector<float>& probs_norecom, std::vector<float>& proportions, const unique_haplotype_block& ref_block)
{
  probs.clear();
  probs.resize(ref_block.cardinalities().size(), 0.f);
  for (std::size_t i = 0; i < probs.size(); ++i)
    probs[i] = 1.f * ref_block.cardinalities()[i];

  probs_norecom = probs;

  proportions.resize(ref_block.expanded_haplotype_size());
  for (std::size_t i = 0; i < proportions.size(); ++i)
    proportions[i] = 1.f / ref_block.cardinalities()[ref_block.unique_map()[i]];
}

void hidden_markov_model::traverse_forward(const std::deque<unique_haplotype_block>& ref_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t hap_idx)
{
  std::size_t n_expanded_haplotypes = ref_haps.front().expanded_haplotype_size();

  precision_jumps_.resize(tar_variants.size());
  junction_prob_proportions_.resize(ref_haps.size(), std::vector<float>(n_expanded_haplotypes));
  forward_probs_.resize(ref_haps.size());
  forward_norecom_probs_.resize(ref_haps.size());
  for (std::size_t b = 0; b < forward_probs_.size(); ++b)
  {
    auto& prob_block = forward_probs_[b];
    auto& norecom_prob_block = forward_norecom_probs_[b];
    const auto& ref_block = ref_haps[b];

    prob_block.resize(ref_block.variant_size());
    norecom_prob_block.resize(ref_block.variant_size());
    for (std::size_t v = 0; v < ref_block.variant_size(); ++v)
    {
      prob_block[v].resize(ref_block.unique_haplotype_size());
      norecom_prob_block[v].resize(ref_block.unique_haplotype_size());
    }
  }

  std::size_t global_idx = 0;
  std::vector<float> temp_row;
  for (std::size_t block_idx = 0; block_idx < ref_haps.size(); ++block_idx,++global_idx)
  {
    const unique_haplotype_block& ref_block = ref_haps[block_idx];
    if (block_idx == 0)
    {
      // Initialize likelihoods at first position
      initialize_likelihoods(forward_probs_[block_idx][0], forward_norecom_probs_[block_idx][0], junction_prob_proportions_[block_idx], ref_block);
    }
    else
    {
      // Transition from previous block
      const unique_haplotype_block& prev_ref_block = ref_haps[block_idx - 1];
      std::vector<float>& prev_junction_proportions = junction_prob_proportions_[block_idx - 1];
      std::vector<float>& prev_last_row = forward_probs_[block_idx - 1].back();
      std::vector<float>& prev_last_row_norecom = forward_norecom_probs_[block_idx - 1].back();
      std::vector<float>& junction_proportions = junction_prob_proportions_[block_idx];
      std::vector<float>& first_row = forward_probs_[block_idx].front();
      std::vector<float>& first_row_norecom = forward_norecom_probs_[block_idx].front();
      temp_row.clear();
      temp_row.resize(first_row.size(), 0.f);
      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        std::size_t prev_uniq_idx = prev_ref_block.unique_map()[i];
        float p = prev_last_row_norecom[prev_uniq_idx] * prev_junction_proportions[i] + (prev_last_row[prev_uniq_idx] - prev_last_row_norecom[prev_uniq_idx]) / prev_ref_block.cardinalities()[prev_uniq_idx];
        junction_proportions[i] = p;
        temp_row[uniq_idx] += p;
        assert((prev_last_row[prev_uniq_idx] - prev_last_row_norecom[prev_uniq_idx]) >= 0.f);
      }

      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        assert(temp_row[uniq_idx] > 0.f);
        junction_proportions[i] = junction_proportions[i] / temp_row[uniq_idx];
        assert(junction_proportions[i] >= 0.f);
        assert(junction_proportions[i] <= 1.f);
      }

      precision_jumps_[global_idx - 1] = transpose(temp_row, first_row, temp_row, first_row_norecom, ref_block.cardinalities(), tar_variants[global_idx - 1].recom, n_expanded_haplotypes);
    }


    const auto& template_variants = ref_block.variants();
    std::size_t n_rows = forward_probs_[block_idx].size();
    std::size_t last_row_idx = n_rows - 1;
    for (std::size_t i = 0; i < last_row_idx; ++i,++global_idx)
    {
      std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
      if (observed >= 0)
        condition(forward_probs_[block_idx][i], forward_norecom_probs_[block_idx][i], template_variants[i].gt, observed, tar_variants[global_idx].err, tar_variants[global_idx].af);
      precision_jumps_[global_idx] = transpose(forward_probs_[block_idx][i], forward_probs_[block_idx][i + 1], forward_norecom_probs_[block_idx][i], forward_norecom_probs_[block_idx][i + 1], ref_block.cardinalities(), tar_variants[global_idx].recom, n_expanded_haplotypes);
    }

    std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
    if (observed >= 0)
      condition(forward_probs_[block_idx][last_row_idx], forward_norecom_probs_[block_idx][last_row_idx], template_variants[last_row_idx].gt, observed, tar_variants[global_idx].err, tar_variants[global_idx].af);
  }
}

void hidden_markov_model::traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t hap_idx,
  std::size_t out_idx,
  const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
  full_dosages_results& output,
  const reduced_haplotypes& full_reference_data)
{
  std::size_t prev_full_ref_block_idx(-1);
  auto full_ref_ritr = --full_reference_data.end();
  const auto full_ref_rend = --full_reference_data.begin();

  std::size_t global_idx = tar_variants.size() - 1;
  std::vector<float> backward;
  std::vector<float> backward_norecom;
  std::vector<float> junction_proportions_backward;
  std::vector<float> extra;
  std::vector<float> constants;

  double prob_sum = std::accumulate(forward_probs_.back().back().begin(), forward_probs_.back().back().end(), 0.);

  int last_block_idx = int(ref_haps.size()) - 1;
  for (int block_idx = last_block_idx; block_idx >= 0; --block_idx)
  {
    const unique_haplotype_block& ref_block = ref_haps[block_idx];
    std::size_t n_expanded_haplotypes = ref_block.expanded_haplotype_size();
    if (block_idx == last_block_idx)
    {
      // Initialize likelihoods at first position
      initialize_likelihoods(backward, backward_norecom, junction_proportions_backward, ref_block);
    }
    else
    {
      // Transition from previous block
      assert(block_idx + 1 < ref_haps.size());
      const unique_haplotype_block& prev_ref_block = ref_haps[block_idx + 1];
      extra.clear();
      extra.resize(ref_block.unique_haplotype_size(), 0.f);
      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        std::size_t prev_uniq_idx = prev_ref_block.unique_map()[i];
        float p = backward_norecom[prev_uniq_idx] * junction_proportions_backward[i] + (backward[prev_uniq_idx] - backward_norecom[prev_uniq_idx]) / prev_ref_block.cardinalities()[prev_uniq_idx];
        junction_proportions_backward[i] = p;
        extra[uniq_idx] += p;
        assert((backward[prev_uniq_idx] - backward_norecom[prev_uniq_idx]) >= 0.f);
      }

      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        assert(extra[uniq_idx] > 0.f);
        junction_proportions_backward[i] /= extra[uniq_idx];
        assert(junction_proportions_backward[i] >= 0.f);
      }

      std::swap(backward, extra);
      backward_norecom = backward;
    }

    best_s3_haps_.clear();
    best_s3_probs_.clear();

    constants.clear();
    constants.resize(reverse_maps[block_idx].size());
    for (std::size_t i = 0; i < reverse_maps[block_idx].size(); ++i)
    {
      for (std::size_t j = 0; j < reverse_maps[block_idx][i].size(); ++j)
        constants[i] += junction_prob_proportions_[block_idx][reverse_maps[block_idx][i][j]] * junction_proportions_backward[reverse_maps[block_idx][i][j]];
    }

    std::size_t best_hap(-1);

    extra.resize(backward.size());

    const auto& template_variants = ref_block.variants();
    std::size_t n_rows = forward_probs_[block_idx].size();
    for (int i = int(n_rows) - 1; i >= 0; --i,--global_idx)
    {
      bool right_jump = transpose(backward, backward, backward_norecom, backward_norecom, ref_block.cardinalities(), tar_variants[global_idx].recom, n_expanded_haplotypes);

      if (global_idx > 0 && precision_jumps_[global_idx - 1])
        auto a = 0;
      if (precision_jumps_[global_idx])
        prob_sum /= jump_fix;
      if (right_jump)
        prob_sum *= jump_fix;

      std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
      impute(prob_sum, best_hap,
        forward_probs_[block_idx][i], backward,
        forward_norecom_probs_[block_idx][i], backward_norecom,
        junction_prob_proportions_[block_idx], junction_proportions_backward,
        constants, ref_block.unique_map(), reverse_maps[block_idx],
        template_variants[i].gt,
        tar_variants,
        global_idx, hap_idx, out_idx,
        output,
        full_ref_ritr, full_ref_rend, prev_full_ref_block_idx);

      if (observed >= 0)
        condition(backward, backward_norecom, template_variants[i].gt, observed, tar_variants[global_idx].err, tar_variants[global_idx].af);
    }
  }
  assert(global_idx == std::size_t(-1));
}

bool hidden_markov_model::transpose(const std::vector<float>& from, std::vector<float>& to, const std::vector<float>& from_norecom, std::vector<float>& to_norecom, const std::vector<std::size_t>& uniq_cardinalities, double recom, std::size_t n_templates)
{
  bool jumped = false;
  assert(from.size() == to.size());

  double sum = 0.;

  for (std::size_t i = 0; i < from.size(); ++i)
  {
    sum += from[i];
  }

  double complement = 1. - recom;

  // prevent probs from getting too small
  if (sum < jump_threshold) // 1e-10
  {
    sum *= jump_fix; // 1e15
    complement *= jump_fix; // 1e15
    jumped = true;
  }

  sum *= (recom / n_templates);

  for (int i = 0; i < to.size(); i++)
  {
    to[i] = from[i] * complement + (uniq_cardinalities[i] * sum);
    to_norecom[i] = from_norecom[i] * complement;
    assert(to[i] >= 0.0f);
    assert(to_norecom[i] >= 0.0f);
    assert(to[i] < 1e18);
    assert(to_norecom[i] < 1e18);
  }

  return jumped;
}

void hidden_markov_model::condition(std::vector<float>& probs, std::vector<float>& probs_norecom, const std::vector<std::int8_t>& template_haps, std::int8_t observed, float err, float af)
{
  float prandom = err * (observed ? af : 1.f - af) + background_error_;
  float pmatch = (1.f - err) + prandom;

  for (std::size_t i = 0; i < probs.size(); ++i)
  {
    if(observed == template_haps[i])
    {
      probs[i] *= pmatch;
      probs_norecom[i] *= pmatch;
    }
    else
    {
      probs[i] *= prandom;
      probs_norecom[i] *= prandom;
    }

    assert(probs[i]>=0.0f);
  }
}

void hidden_markov_model::impute_typed_site(double& prob_sum, std::size_t& prev_best_hap,
  const std::vector<float>& left_probs,
  const std::vector<float>& right_probs,
  const std::vector<float>& left_probs_norecom,
  const std::vector<float>& right_probs_norecom,
  const std::vector<float>& left_junction_proportions,
  const std::vector<float>& right_junction_proportions,
  const std::vector<float>& constants,
  const std::vector<std::vector<std::size_t>>& reverse_map,
  const std::vector<std::int8_t>& template_haps,
  std::int8_t observed, float err, float af,
  std::vector<std::uint32_t>& best_unique_haps, std::vector<float>& best_unique_probs, float& dose, float& loo_dose)
{
  assert(left_probs.size() == right_probs.size());

  if (prev_best_hap < constants.size())
  {
    std::size_t i = prev_best_hap;
    double left_nr = left_probs_norecom[i];
    double right_nr = right_probs_norecom[i];
    double n = reverse_map[i].size();
    double prob = constants[i] * left_nr * right_nr + (double(left_probs[i]) * double(right_probs[i]) - left_nr * right_nr) / n;
    prob /= prob_sum;
    if (prob > (1. - prob_threshold_))
    {
      best_unique_probs.push_back(prob);
      best_unique_haps.push_back(i);
//      if (template_haps[i])
//        p_alt = std::min(1.f, std::max(0.f, prob)); // TODO: + (1. - p_alt) * AF_other to support larger thresholds
      dose = template_haps[i] ? 1.f : 0.f;
      loo_dose = savvy::typed_value::is_end_of_vector(observed) ? savvy::typed_value::end_of_vector_value<float>() : dose;
      return;
    }
  }

  assert(best_unique_haps.empty() && best_unique_probs.empty());
  double p_ref = 0.;
  double p_alt = 0.;
  double denorm_threshold = prob_threshold_ * prob_sum;
  for (std::size_t i = 0; i < constants.size(); ++i)
  {
    double left_nr = left_probs_norecom[i];
    double right_nr = right_probs_norecom[i];
    double n = reverse_map[i].size();
    double prob = constants[i] * left_nr * right_nr + (double(left_probs[i]) * double(right_probs[i]) - left_nr * right_nr) / n;
    if (template_haps[i])
      p_alt += prob;
    else
      p_ref += prob;

    if (prob > denorm_threshold)
    {
      best_unique_probs.push_back(prob / prob_sum);
      best_unique_haps.push_back(i);
      assert(prob / prob_sum > prob_threshold_); // TODO: remove since rounding error could make this assert fail.
    }
  }
  prob_sum = p_alt + p_ref;
  dose = std::min(1.f, std::max(0.f, float(p_alt / prob_sum)));
  dose = float(std::int16_t(dose * bin_scalar_ + 0.5f)) / bin_scalar_; // bin

  if (savvy::typed_value::is_end_of_vector(observed))
  {
    loo_dose = savvy::typed_value::end_of_vector_value<float>();
  }
  else
  {
    double loo_p_alt = p_alt;
    double loo_p_ref = p_ref;

    float fmismatch = err * (observed ? af : 1.f - af) + background_error_;
    float fmatch = 1.f - err + fmismatch;
    assert(fmismatch > 0.f && fmatch > 0.f);

    if (observed == 1)
    {
      loo_p_alt /= fmatch;
      loo_p_ref /= fmismatch;
    }
    else
    {
      loo_p_ref /= fmatch;
      loo_p_alt /= fmismatch;
    }

    loo_dose = float(loo_p_alt / (loo_p_alt + loo_p_ref));
    loo_dose = float(std::int16_t(loo_dose * bin_scalar_ + 0.5f)) / bin_scalar_; // bin
  }

  prev_best_hap = best_unique_haps.size() == 1 ? best_unique_haps.front() : std::numeric_limits<std::size_t>::max();
}

inline bool sites_match(const target_variant& t, const reference_site_info& r)
{
  return t.pos == r.pos && t.alt == r.alt && t.ref == r.ref;
}

void hidden_markov_model::s3_to_s1_probs(
  const std::vector<float>& left_probs, const std::vector<float>& right_probs,
  const std::vector<float>& left_probs_norecom, const std::vector<float>& right_probs_norecom,
  const std::vector<float>& left_junction_proportions, const std::vector<float>& right_junction_proportions,
  const std::vector<std::vector<std::size_t>>& s3_reverse_map, double prob_sum)
{
  best_s1_haps_.clear();
  best_s1_probs_.clear();
  std::size_t n_templates = left_junction_proportions.size();
  double denorm_threshold = prob_sum * std::min(prob_threshold_, s1_prob_threshold_ < 0.f ? 1.f / n_templates : s1_prob_threshold_);
  for (std::size_t i = 0; i < best_s3_haps_.size(); ++i)
  {
    std::size_t uniq_idx = best_s3_haps_[i];
    assert(uniq_idx < s3_reverse_map.size());
    assert(uniq_idx < left_probs.size());
    assert(uniq_idx < left_probs_norecom.size());
    assert(uniq_idx < right_probs.size());
    assert(uniq_idx < right_probs_norecom.size());
    std::size_t cardinality = s3_reverse_map[uniq_idx].size();
    double left_full_prob = left_probs[uniq_idx];
    double left_norecom_prob = left_probs_norecom[uniq_idx];
    double left_recom_prob = (left_full_prob - left_norecom_prob) / cardinality;

    double right_full_prob = right_probs[uniq_idx];
    double right_norecom_prob = right_probs_norecom[uniq_idx];
    double right_recom_prob = (right_full_prob - right_norecom_prob) / cardinality;

    for (std::size_t j = 0; j < cardinality; ++j)
    {
      assert(uniq_idx < s3_reverse_map.size());
      assert(j < s3_reverse_map[uniq_idx].size());
      std::size_t expanded_idx = s3_reverse_map[uniq_idx][j];
      assert(expanded_idx < left_junction_proportions.size());
      assert(expanded_idx < right_junction_proportions.size());
      double orig_prob = (left_norecom_prob * left_junction_proportions[expanded_idx] + left_recom_prob) * (right_norecom_prob * right_junction_proportions[expanded_idx] + right_recom_prob);
      if (orig_prob > denorm_threshold)
      {
        best_s1_probs_.push_back(float(orig_prob / prob_sum));
        best_s1_haps_.push_back(expanded_idx);
      }
    }
  }
}

void hidden_markov_model::s1_to_s2_probs(std::vector<std::size_t>& cardinalities, const std::vector<std::int64_t>& unique_map, std::size_t s2_size)
{
  cardinalities.clear();
  cardinalities.resize(s2_size);
  s2_probs_.clear();
  s2_probs_.resize(s2_size);
  for (std::size_t i = 0; i < best_s1_haps_.size(); ++i)
  {
    std::size_t s2_idx = unique_map[best_s1_haps_[i]];
    ++cardinalities[s2_idx];
    //full_ref_prob_hash[s2_idx] += best_s1_probs_[i];
    s2_probs_[s2_idx] += best_s1_probs_[i];
  }

  best_s2_haps_.clear();
  for (std::size_t i = 0; i < cardinalities.size(); ++i)
  {
    if (cardinalities[i] > 0)
      best_s2_haps_.push_back(i);
  }
}

void hidden_markov_model::impute(double& prob_sum, std::size_t& prev_best_typed_hap,
  const std::vector<float>& left_probs,
  const std::vector<float>& right_probs,
  const std::vector<float>& left_probs_norecom,
  const std::vector<float>& right_probs_norecom,
  const std::vector<float>& left_junction_proportions,
  const std::vector<float>& right_junction_proportions,
  const std::vector<float>& constants,
  const std::vector<std::int64_t>& uniq_map,
  const std::vector<std::vector<std::size_t>>& reverse_map,
  const std::vector<std::int8_t>& template_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t row, std::size_t column, std::size_t out_column,
  full_dosages_results& output,
  reduced_haplotypes::iterator& full_ref_ritr, const reduced_haplotypes::iterator& full_ref_rend,
  std::size_t& prev_block_idx)
{
  assert(row == 0 || tar_variants[row].pos >= tar_variants[row - 1].pos);
  std::size_t mid_point = row == 0 ? 1 : std::max<std::int32_t>(1, std::int32_t(tar_variants[row].pos) - std::int32_t(tar_variants[row].pos - tar_variants[row - 1].pos) / 2);
  if (full_ref_ritr == full_ref_rend || full_ref_ritr->pos < mid_point) // TODO: stop traverse_backward at beginning of region.
    return;

  float typed_dose = std::numeric_limits<float>::quiet_NaN();
  float typed_loo_dose = std::numeric_limits<float>::quiet_NaN();
  std::vector<std::uint32_t> best_typed_haps;
  std::vector<float> best_typed_probs;
  impute_typed_site(prob_sum, prev_best_typed_hap,
    left_probs,
    right_probs,
    left_probs_norecom,
    right_probs_norecom,
    left_junction_proportions,
    right_junction_proportions,
    constants,
    reverse_map,
    template_haps,
    tar_variants[row].gt[column], tar_variants[row].err, tar_variants[row].af, best_typed_haps, best_typed_probs,
    typed_dose, typed_loo_dose);



  //  if (full_ref_ritr == full_ref_rend || full_ref_ritr->pos <= mid_point) // TODO: stop traverse_backward at beginning of region.
  //    return;

  if (best_s3_haps_.empty()
    ||  best_typed_haps != best_s3_haps_
    || std::inner_product(best_typed_probs.begin(), best_typed_probs.end(), best_s3_probs_.begin(), 0.f, std::plus<float>(), [](float l, float r) { return std::abs(l - r); }) > diff_threshold_)
  {
    best_s3_haps_ = best_typed_haps;
    best_s3_probs_ = best_typed_probs;

    s3_to_s1_probs(
      left_probs,
      right_probs,
      left_probs_norecom,
      right_probs_norecom,
      left_junction_proportions,
      right_junction_proportions,
      reverse_map, prob_sum);

    prev_block_idx = std::size_t(-1);
  }

  // vvvvvvvvvvvvvvvv TODO vvvvvvvvvvvvvvvv //
  double best_sum = std::accumulate(best_typed_probs.begin(), best_typed_probs.end(), 0.);
  std::size_t n_templates = left_junction_proportions.size();
  for ( ; full_ref_ritr != full_ref_rend && full_ref_ritr->pos >= mid_point; --full_ref_ritr)
  {
    if (sites_match(tar_variants[row], *full_ref_ritr))
    {
      output.dosage(full_ref_ritr.global_idx(), out_column) = typed_dose;
      output.loo_dosage(row, out_column) = typed_loo_dose;
      if (tar_variants[row].pos == mid_point)
        ++mid_point; // equivalent to breaking and decrementing full_ref_ritr
    }
    else
    {
      if (full_ref_ritr.block_idx() != prev_block_idx)
      {
        s1_to_s2_probs(s2_cardinalities_, full_ref_ritr.unique_map(), full_ref_ritr.cardinalities().size());
        prev_block_idx = full_ref_ritr.block_idx();
      }

      double p_alt = 0;
      std::size_t ac = 0;
      std::size_t an = 0;
      for (std::size_t i = 0; i < best_s2_haps_.size(); ++i)
      {
        std::size_t card = s2_cardinalities_[best_s2_haps_[i]];

        if (full_ref_ritr->gt[best_s2_haps_[i]])
        {
          p_alt += s2_probs_[best_s2_haps_[i]]; //best_s2_probs_[i];
          ac += card;
        }
        an += card;
      }

      if (p_alt > 1. || p_alt < 0.)
      {
        auto a = 0;
      }

      assert(full_ref_ritr->ac >= ac);
      assert(n_templates >= an);
      if (n_templates - an > 0)
        p_alt += (1. - best_sum) * (double(full_ref_ritr->ac - ac) / double(n_templates - an));

      if (p_alt > 1.0f)
      {
        auto a = 0;
      }

      p_alt = std::max(0., std::min(1., p_alt));
      output.dosage(full_ref_ritr.global_idx(), out_column) = float(std::int16_t(p_alt * bin_scalar_ + 0.5f)) / bin_scalar_;
    }
  }
}
