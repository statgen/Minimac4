#include "hidden_markov_model.hpp"
#include "variant.hpp"

#include <cassert>
#include <numeric>

hidden_markov_model::hidden_markov_model(float background_error) :
  background_error_(background_error)
{
}

void initialize_likelihoods(std::vector<float>& start_row, const std::vector<std::size_t>& cardinalities)
{
  start_row.resize(cardinalities.size());
  for (std::size_t i = 0; i < start_row.size(); ++i)
    start_row[i] = 1.f * cardinalities[i];
}

void hidden_markov_model::traverse_forward(const std::deque<unique_haplotype_block>& ref_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t hap_idx,
  const std::vector<float>& recom,
  const std::vector<float>& err,
  const std::vector<float>& freq)
{
  std::size_t n_expanded_haplotypes = ref_haps.front().expanded_haplotype_size();

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
  std::vector<float> expanded_probs(n_expanded_haplotypes);
  for (std::size_t block_idx = 0; block_idx < ref_haps.size(); ++block_idx,++global_idx)
  {
    const unique_haplotype_block& ref_block = ref_haps[block_idx];
    if (block_idx == 0)
    {
      // Initialize likelihoods at first position
      initialize_likelihoods(forward_probs_[block_idx][0], ref_block.cardinalities());
      forward_norecom_probs_[block_idx][0] = forward_probs_[block_idx][0];
      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
        junction_prob_proportions_[block_idx][i] =  1.f / ref_block.cardinalities()[ref_block.unique_map()[i]];
    }
    else
    {
      // Transition from previous block
      const unique_haplotype_block& prev_ref_block = ref_haps[block_idx - 1];
      std::vector<float>& prev_junction_proportions = junction_prob_proportions_[block_idx - 1];
      std::vector<float>& prev_last_row = forward_probs_[block_idx - 1].back();
      std::vector<float>& prev_last_row_norecom = forward_norecom_probs_[block_idx - 1].back();
      std::vector<float>& first_row = forward_probs_[block_idx].front();
      std::vector<float>& first_row_norecom = forward_norecom_probs_[block_idx].front();
      temp_row.clear();
      temp_row.resize(first_row.size(), 0.f);
      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        std::size_t prev_uniq_idx = prev_ref_block.unique_map()[i];
        float p = prev_last_row_norecom[prev_uniq_idx] * prev_junction_proportions[i] + (prev_last_row[prev_uniq_idx] - prev_last_row_norecom[prev_uniq_idx]) / prev_ref_block.cardinalities()[prev_uniq_idx];
        expanded_probs[i] = p;
        temp_row[uniq_idx] += p;
        //temp_row[uniq_idx] += prev_last_row[prev_uniq_idx] / prev_ref_block.cardinalities()[prev_uniq_idx];
        assert((prev_last_row[prev_uniq_idx] - prev_last_row_norecom[prev_uniq_idx]) >= 0.f);
      }

      std::vector<float>& junction_proportions = junction_prob_proportions_[block_idx];
      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        assert(temp_row[uniq_idx] > 0.f);
        junction_proportions[i] = expanded_probs[i] / temp_row[uniq_idx];
        assert(junction_proportions[i] >= 0.f);
        assert(junction_proportions[i] <= 1.f);
        //temp_row[uniq_idx] += prev_last_row[prev_uniq_idx] / prev_ref_block.cardinalities()[prev_uniq_idx];
      }
#ifndef NDEBUG
      float s1 = std::accumulate(prev_last_row.begin(), prev_last_row.end(), 0.f);
      float s2 = std::accumulate(temp_row.begin(), temp_row.end(), 0.f);
      float sdiff = s2 - s1;
#endif

      transpose(temp_row, first_row, temp_row, first_row_norecom, ref_block.cardinalities(), recom[global_idx], n_expanded_haplotypes);
    }


    const auto& template_haps = ref_block.compressed_matrix();
    std::size_t n_rows = forward_probs_[block_idx].size();
    std::size_t last_row_idx = n_rows - 1;
    for (std::size_t i = 0; i < last_row_idx; ++i,++global_idx)
    {
      std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
      if (observed >= 0)
        condition(forward_probs_[block_idx][i], forward_norecom_probs_[block_idx][i], template_haps[i], observed, err[global_idx], freq[global_idx]);
      transpose(forward_probs_[block_idx][i], forward_probs_[block_idx][i + 1], forward_norecom_probs_[block_idx][i], forward_norecom_probs_[block_idx][i + 1], ref_block.cardinalities(), recom[global_idx], n_expanded_haplotypes);
    }

    std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
    if (observed >= 0)
      condition(forward_probs_[block_idx][last_row_idx], forward_norecom_probs_[block_idx][last_row_idx], template_haps[last_row_idx], observed, err[global_idx], freq[global_idx]);
  }
}

void hidden_markov_model::traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t hap_idx,
  const std::vector<float>& recom,
  const std::vector<float>& err,
  const std::vector<float>& freq,
  const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
  best_templates_results& output)
{
  std::size_t global_idx = tar_variants.size() - 1;
  std::vector<float> backward;
  std::vector<float> backward_norecom;
  std::vector<float> junction_proportions_backward;
  std::vector<float> extra;
  std::vector<float> constants;

  int last_block_idx = int(ref_haps.size()) - 1;
  for (int block_idx = last_block_idx; block_idx >= 0; --block_idx)
  {
    const unique_haplotype_block& ref_block = ref_haps[block_idx];
    std::size_t n_expanded_haplotypes = ref_block.expanded_haplotype_size();
    if (block_idx == last_block_idx)
    {
      // Initialize likelihoods at first position
      initialize_likelihoods(backward, ref_block.cardinalities());
      backward_norecom = backward;
      junction_proportions_backward.resize(ref_block.expanded_haplotype_size());
      for (std::size_t i = 0; i < n_expanded_haplotypes; ++i)
        junction_proportions_backward[i] =  1.f / ref_block.cardinalities()[ref_block.unique_map()[i]];
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
        extra[uniq_idx] += p; //backward[prev_uniq_idx] / prev_ref_block.cardinalities()[prev_uniq_idx];
        assert((backward[prev_uniq_idx] - backward_norecom[prev_uniq_idx]) >= 0.f);
      }

      for (std::size_t i = 0; i < ref_block.expanded_haplotype_size(); ++i)
      {
        std::size_t uniq_idx = ref_block.unique_map()[i];
        assert(extra[uniq_idx] > 0.f);
        junction_proportions_backward[i] /= extra[uniq_idx];
        assert(junction_proportions_backward[i] >= 0.f);
//        assert(junction_proportions_backward[i] <= 1.f);
        if (junction_proportions_backward[i] > 1.f)
        {
          auto a = junction_proportions_backward[i];
          auto a2 = a;
        }
      }

#ifndef NDEBUG
      float s1 = std::accumulate(backward.begin(), backward.end(), 0.f);
      float s2 = std::accumulate(extra.begin(), extra.end(), 0.f);
      float sdiff = s2 - s1;
#endif

      std::swap(backward, extra);
      backward_norecom = backward;
    }

    constants.clear();
    constants.resize(reverse_maps[block_idx].size());
    for (std::size_t i = 0; i < reverse_maps[block_idx].size(); ++i)
    {
      for (std::size_t j = 0; j < reverse_maps[block_idx][i].size(); ++j)
        constants[i] += junction_prob_proportions_[block_idx][reverse_maps[block_idx][i][j]] * junction_proportions_backward[reverse_maps[block_idx][i][j]];
    }

    extra.resize(backward.size());

    const auto& template_haps = ref_block.compressed_matrix();
    std::size_t n_rows = forward_probs_[block_idx].size();
    for (int i = int(n_rows) - 1; i >= 0; --i,--global_idx)
    {
      transpose(backward, backward, backward_norecom, backward_norecom, ref_block.cardinalities(), recom[global_idx], n_expanded_haplotypes);

      std::vector<float>& forward = forward_probs_[block_idx][i];
      assert(forward.size() == backward.size());
      for (int j = 0; j < forward.size(); j++)
        extra[j] = backward[j] * forward[j];

      std::int8_t observed = tar_variants[global_idx].gt[hap_idx];
      impute2(extra, forward_probs_[block_idx][i], backward, forward_norecom_probs_[block_idx][i], backward_norecom, junction_prob_proportions_[block_idx], junction_proportions_backward, constants, reverse_maps[block_idx], template_haps[i], observed, err[global_idx], freq[global_idx], output, global_idx, hap_idx);

      if (observed >= 0)
        condition(backward, backward_norecom, template_haps[i], observed, err[global_idx], freq[global_idx]);
      if (global_idx == 0)
      {
        auto a = 0;
      }
    }
  }
  assert(global_idx == std::size_t(-1));
}


void hidden_markov_model::impute_old(const std::vector<float>& probs, const std::vector<std::int8_t>& template_haps, std::int8_t observed, float err, float af, best_templates_results& output, std::size_t row, std::size_t column)
{
  const float threshold = 0.01;
  float p_alt = 0.f;
  float sum = std::accumulate(probs.begin(), probs.end(), 0.f);
  std::vector<float>& best_probs = output.best_probs(row, column);
  std::vector<std::uint32_t>& best_haps = output.best_templates(row, column);
  float denorm_threshold = threshold * sum;
  for (std::size_t i = 0; i < probs.size(); ++i)
  {
    if (template_haps[i])
      p_alt += probs[i];

    if (probs[i] > denorm_threshold)
    {
      best_probs.push_back(probs[i] / sum);
      best_haps.push_back(i);
      assert(probs[i] / sum > threshold); // TODO: remove since rounding error could make this assert fail.
    }

    //      float normalized_prob = probs[i] / sum;
    //      if (normalized_prob > threshold)
    //      {
    //        best_probs.push_back(normalized_prob);
    //        best_haps.push_back(i);
    //      }
    //
    //      if (template_haps[i])
    //        p_alt += normalized_prob;
  }

  p_alt /= sum;

  output.dosage(row, column) = p_alt;

  if (observed < 0)
  {
    output.loo_dosage(row, column) = std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    float loo_p_alt = p_alt;
    float loo_p_ref = 1.f - p_alt;

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

    output.loo_dosage(row, column) = loo_p_alt / (loo_p_alt + loo_p_ref);
  }
}

void hidden_markov_model::impute2(const std::vector<float>& probs_not_used,
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
  best_templates_results& output, std::size_t row, std::size_t column)
{
  const float threshold = 0.01;
  float p_alt = 0.f;
  std::vector<double> probs(constants.size());
  std::vector<double> probs_broken(constants.size());
  for (std::size_t i = 0; i < constants.size(); ++i)
  {
    double lr = double(left_probs[i]) - double(left_probs_norecom[i]);
    double rr = double(right_probs[i]) - double(right_probs_norecom[i]);
    std::size_t n = reverse_map[i].size();
    double broken = double(constants[i]) * double(left_probs_norecom[i]) * double(right_probs_norecom[i]) + (double(left_probs[i]) * double(right_probs[i]) - double(left_probs_norecom[i]) * double(right_probs_norecom[i])) / reverse_map[i].size();
    double fixed = double(constants[i]) * double(left_probs_norecom[i]) * double(right_probs_norecom[i]) + ((double(left_probs_norecom[i]) * rr + double(right_probs_norecom[i]) * lr + (lr * rr) / n) / n);
    probs[i] = fixed;
    probs_broken[i] = broken;
    if (std::abs(broken - fixed) > 0.f)
    {
      auto a = 0;
    }
    assert(left_probs[i] * right_probs[i] - left_probs_norecom[i] * right_probs_norecom[i] >= 0.f);
  }

  double sum = std::accumulate(probs.begin(), probs.end(), 0.);
  double sum_broken = std::accumulate(probs_broken.begin(), probs_broken.end(), 0.);
  double reduced_prob_sum = std::accumulate(probs_not_used.begin(), probs_not_used.end(), 0.);
  std::vector<std::size_t> best_unique_haps;
  std::vector<double> best_unique_probs;
  double denorm_threshold = threshold * sum;
  for (std::size_t i = 0; i < probs.size(); ++i)
  {
    if (template_haps[i])
      p_alt += probs[i];

    if (probs[i] > denorm_threshold)
    {
      best_unique_probs.push_back(probs[i] / sum);
      best_unique_haps.push_back(i);
      assert(probs[i] / sum > threshold); // TODO: remove since rounding error could make this assert fail.
    }

    //      float normalized_prob = probs[i] / sum;
    //      if (normalized_prob > threshold)
    //      {
    //        best_probs.push_back(normalized_prob);
    //        best_haps.push_back(i);
    //      }
    //
    //      if (template_haps[i])
    //        p_alt += normalized_prob;
  }

  p_alt /= sum;

  output.dosage(row, column) = p_alt;

  if (observed < 0)
  {
    output.loo_dosage(row, column) = std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    float loo_p_alt = p_alt;
    float loo_p_ref = 1.f - p_alt;

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

    output.loo_dosage(row, column) = loo_p_alt / (loo_p_alt + loo_p_ref);
  }

  double best_orig_prob = 0.f;
  std::vector<float>& best_probs = output.best_probs(row, column);
  std::vector<std::uint32_t>& best_haps = output.best_templates(row, column);
  std::vector<bool>& best_uniq_flags = output.best_uniq_flags(row, column);
  for (std::size_t i = 0; i < best_unique_haps.size(); ++i)
  {
    std::size_t uniq_idx = best_unique_haps[i];
    std::size_t cardinality = reverse_map[uniq_idx].size();
    double left_full_prob = left_probs[uniq_idx];
    double left_norecom_prob = left_probs_norecom[uniq_idx];
    double left_recom_prob = (left_full_prob - left_norecom_prob) / cardinality;

    double right_full_prob = right_probs[uniq_idx];
    double right_norecom_prob = right_probs_norecom[uniq_idx];
    double right_recom_prob = (right_full_prob - right_norecom_prob) / cardinality;

    std::size_t old_size = best_haps.size();
    bool found = false;
    double local_sum = 0.f;
    double local_sum_used = 0.f;
    for (std::size_t j = 0; j < cardinality; ++j)
    { assert(uniq_idx < reverse_map.size()); assert(j < reverse_map[uniq_idx].size());
      std::size_t expanded_idx = reverse_map[uniq_idx][j];
      double orig_prob = (left_norecom_prob * left_junction_proportions[expanded_idx] + left_recom_prob) * (right_norecom_prob * right_junction_proportions[expanded_idx] + right_recom_prob);
      if (orig_prob >= denorm_threshold)
      {
        best_haps.push_back(expanded_idx);
        best_probs.push_back(orig_prob / sum);
        best_uniq_flags.push_back(false);
        found = true;
        local_sum_used += orig_prob / sum;
      }
      if (orig_prob > best_orig_prob)
        best_orig_prob = orig_prob;

      local_sum += orig_prob / sum;
    }
    //local_sum /= sum;
    auto diff = best_unique_probs[i] - local_sum;
    auto a = 0;

    //double lost_prob = best_unique_probs[i] - local_sum_used;
    if (/*lost_prob > 0.05 || */!found)
    {
#if 1
      best_haps.push_back(best_unique_haps[i]);
      best_probs.push_back(best_unique_probs[i]);
      best_uniq_flags.push_back(true);
#else
      best_haps.resize(old_size);
      best_probs.resize(old_size);
      best_uniq_flags.resize(old_size);
      for (std::size_t j = 0; j < cardinality; ++j)
      {
        std::size_t expanded_idx = reverse_map[uniq_idx][j];
        double orig_prob = (left_norecom_prob * left_junction_proportions[expanded_idx] + left_recom_prob) * (right_norecom_prob * right_junction_proportions[expanded_idx] + right_recom_prob);

        best_haps.push_back(expanded_idx);
        best_probs.push_back(orig_prob / sum);
        best_uniq_flags.push_back(false);
      }
#endif
    }

  }
  assert(best_haps.size() > 0);

#ifndef NDEBUG
  double expanded_sum = 0.f;
  for (std::size_t i = 0; i < reverse_map.size(); ++i)
  {
    std::size_t uniq_idx = i;
    std::size_t cardinality = reverse_map[uniq_idx].size();

    double left_full_prob = left_probs[uniq_idx];
    double left_norecom_prob = left_probs_norecom[uniq_idx];
    double left_recom_prob = (left_full_prob - left_norecom_prob) / cardinality;

    double right_full_prob = right_probs[uniq_idx];
    double right_norecom_prob = right_probs_norecom[uniq_idx];
    double right_recom_prob = (right_full_prob - right_norecom_prob) / cardinality;

    for (std::size_t j = 0; j < cardinality; ++j)
    {
      std::size_t expanded_idx = reverse_map[uniq_idx][j];
      double expanded_prob = (left_norecom_prob * left_junction_proportions[expanded_idx] + left_recom_prob) * (right_norecom_prob * right_junction_proportions[expanded_idx] + right_recom_prob);
      expanded_sum += expanded_prob;
    }
  }

  double sum_diff = expanded_sum - sum;
  double sum_diff_broken = expanded_sum - sum_broken;
  auto fix_broken_diff = std::abs(sum_diff_broken) - std::abs(sum_diff);
  auto a = 0;
#endif
}

void hidden_markov_model::impute(const std::vector<float>& probs,
  const std::vector<float>& left_probs,
  const std::vector<float>& right_probs,
  const std::vector<float>& left_probs_norecom,
  const std::vector<float>& right_probs_norecom,
  const std::vector<float>& left_junction_proportions,
  const std::vector<float>& right_junction_proportions,
  const std::vector<std::vector<std::size_t>>& reverse_map,
  const std::vector<std::int8_t>& template_haps,
  std::int8_t observed, float err, float af,
  best_templates_results& output, std::size_t row, std::size_t column)
{
  const float threshold = 0.01;
  float p_alt = 0.f;
  float sum = std::accumulate(probs.begin(), probs.end(), 0.f);
  std::vector<float>& best_probs = output.best_probs(row, column);
  std::vector<std::uint32_t>& best_haps = output.best_templates(row, column);
  float denorm_threshold = threshold * sum;
  for (std::size_t i = 0; i < probs.size(); ++i)
  {
    if (template_haps[i])
      p_alt += probs[i];

    if (probs[i] > denorm_threshold)
    {
      best_probs.push_back(probs[i] / sum);
      best_haps.push_back(i);
      assert(probs[i] / sum > threshold); // TODO: remove since rounding error could make this assert fail.
    }

    //      float normalized_prob = probs[i] / sum;
    //      if (normalized_prob > threshold)
    //      {
    //        best_probs.push_back(normalized_prob);
    //        best_haps.push_back(i);
    //      }
    //
    //      if (template_haps[i])
    //        p_alt += normalized_prob;
  }

  p_alt /= sum;

  output.dosage(row, column) = p_alt;

  if (observed < 0)
  {
    output.loo_dosage(row, column) = std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    float loo_p_alt = p_alt;
    float loo_p_ref = 1.f - p_alt;

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

    output.loo_dosage(row, column) = loo_p_alt / (loo_p_alt + loo_p_ref);
  }

  float global_sum = 0.f;
  for (std::size_t i = 0; i < best_haps.size(); ++i)
  {
    std::size_t uniq_idx = best_haps[i];
    float local_sum = 0.f;
    std::size_t cardinality = reverse_map[uniq_idx].size();
    float left_full_prob = left_probs[uniq_idx];
    float left_norecom_prob = left_probs_norecom[uniq_idx];
    float left_recom_prob = (left_full_prob - left_norecom_prob) / cardinality;

    float right_full_prob = right_probs[uniq_idx];
    float right_norecom_prob = right_probs_norecom[uniq_idx];
    float right_recom_prob = (right_full_prob - right_norecom_prob) / cardinality;
    float l = 0.f, r = 0.f;
    for (std::size_t j = 0; j < cardinality; ++j)
    {
      float orig_prob = (left_norecom_prob * left_junction_proportions[reverse_map[uniq_idx][j]] + left_recom_prob) * (right_norecom_prob * right_junction_proportions[reverse_map[uniq_idx][j]] + right_recom_prob);
      l += (left_norecom_prob * left_junction_proportions[reverse_map[uniq_idx][j]] + left_recom_prob);
      r += (right_norecom_prob * right_junction_proportions[reverse_map[uniq_idx][j]] + right_recom_prob);
      local_sum += orig_prob;
      if (orig_prob / sum > 0.01)
      {
        auto a = 0;
      }
    }
    local_sum /= sum;
    float l_r = l * r;
    l_r /= sum;
    auto diff = best_probs[i] - local_sum;
    global_sum += local_sum;
  }
  auto b = 0;
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

void hidden_markov_model::transpose(const std::vector<float>& from, std::vector<float>& to, const std::vector<float>& from_norecom, std::vector<float>& to_norecom, const std::vector<std::size_t>& uniq_cardinalities, double recom, std::size_t n_templates)
{
  assert(from.size() == to.size());

  double sum = 0.;

  for (std::size_t i = 0; i < from.size(); ++i)
  {
    sum += from[i];
  }

  //sum *= (recom / n_templates);
  double complement = 1. - recom;

  // avoid underflows
  if (sum < 1e-10)
  {
    sum *= 1e15;
    complement *= 1e15;
//    for(int i=0;i<noReducedStatesCurrent;i++)
//      noRecomProb[i]*=JumpFix;
  }

  sum *= (recom / n_templates);


  for (int i = 0; i < to.size(); i++)
  {
    to[i] = from[i] * complement + (uniq_cardinalities[i] * sum);
    to_norecom[i] = from_norecom[i] * complement;
    assert(to[i]>=0.0f);
    //assert(noRecomProb[i]>=0.0f);
    assert(to[i]<1e18);
    //assert(noRecomProb[i]<1e18);
  }
}
