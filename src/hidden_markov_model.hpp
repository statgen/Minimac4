#ifndef MINIMAC4_HIDDEN_MARKOV_MODEL_HPP
#define MINIMAC4_HIDDEN_MARKOV_MODEL_HPP

#include "unique_haplotype.hpp"
#include "variant.hpp"

#include <array>
#include <cassert>
#include <numeric>

class full_dosages_results
{
public:
  std::vector<std::vector<float>> dosages_;
  std::vector<std::vector<float>> loo_dosages_;
public:
  void resize(std::size_t n_rows, std::size_t n_loo_rows, std::size_t n_columns)
  {
    dosages_.resize(n_rows, std::vector<float>(n_columns, savvy::typed_value::end_of_vector_value<float>()));
    loo_dosages_.resize(n_loo_rows, std::vector<float>(n_columns, savvy::typed_value::end_of_vector_value<float>()));
  }

  std::array<std::size_t, 2> dimensions() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }
  std::array<std::size_t, 2> dimensions_loo() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }

  float& dosage(std::size_t i, std::size_t j) { return dosages_[i][j]; }
  const float& dosage(std::size_t i, std::size_t j) const  { return dosages_[i][j]; }
  float& loo_dosage(std::size_t i, std::size_t j) { return loo_dosages_[i][j]; }
  const float& loo_dosage(std::size_t i, std::size_t j) const { return loo_dosages_[i][j]; }
};

class hidden_markov_model
{
private:
  std::deque<std::vector<std::vector<float>>> forward_probs_;
  std::deque<std::vector<std::vector<float>>> forward_norecom_probs_;
  std::vector<std::vector<float>> junction_prob_proportions_;
  std::vector<bool> precision_jumps_;
  float prob_threshold_ = 0.01f;
  float s1_prob_threshold_ = -1.f;
  float diff_threshold_ = 0.01f;
  float background_error_ = 1e-5f;
  static constexpr float jump_fix = 1e15f;
  static constexpr float jump_threshold = 1e-10f;
  const std::int16_t bin_scalar_ = 1000;

  std::vector<std::uint32_t> best_s1_haps_;
  std::vector<std::uint32_t> best_s2_haps_;
  std::vector<std::uint32_t> best_s3_haps_;
  std::vector<float> best_s1_probs_;
  std::vector<float> best_s2_probs_;
  std::vector<float> s2_probs_;
  std::vector<std::size_t> s2_cardinalities_;
  std::vector<float> best_s3_probs_;

public:
  hidden_markov_model(float s3_prob_threshold, float s1_prob_threshold, float diff_threshold, float background_error);

  void traverse_forward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx);

  void traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx,
    std::size_t out_idx,
    const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
    full_dosages_results& output,
    const reduced_haplotypes& full_reference_data);
private:
  void condition(std::vector<float>& probs, std::vector<float>& probs_norecom, const std::vector<std::int8_t>& template_haps, std::int8_t observed, float err, float freq);
  bool transpose(const std::vector<float>& from, std::vector<float>& to, const std::vector<float>& from_norecom, std::vector<float>& to_norecom, const std::vector<std::size_t>& uniq_cardinalities, double recom, std::size_t n_templates);

  void impute_typed_site(double& prob_sum, std::size_t& prev_best_hap,
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
    std::vector<std::uint32_t>& best_uniq_haps, std::vector<float>& best_uniq_probs, float& dose, float& loo_dose);

  void impute(double& prob_sum, std::size_t& prev_best_expanded_hap,
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
    std::size_t& prev_block_idx);

  void initialize_likelihoods(std::vector<float>& probs, std::vector<float>& probs_norecom, std::vector<float>& proportions, const unique_haplotype_block& ref_block);

  void s3_to_s1_probs(
    const std::vector<float>& left_probs, const std::vector<float>& right_probs,
    const std::vector<float>& left_probs_norecom, const std::vector<float>& right_probs_norecom,
    const std::vector<float>& left_junction_proportions, const std::vector<float>& right_junction_proportions,
    const std::vector<std::vector<std::size_t>>& s3_reverse_map, double prob_sum);
  void s1_to_s2_probs(std::vector<std::size_t>& cardinalities, const std::vector<std::int64_t>& uniq_map, std::size_t s2_size);
};

#endif // MINIMAC4_HIDDEN_MARKOV_MODEL_HPP