#ifndef MINIMAC4_HIDDEN_MARKOV_MODEL_HPP
#define MINIMAC4_HIDDEN_MARKOV_MODEL_HPP

#include "unique_haplotype.hpp"
#include "variant.hpp"

#include <array>

class best_templates_results
{
private:
  std::vector<std::vector<float>> dosages_;
  std::vector<std::vector<float>> loo_dosages_;
  std::vector<std::vector<std::vector<std::uint32_t>>> best_templates_;
  std::vector<std::vector<std::vector<float>>> best_templates_probs_;
  std::vector<std::vector<std::vector<bool>>> best_templates_uniq_;
public:
  void resize(std::size_t n_rows, std::size_t n_columns)
  {
    dosages_.resize(n_rows, std::vector<float>(n_columns));
    loo_dosages_.resize(n_rows, std::vector<float>(n_columns));
    best_templates_.resize(n_rows, std::vector<std::vector<std::uint32_t>>(n_columns));
    best_templates_probs_.resize(n_rows, std::vector<std::vector<float>>(n_columns));
    best_templates_uniq_.resize(n_rows, std::vector<std::vector<bool>>(n_columns));
  }

  std::array<std::size_t, 2> dimensions() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }

//  void swap(std::size_t row, std::size_t column, std::vector<std::uint32_t>& best_haps, std::vector<float>& best_probs, float dosage, float loo_dosage)
//  {
//    std::swap(best_templates_[row][column], best_haps);
//    std::swap(best_templates_probs_[row][column], best_probs);
//    std::swap(loo_dosages_[row][column], loo_dosage);
//    std::swap(dosages_[row][column], dosage);
//  }

  float& dosage(std::size_t i, std::size_t j) { return dosages_[i][j]; }
  const float& dosage(std::size_t i, std::size_t j) const  { return dosages_[i][j]; }
  float& loo_dosage(std::size_t i, std::size_t j) { return loo_dosages_[i][j]; }
  const float& loo_dosage(std::size_t i, std::size_t j) const { return loo_dosages_[i][j]; }
  std::vector<std::uint32_t>& best_templates(std::size_t i, std::size_t j) { return best_templates_[i][j]; }
  const std::vector<std::uint32_t>& best_templates(std::size_t i, std::size_t j) const { return best_templates_[i][j]; }
  std::vector<float>& best_probs(std::size_t i, std::size_t j) { return best_templates_probs_[i][j]; }
  const std::vector<float>& best_probs(std::size_t i, std::size_t j) const { return best_templates_probs_[i][j]; }
  std::vector<bool>& best_uniq_flags(std::size_t i, std::size_t j) { return best_templates_uniq_[i][j]; }
  const std::vector<bool>& best_uniq_flags(std::size_t i, std::size_t j) const { return best_templates_uniq_[i][j]; }
};

class hidden_markov_model
{
private:
  std::deque<std::vector<std::vector<float>>> forward_probs_;
  std::deque<std::vector<std::vector<float>>> forward_norecom_probs_;
  std::vector<std::vector<float>> junction_prob_proportions_;
  std::vector<bool> precision_jumps_;
  float background_error_;
  static constexpr float jump_fix = 1e15f;
public:
  hidden_markov_model(float background_error = 1e-5f);

  void traverse_forward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx,
    const std::vector<float>& recom,
    const std::vector<float>& err,
    const std::vector<float>& freq);

  void traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx,
    const std::vector<float>& recom,
    const std::vector<float>& err,
    const std::vector<float>& freq,
    const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
    best_templates_results& output);
private:
  void condition(std::vector<float>& probs, std::vector<float>& probs_norecom, const std::vector<std::int8_t>& template_haps, std::int8_t observed, float err, float freq);
  bool transpose(const std::vector<float>& from, std::vector<float>& to, const std::vector<float>& from_norecom, std::vector<float>& to_norecom, const std::vector<std::size_t>& uniq_cardinalities, double recom, std::size_t n_templates);

  void impute(double& prob_sum, std::size_t& prev_best_hap,
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
    best_templates_results& output, std::size_t row, std::size_t column);

  void initialize_likelihoods(std::vector<float>& start_row, const std::vector<std::size_t>& cardinalities);
};

#endif // MINIMAC4_HIDDEN_MARKOV_MODEL_HPP