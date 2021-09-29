#ifndef MINIMAC4_HIDDEN_MARKOV_MODEL_HPP
#define MINIMAC4_HIDDEN_MARKOV_MODEL_HPP

#include "unique_haplotype.hpp"
#include "variant.hpp"

#include <array>
#include <cassert>
#include <numeric>

class full_dosages_results
{
private:
  std::vector<std::vector<float>> dosages_;
  std::vector<std::vector<float>> loo_dosages_;
public:
  void resize(std::size_t n_rows, std::size_t n_columns)
  {
    dosages_.resize(n_rows, std::vector<float>(n_columns));
    loo_dosages_.resize(n_rows, std::vector<float>(n_columns));
  }

  std::array<std::size_t, 2> dimensions() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }

  float& dosage(std::size_t i, std::size_t j) { return dosages_[i][j]; }
  const float& dosage(std::size_t i, std::size_t j) const  { return dosages_[i][j]; }
  float& loo_dosage(std::size_t i, std::size_t j) { return loo_dosages_[i][j]; }
  const float& loo_dosage(std::size_t i, std::size_t j) const { return loo_dosages_[i][j]; }
};

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
    std::size_t hap_idx);

  template <typename ... Args>
  void traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx,
    const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
    Args& ... args); //best_templates_results& output);
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
    std::vector<std::uint32_t>& best_uniq_haps, std::vector<double>& best_uniq_probs, float& dose, float& loo_dose);

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
    const std::vector<target_variant>& tar_variants,
    std::size_t row, std::size_t column, best_templates_results& output);

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
    const std::vector<target_variant>& tar_variants,
    std::size_t row, std::size_t column, full_dosages_results& output, reduced_haplotypes::iterator& full_ref_ritr, const reduced_haplotypes::iterator& full_ref_rend);

  void initialize_likelihoods(std::vector<float>& probs, std::vector<float>& probs_norecom, std::vector<float>& proportions, const unique_haplotype_block& ref_block);
  void typed_to_full_probs(
    std::vector<std::uint32_t>& full_haps,
    std::vector<float>& full_probs,
    const std::vector<std::uint32_t>& typed_haps,
    const std::vector<float>& left_probs, const std::vector<float>& right_probs,
    const std::vector<float>& left_probs_norecom, const std::vector<float>& right_probs_norecom,
    const std::vector<float>& left_junction_proportions, const std::vector<float>& right_junction_proportions,
    const std::vector<std::vector<std::size_t>>& typed_reverse_map,
    const std::vector<std::size_t>& full_map,
    double prob_sum,
    std::size_t full_hap_count);
};

template <typename ... Args>
void hidden_markov_model::traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
  const std::vector<target_variant>& tar_variants,
  std::size_t hap_idx,
  const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
  Args& ... args)
{
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
//      backward_norecom = backward;
//      junction_proportions_backward.resize(ref_block.expanded_haplotype_size());
//      for (std::size_t i = 0; i < n_expanded_haplotypes; ++i)
//        junction_proportions_backward[i] =  1.f / ref_block.cardinalities()[ref_block.unique_map()[i]];
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
      impute(prob_sum, best_hap, forward_probs_[block_idx][i], backward, forward_norecom_probs_[block_idx][i], backward_norecom, junction_prob_proportions_[block_idx], junction_proportions_backward, constants, reverse_maps[block_idx], template_variants[i].gt, tar_variants, global_idx, hap_idx, args...);

      if (observed >= 0)
        condition(backward, backward_norecom, template_variants[i].gt, observed, tar_variants[global_idx].err, tar_variants[global_idx].af);
    }
  }
  assert(global_idx == std::size_t(-1));
}

#endif // MINIMAC4_HIDDEN_MARKOV_MODEL_HPP