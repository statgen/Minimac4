#ifndef GETOPT_WRAPPER_HPP
#define GETOPT_WRAPPER_HPP

#include <cstdint>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <getopt.h>

class getopt_wrapper
{
public:
  struct option_with_desc : public ::option
  {
    const char* description;
    option_with_desc(const char* _name, int _has_arg, int* _flag, int _val, const char* _description)
    {
      name = _name;
      has_arg = _has_arg;
      flag = _flag;
      val = _val;
      description = _description;
    }
  };
protected:
  std::vector<option> long_options_;
  std::vector<option_with_desc> opts_;
  std::string usage_str_;
  std::string short_opt_string_;
  std::size_t max_long_opt_length_ = 0;
public:
  getopt_wrapper(std::string usage_str, std::vector<option_with_desc> long_opts) :
    usage_str_(std::move(usage_str)),
    opts_(std::move(long_opts))
  {

    long_options_.resize(opts_.size() + 1, {0, 0, 0, 0});
    auto lit = long_options_.begin();
    for (auto it = opts_.begin(); it != opts_.end(); ++it)
    {
      *(lit++) = *it;
      max_long_opt_length_ = std::max(max_long_opt_length_, it->name ? std::strlen(it->name) : 0);
      if (it->val)
      {
        short_opt_string_ += (char)it->val;
        if (it->has_arg == required_argument)
          short_opt_string_ += ':';
      }
    }
  }

  void print_usage(std::ostream& os)
  {
    os << usage_str_ << '\n';
    os << '\n';
    for (auto it = opts_.begin(); it != opts_.end(); ++it)
    {
      if (!it->description)
        continue;

      if (std::isprint(it->val))
      {
        if (it->name)
          os << " -" << (char)it->val << ", ";
        else
          os << " -" << (char)it->val << "  ";
      }
      else
        os << "     ";

      std::size_t n_spaces = 2;
      if (it->name)
        os << "--" << it->name;
      else
        n_spaces += 2;

      n_spaces += max_long_opt_length_ - std::strlen(it->name);
      for (std::size_t i = 0; i < n_spaces; ++i)
        os.put(' ');
      os << it->description << '\n';
    }

    os << std::flush;
  }
};

#endif // GETOPT_WRAPPER_HPP
