
#include "getopt_wrapper.hpp"
#include "unique_haplotype.hpp"
#include "hidden_markov_model.hpp"
#include "recombination.hpp"

#include "MyVariables.h"
#include "Parameters.h"
#include "StringBasics.h"
#include "Analysis.h"

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>
#include <omp.hpp>

#include <stdio.h>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <fstream>
#include <ostream>


#if 0
void Minimac4Version();
void helpFile();

int main_old(int argc, char ** argv)
{
	// Parameter Options

    String refHaps = "";
	String haps = "";
	String recFile = "", errFile = "";
    bool log = false, help = false, params = false;

    AllVariable MyAllVariable;
    OutputFormatVariable &MyOutFormat=MyAllVariable.myOutFormat;
    ModelVariable &MyModelVariables=MyAllVariable.myModelVariables;
    HaplotypeDataVariables &MyHapDataVariables=MyAllVariable.myHapDataVariables;



	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Reference Haplotypes")
		LONG_STRINGPARAMETER("refHaps", &refHaps)
		LONG_PARAMETER("passOnly", &MyHapDataVariables.passOnly)
		LONG_PARAMETER("rsid", &MyOutFormat.RsId)
		LONG_PARAMETER("referenceEstimates", &MyModelVariables.referenceEstimates)
		LONG_STRINGPARAMETER("mapFile", &MyHapDataVariables.mapFile)



		LONG_PARAMETER_GROUP("Target Haplotypes")
		LONG_STRINGPARAMETER("haps", &haps)


		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &MyOutFormat.OutPrefix)
		LONG_PARAMETER("estimate" , &MyModelVariables.processReference)
		LONG_PARAMETER("nobgzip", &MyOutFormat.nobgzip)
		LONG_INTPARAMETER("vcfBuffer", &MyOutFormat.vcfBuffer)
		LONG_STRINGPARAMETER("format", &MyOutFormat.formatString)
		LONG_PARAMETER("allTypedSites", &MyOutFormat.TypedOnly)
		LONG_PARAMETER("meta", &MyOutFormat.meta)
		LONG_PARAMETER("memUsage", &MyOutFormat.memUsage)


		LONG_PARAMETER_GROUP("Chunking Parameters")
		LONG_DOUBLEPARAMETER("ChunkLengthMb", &MyHapDataVariables.ChunkLengthMb)
		LONG_DOUBLEPARAMETER("ChunkOverlapMb", &MyHapDataVariables.ChunkOverlapMb)


		LONG_PARAMETER_GROUP("Subset Parameters")
		LONG_STRINGPARAMETER("chr", &MyHapDataVariables.chr)
		LONG_INTPARAMETER("start", &MyHapDataVariables.start)
		LONG_INTPARAMETER("end", &MyHapDataVariables.end)
		LONG_INTPARAMETER("window", &MyHapDataVariables.window)
		//LONG_INTPARAMETER("block", &max_block)



//		LONG_PARAMETER_GROUP("Starting Parameters")
//		LONG_STRINGPARAMETER("rec", &recFile)
//		LONG_STRINGPARAMETER("err", &errFile)

		LONG_PARAMETER_GROUP("Approximation Parameters")
//		LONG_INTPARAMETER("rounds", &MyModelVariables.rounds)
//		LONG_INTPARAMETER("states", &MyModelVariables.states)
		LONG_PARAMETER("minimac3", &MyModelVariables.minimac3)
		LONG_DOUBLEPARAMETER("probThreshold", &MyModelVariables.probThreshold)
		LONG_DOUBLEPARAMETER("diffThreshold", &MyModelVariables.diffThreshold)
		LONG_DOUBLEPARAMETER("topThreshold", &MyModelVariables.topThreshold)



		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_PARAMETER("log", &log)
//		LONG_PARAMETER("lowMemory", &MyModelVariables.lowMemory)
		LONG_PARAMETER("help", &help)
		LONG_INTPARAMETER("cpus", &MyModelVariables.cpus)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)



		BEGIN_LEGACY_PARAMETERS()
		LONG_STRINGPARAMETER("intermediate", &MyModelVariables.intermediate)
		LONG_PARAMETER("reEstimate", &MyModelVariables.reEstimate)
		LONG_DOUBLEPARAMETER("constantParam", &MyModelVariables.constantParam)
		LONG_INTPARAMETER("printBuffer", &MyOutFormat.PrintBuffer)
		LONG_PARAMETER("verbose", &MyOutFormat.verbose)
		LONG_DOUBLEPARAMETER("minRatioPercent", &MyHapDataVariables.minRatio)
		LONG_STRINGPARAMETER("MyChromosome", &MyHapDataVariables.MyChromosome)
		LONG_INTPARAMETER("ignoreDuplicates", &MyHapDataVariables.ignoreDuplicates)
		LONG_INTPARAMETER("transFactor", &MyModelVariables.transFactor)
		LONG_INTPARAMETER("cisFactor", &MyModelVariables.cisFactor)
		LONG_PARAMETER("allowRefDuplicates", &MyHapDataVariables.allowRefDuplicates)
		LONG_PARAMETER("unphasedOutput", &MyOutFormat.unphasedOutput)
		LONG_PARAMETER("longZero", &MyOutFormat.longZero)
		END_LONG_PARAMETERS();

	//MyHapDataVariables.GetMapFileLocation(argc,argv);
    MyOutFormat.CreateCommandLine(argc,argv);

	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));

    String compStatus;
	inputParameters.Read(argc, &(argv[0]));

    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(MyOutFormat.OutPrefix+".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));


    Minimac4Version();
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();

	int start_time = time(0);


    Analysis MyAnalysis;
    String MySuccessStatus="Error";


//
    MySuccessStatus = MyAnalysis.AnalyzeExperiment(refHaps, haps, recFile, errFile, MyAllVariable);

    if(MySuccessStatus!="Success")
    {
        compStatus=MySuccessStatus;
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }

	int time_tot = time(0) - start_time;

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using Minimac4 !!! "<<endl<<endl;

    if(log)
        fclose (LogFile);


    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void Minimac4Version()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          Minimac4 - Fast Imputation Based on State Space Reduction HMM\n");
	printf(" --------------------------------------------------------------------------------\n");
    printf("           (c) 2014 - Sayantan Das, Christian Fuchsberger, David Hinds\n");
    printf("                             Mary Kate Wing, Goncalo Abecasis \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
    printf("\n URL = http://genome.sph.umich.edu/wiki/Minimac4\n");



}

void helpFile()
{

    printf("\n\n\t  Minimac4 is a lower memory and more computationally efficient implementation of \"minimac2/3\".\n");


    printf("\t It is an algorithm for genotypic imputation that works on phased genotypes (say from MaCH).\n");
    printf("\t Minimac4 is designed to handle very large reference panels in a more computationally efficient \n");
    printf("\t way with no loss of accuracy. This algorithm analyzes only the unique sets of haplotypes in \n");
    printf("\t small genomic segments, thereby saving on time-complexity, computational memory but no loss\n");
    printf("\t in degree of accuracy.\n");

printf("\n\n ----------------------------------------------------------------------------------------- \n");
	printf("                            Minimac4 - List of Usage Options \n");
	printf(" -----------------------------------------------------------------------------------------\n\n");

    printf(" --------- Reference Haplotypes --------- \n");
  printf("\n              --refHaps   : M3VCF file containing haplotype data for reference panel.\n");
    printf("             --passOnly   : This option only imports variants with FILTER = PASS.\n");
    printf("                 --rsid   : This option only imports RS ID of variants from ID column (if available).\n");



  printf("\n --------- GWAS Haplotypes --------- \n");
  printf("\n                 --haps   : File containing haplotype data for target (gwas) samples. Must be VCF \n");
    printf("                            file. Zipped versions allowed.\n");

  printf("\n --------- Output Parameters --------- \n");
  printf("\n               --prefix   : Prefix for all output files generated. By default: [Minimac4.Output]\n");
    printf("     --processReference   : This option will only convert an input VCF file to M3VCF format\n");
    printf("                            (currently de-activated in minimac4). If this option is ON, \n");
    printf("                            no imputation would be performed.\n");
    printf("              --nobgzip   : If ON, output files will NOT be gzipped.\n");
    printf("               --format   : Specifies which fields to output for the FORMAT field in output \n");
    printf("                            VCF file. Available handles: GT,DS,HDS,GP [Default: GT,DS].\n");
    printf("        --allTypedSites   : If ON, sites available ONLY in GWAS panel will also be output [Default: OFF]. \n");



  printf("\n --------- Subset Parameters --------- \n");
  printf("\n                  --chr   : Chromosome number for which to carry out imputation.\n");
    printf("                --start   : Start position for imputation by chunking.\n");
    printf("                  --end   : End position for imputation by chunking. \n");
    printf("               --window   : Length of buffer region on either side of --start and --end.\n");

  printf("\n --------- Other Parameters --------- \n");
  printf("\n                  --log   : If ON, log will be written to $prefix.logfile.\n");
    printf("                 --help   : If ON, detailed help on options and usage.\n");
    printf("                 --cpus   : Number of cpus for parallel computing. Works only with Minimac4-omp.\n\n");


  printf("\n Please visit <http://genome.sph.umich.edu/wiki/Minimac4> for detailed documentation ...\n\n");
    cout<<endl;

	return;
}
#endif
class prog_args : public getopt_wrapper
{
private:
  std::string ref_path_;
  std::string tar_path_;
  std::string map_path_;
  std::string out_path_ = "/dev/stdout";
  savvy::file::format out_format_ = savvy::file::format::bcf;
  std::uint8_t out_compression_ = 6;
  std::vector<std::string> fmt_fields_ = {"GT","HDS","DS"};
  savvy::genomic_region reg_ = {""};
  std::int64_t overlap_;
  std::int16_t threads_ = 1;
  bool deferred_interpolation_ = false;
  bool help_ = false;

public:
  bool help_is_set() const
  {
    return help_;
  }

  const std::string& ref_path() const { return ref_path_; }
  const std::string& tar_path() const { return tar_path_; }
  const std::string& map_path() const { return map_path_; }
  const std::string& out_path() const { return out_path_; }
  savvy::file::format out_format() const { return out_format_; }
  std::uint8_t out_compression() const { return out_compression_; }
  const std::vector<std::string>& fmt_fields() const { return fmt_fields_; }
  const savvy::genomic_region& region() const { return reg_; }
  std::int64_t overlap() const { return overlap_; }
  std::int16_t threads() const { return threads_; }
  bool deferred_interpolation() const { return deferred_interpolation_; }

  prog_args() :
    getopt_wrapper(
      "Usage: minimac4 [opts ...] <reference.m3vcf.gz> <target.vcf.gz>",
      {
        {"help", no_argument, 0, 'h', "Print usage"},
        {"map", required_argument, 0, 'm', "Genetic map file"},
        {"output", required_argument, 0, 'o', "Output path (default: /dev/stdout)"},
        {"output-format", required_argument, 0, 'O', "Output file format (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: bcf)"},
        {"region", required_argument, 0, 'r', "Genomic region to impute"},
        {"threads", required_argument, 0, 't', "Number of threads (default: 1)"},
        {"overlap", required_argument, 0, 'w', "Size (in Bp) of overlap before and after region to use as input to HMM (default: 1000000)"},
        {"deferred-interpolation", no_argument, 0, '\x01', "Enables experimental deferred interpolation algorithm"},
      })
  {
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, short_opt_string_.c_str(), long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      case 'm':
        map_path_ = optarg ? optarg : "";
        break;
      case 'o':
        out_path_ = optarg ? optarg : "";
        break;
      case 'O':
        {
          using fmt = savvy::file::format;
          std::string ot = optarg ? optarg : "";
          if (ot == "vcf")
          {
            out_format_ = fmt::vcf;
            out_compression_ = 0;
          }
          else if (ot == "vcf.gz")
          {
            out_format_ = fmt::vcf;
          }
          else if (ot == "bcf")
          {
            out_format_ = fmt::bcf;
          }
          else if (ot == "ubcf")
          {
            out_format_ = fmt::bcf;
            out_compression_ = 0;
          }
          else if (ot == "sav")
          {
            out_format_ = fmt::sav;
          }
          else if (ot == "usav")
          {
            out_format_ = fmt::sav;
            out_compression_ = 0;
          }
          else
          {
            std::cerr << "Invalid --output-format: " << ot << std::endl;
            return false;
          }
          break;
        }
      case 'r':
        reg_ = string_to_region(optarg ? optarg : "");
        break;
      case 't':
        threads_ = atoi(optarg ? optarg : "");
        break;
      case 'w':
        overlap_ = std::atoll(optarg ? optarg : "");
        break;
      case '\x01':
        if (std::string(long_options_[long_index].name) =="deferred-interpolation")
        {
          deferred_interpolation_ = true;
          break;
        } // else pass through to default
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 2)
    {
      ref_path_ = argv[optind];
      tar_path_ = argv[optind + 1];
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
private:
  savvy::genomic_region string_to_region(const std::string& s)
  {
    const std::size_t colon_pos = s.find(':');
    if (colon_pos == std::string::npos)
    {
      return savvy::genomic_region(s);
    }
    else
    {
      std::string chr = s.substr(0, colon_pos);
      const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
      if (hyphen_pos == std::string::npos)
      {
        std::string slocus = s.substr(colon_pos + 1);
        std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
        return savvy::genomic_region(chr, ilocus, ilocus);
      }
      else
      {
        std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
        std::string send = s.substr(hyphen_pos + 1);
        if (send.empty())
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
        }
        else
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
        }
      }
    }
  }
};

bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids)
{
  savvy::reader input(file_path);
  sample_ids = input.samples();
  input.reset_bounds(reg, savvy::bounding_point::any);
  savvy::variant var;
  std::vector<std::int8_t> tmp_geno;
  while (input >> var)
  {
    var.get_format("GT", tmp_geno);
    for (std::size_t i = 0; i < var.alts().size(); ++i)
    {
      std::size_t allele_idx = i + 1;
      target_sites.push_back({var.chromosome(), var.position(), var.ref(), var.alts()[i], true, false, std::numeric_limits<float>::quiet_NaN(), 0.00999, recombination::recom_min, 0, {}});
      if (var.alts().size() == 1)
        tmp_geno.swap(target_sites.back().gt);
      else
      {
        target_sites.back().gt.resize(tmp_geno.size());
        std::int8_t* dest_ptr = target_sites.back().gt.data();
        for (std::size_t j = 0; j < tmp_geno.size(); ++j)
          dest_ptr[j] = std::int8_t(tmp_geno[j] == allele_idx);
      }
    }
  }

  return !input.bad();
}

bool load_reference_haplotypes(const std::string& file_path, const savvy::genomic_region& extended_reg, const savvy::genomic_region& impute_reg, std::vector<target_variant>& target_sites, reduced_haplotypes& typed_only_reference_data, reduced_haplotypes* full_reference_data)
{
  savvy::reader input(file_path);

  if (input)
  {
    input.reset_bounds(extended_reg);
    savvy::variant var;
    std::vector<std::int8_t> tmp_geno;

    if (full_reference_data)
    {
      auto tar_it = target_sites.begin();
      unique_haplotype_block block;
      std::size_t ref_cnt = 0;
      while (block.deserialize(input))
      {
        if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
          break;

        for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
        {
          while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
          {
            if (ref_cnt)
            {
              tar_it->ref_cnt = ref_cnt;
              ref_cnt = 0;
            }
            ++tar_it;
          }

          for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
          {
            if (it->ref == ref_it->ref && it->alt == ref_it->alt)
            {
              tmp_geno.resize(block.unique_map().size());
              for (std::size_t i = 0; i < tmp_geno.size(); ++i)
                tmp_geno[i] = ref_it->gt[block.unique_map()[i]];
              typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
              it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
              it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
              it->in_ref = true;
              if (it != tar_it)
                std::swap(*it, *tar_it);
              if (ref_cnt)
              {
                tar_it->ref_cnt = ref_cnt;
                ref_cnt = 0;
              }
              ++tar_it;
              break;
            }
          }

          if (ref_it->pos >= extended_reg.from())
            ++ref_cnt;
        }

        if (full_reference_data)
        {
          // TODO: remove redundant variants
          block.trim(impute_reg.from(), impute_reg.to());
          if (!block.variants().empty())
            full_reference_data->append_block(block);
        }
      }
    }
    else
    {
      auto tar_it = target_sites.begin();
      while (input >> var)
      {
        while (tar_it != target_sites.end() && tar_it->pos < var.position())
          ++tar_it;

        for (auto it = tar_it; it != target_sites.end() && it->pos == var.position(); ++it)
        {
          if (it->ref == var.ref() && it->alt == (var.alts().size() ? var.alts()[0] : ""))
          {
            var.get_format("GT", tmp_geno);
            typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
            // freq.push_back(std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size());
            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size(); // TODO; remove
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;
            if (it != tar_it)
              std::swap(*it, *tar_it);
            ++tar_it;
            break;
          }
        }
      }
    }

    return !input.bad();
  }
  else
  {
    shrinkwrap::gz::istream input_file(file_path);
    std::string line;

    std::uint8_t m3vcf_version = 0;
    const std::string m3vcf_version_line = "##fileformat=M3VCF";
    while (std::getline(input_file, line))
    {
      if (line.substr(0, m3vcf_version_line.size()) == m3vcf_version_line)
      {
        if (line == "##fileformat=M3VCFv2.0")
          m3vcf_version = 2;
        else
          m3vcf_version = 1;
        break;
      }

      if (line.size() < 2 || line[1] != '#')
      {
        std::cerr << "Error: invalid reference file" << std::endl;
        return false;
      }
    }

    std::size_t n_samples = 0;
    while (std::getline(input_file, line))
    {
      if (line.size() < 2 || line[1] != '#')
      {
        n_samples = std::count(line.begin(), line.end(), '\t') - 8;
        break;
      }
    }

    std::size_t ref_cnt = 0;

    auto tar_it = target_sites.begin();
    unique_haplotype_block block;
    std::vector<std::int8_t> tmp_geno;
    while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples))
    {
      if (block.variants().empty() || block.variants().front().pos > extended_reg.to())
        break;

      for (auto ref_it = block.variants().begin(); ref_it != block.variants().end(); ++ref_it)
      {
        while (tar_it != target_sites.end() && tar_it->pos < ref_it->pos)
        {
          if (ref_cnt)
          {
            tar_it->ref_cnt = ref_cnt;
            ref_cnt = 0;
          }
          ++tar_it;
        }

        for (auto it = tar_it; it != target_sites.end() && it->pos == ref_it->pos; ++it)
        {
          if (it->ref == ref_it->ref && it->alt == ref_it->alt)
          {
            tmp_geno.resize(block.unique_map().size());
            for (std::size_t i = 0; i < tmp_geno.size(); ++i)
              tmp_geno[i] = ref_it->gt[block.unique_map()[i]];
            typed_only_reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
            it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size();
            it->af = float((--typed_only_reference_data.end())->ac) / tmp_geno.size();
            it->in_ref = true;
            if (it != tar_it)
              std::swap(*it, *tar_it);
            if (ref_cnt)
            {
              tar_it->ref_cnt = ref_cnt;
              ref_cnt = 0;
            }
            ++tar_it;
            break;
          }
        }

        if (ref_it->pos >= extended_reg.from())
          ++ref_cnt;
      }

      if (full_reference_data)
      {
        // TODO: remove redundant variants
        block.trim(impute_reg.from(), impute_reg.to());
        if (!block.variants().empty())
          full_reference_data->append_block(block);
      }
    }

    if (ref_cnt && tar_it != target_sites.end())
    {
      tar_it->ref_cnt = ref_cnt;
      ref_cnt = 0;
    }
  }

  return !input.bad();
}

std::vector<target_variant> separate_target_only_variants(std::vector<target_variant>& target_sites)
{
  std::vector<target_variant> target_only_sites;
  std::size_t shift_idx = 0;
  for (std::size_t i = 0; i < target_sites.size(); ++i)
  {
    if (!target_sites[i].in_ref)
    {
      target_only_sites.emplace_back();
      std::swap(target_only_sites.back(), target_sites[i]);
    }
    else
    {
      if (i != shift_idx)
      {
        std::swap(target_sites[i], target_sites[shift_idx]);
      }
      ++shift_idx;
    }
  }

  target_sites.resize(shift_idx);
  return target_only_sites;
}

class haplotype_interpolator
{
private:
  savvy::writer out_file_;
  savvy::file::format file_format_;
  std::vector<std::string> fmt_fields_;
  std::size_t n_samples_ = 0;
  const std::int16_t bin_scalar_ = 100; //256;
  struct variant_update_ctx
  {
    savvy::compressed_vector<std::int8_t> sparse_gt;
    savvy::compressed_vector<float> sparse_dosages;
    std::vector<float> gp_vec;
  };
public:
  haplotype_interpolator(const std::string& file_path, savvy::file::format file_format, std::uint8_t out_compression, const std::vector<std::string>& sample_ids, const std::vector<std::string>& fmt_fields, const std::string& chromosome) :
    out_file_(file_path, file_format, gen_headers(fmt_fields, chromosome), sample_ids, out_compression),
    file_format_(file_format),
    fmt_fields_(fmt_fields),
    n_samples_(sample_ids.size())
  {

  }

  static std::vector<std::pair<std::string, std::string>> gen_headers(const std::vector<std::string>& fmt_fields, const std::string& chromosome)
  {
    std::time_t t = std::time(nullptr);
    char datestr[11];
    assert(std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));

    std::vector<std::pair<std::string, std::string>> headers = {
      {"fileformat","VCFv4.2"},
      {"filedate", datestr},
      {"source","Minimac v" + std::string(VERSION)},
      {"phasing","full"},
      {"contig","<ID=" + std::string(chromosome) + ">"},
      {"INFO","<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">"},
      {"INFO","<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">"},
      {"INFO","<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">"},
      {"INFO","<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">"},
      {"INFO","<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">"},
      {"INFO","<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">"},
      {"INFO","<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Marker was genotyped but NOT imputed\">"}};

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

  bool sites_match(const target_variant& t, const reference_site_info& r)
  {
    return t.pos == r.pos && t.alt == r.alt && t.ref == r.ref;
  }

  bool merge_temp_files(const std::vector<std::string>& temp_files_paths)
  {
    if (temp_files_paths.empty())
      return false;

    std::list<savvy::reader> temp_files;
    for (auto it = temp_files_paths.begin(); it != temp_files_paths.end(); ++it)
    {
      temp_files.emplace_back(*it);
      std::remove(it->c_str());
    }


    savvy::variant out_var;
    savvy::compressed_vector<float> pasted_hds;
    savvy::compressed_vector<float> partial_hds;

    int good_count = temp_files.size();
    while (good_count == temp_files.size())
    {
      pasted_hds.clear();
      good_count = 0;
      for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
      {
        good_count += (int)(*it >> out_var).good();
        out_var.get_format("HDS", partial_hds);
        std::size_t old_size = pasted_hds.size();
        pasted_hds.resize(old_size + partial_hds.size());
        for (auto jt = partial_hds.begin(); jt != partial_hds.end(); ++jt)
          pasted_hds[old_size + jt.offset()] = *jt;
      }

      if (good_count)
      {
        if (good_count < temp_files.size())
          return std::cerr << "Error: record mismatch in temp files" << std::endl, false;

        std::size_t n = pasted_hds.size();
        float s_x = std::accumulate(pasted_hds.begin(), pasted_hds.end(), 0.f);
        float s_xx = std::inner_product(pasted_hds.begin(), pasted_hds.end(), pasted_hds.begin(), 0.f);
        float af = s_x / n;

        float r2 = savvy::typed_value::missing_value<float>();
        if (af > 0.f && af < 1.f)
          r2 = ((s_xx - s_x * s_x / n) / n) / (af * (1.f - af));

        out_var.set_info("AF", af);
        out_var.set_info("MAF", af > 0.5f ? 1.f - af : af);
        out_var.set_info("R2", r2);

        out_var.set_format("HDS", pasted_hds);
        out_file_ << out_var;
      }
    }

    int bad_count = 0;
    for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
      bad_count += (int)it->bad();

    if (bad_count || !out_file_.good())
      return std::cerr << "Error: I/O failed during merging" << std::endl, false;
    return true;
  }

  bool write_dosages(const full_dosages_results& hmm_results, const std::vector<target_variant>& tar_variants, std::pair<std::size_t, std::size_t> oberved_range, const reduced_haplotypes& full_reference_data, const std::string& output_path, omp::internal::thread_pool2& tpool)
  {
    assert(hmm_results.dimensions()[0] == full_reference_data.variant_size());

    std::vector<variant_update_ctx> update_contexts(tpool.thread_count());
    std::vector<savvy::variant> out_vars(tpool.thread_count());

    assert(!tar_variants.empty());
    std::vector<std::vector<float>> dosage_vecs(tpool.thread_count(), std::vector<float>(tar_variants[0].gt.size()));
    std::vector<std::vector<float>> loo_dosage_vecs(tpool.thread_count(), std::vector<float>(tar_variants[0].gt.size()));
    const std::vector<std::int8_t> empty_gt_vec;
    const std::vector<float> empty_loo_vec;

    std::vector<std::size_t> tar_indices(tpool.thread_count());
    std::size_t tar_idx = 0;
    std::size_t i = 0;
#if 0
    auto ref_it = full_reference_data.begin();
    auto ref_it_end = full_reference_data.end();
    while (i < full_reference_data.variant_size())
    {
      assert(i == ref_it.global_idx());
      std::fill(tar_indices.begin(), tar_indices.end(), tar_idx);
      tpool([this, i, tar_idx, ref_it, ref_it_end, &tar_indices, &tar_variants, &hmm_results, &dosage_vecs, &loo_dosage_vecs, &out_vars, &update_contexts, &empty_loo_vec, &empty_gt_vec](std::size_t thread_idx)
        {
          std::size_t ref_idx = i + thread_idx;
          std::size_t local_tar_idx = tar_idx;
          auto local_ref_it = ref_it;
          auto& dosages = dosage_vecs[thread_idx];
          auto& loo_dosages = loo_dosage_vecs[thread_idx];

          for (std::size_t j = 0; j < thread_idx && local_ref_it != ref_it_end; ++j)
            ++local_ref_it;

          if (local_ref_it == ref_it_end)
            return;

          while (local_tar_idx < tar_variants.size() && tar_variants[local_tar_idx].pos < local_ref_it->pos)
            ++local_tar_idx;

          tar_indices[thread_idx] = local_tar_idx;

          bool ref_matches_tar = local_tar_idx < tar_variants.size() && sites_match(tar_variants[local_tar_idx], *local_ref_it);

//          for (std::size_t j = 0; j < dosages.size(); ++j)
//            dosages[j] = float(std::int16_t(hmm_results.dosage(ref_idx, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
//
//          if (ref_matches_tar)
//          {
//            for (std::size_t j = 0; j < loo_dosages.size(); ++j)
//              loo_dosages[j] = float(std::int16_t(hmm_results.loo_dosage(local_tar_idx, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
//          }


          out_vars[thread_idx] = savvy::site_info(local_ref_it->chrom, local_ref_it->pos, local_ref_it->ref, {local_ref_it->alt}, ""/*ref_var.id()*/);
          prepare_output_variant(out_vars[thread_idx], update_contexts[thread_idx], hmm_results.dosages_[ref_idx], ref_matches_tar ? hmm_results.loo_dosages_[local_tar_idx] : empty_loo_vec, ref_matches_tar ? tar_variants[local_tar_idx].gt : empty_gt_vec);
        });

      for (std::size_t j = 0; j < tpool.thread_count() && ref_it != ref_it_end; ++j,++ref_it,++i)
      {
        out_file_ << out_vars[j];
        if (tar_indices[j] > tar_idx)
          tar_idx = tar_indices[j];
      }
    }
#else
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



      out_vars[0] = savvy::site_info(ref_it->chrom, ref_it->pos, ref_it->ref, {ref_it->alt}, ""/*ref_var.id()*/);
      prepare_output_variant(out_vars[0], update_contexts[0], hmm_results.dosages_[i],
        ref_matches_tar ? hmm_results.loo_dosages_[tar_idx] : empty_loo_vec,
        ref_matches_tar ? std::vector<std::int8_t>(tar_variants[tar_idx].gt.begin() + oberved_range.first, tar_variants[tar_idx].gt.begin() + oberved_range.second) : empty_gt_vec);
      out_file_ << out_vars[0];
    }
#endif
    return out_file_.good();
  }

  bool piecewise_constant( const best_templates_results& hmm_results, const std::vector<target_variant>& typed_variants, const reduced_haplotypes& reduced_reference, const std::string& reference_path, const savvy::genomic_region& reg, const std::string& output_path)
  {
    auto dims = hmm_results.dimensions();
    savvy::variant ref_var, out_var;
    variant_update_ctx update_ctx;
    savvy::reader ref_file(reference_path);
    ref_file.reset_bounds(reg);

    //savvy::writer out_file(output_path, savvy::file::format::bcf, {}, {});

    std::vector<float> dosages(dims[1]);
    std::vector<float> loo_dosages(dims[1]);
    const std::vector<float> empty_loo_vec;
    const std::vector<std::int8_t> empty_gt_vec;

    savvy::compressed_vector<std::int8_t> genos;
    std::vector<std::size_t> allele_counts;
    std::vector<std::int8_t> dense_genos;

    auto sites_match = [](const target_variant& t, const savvy::variant& r)
    {
      return t.pos == r.position() && t.alt == r.alts()[0] && t.ref == r.ref();
    };

    ref_file >> ref_var;

    std::size_t global_idx = 0, block_idx = 0, local_idx = 0;
    while (global_idx < typed_variants.size() && ref_file)
    {
      bool ref_matches_typed = sites_match(typed_variants[global_idx], ref_var);

      const auto& ref_block = reduced_reference.blocks()[block_idx];

      std::size_t next_global_idx = global_idx + 1, next_block_idx = block_idx, next_local_idx = local_idx + 1;
      if (next_local_idx >= ref_block.variant_size())
      {
        ++next_block_idx;
        next_local_idx = 0;
      }
      const auto& next_ref_block = reduced_reference.blocks()[next_block_idx];

      std::size_t mid_pos = (typed_variants[next_global_idx].pos - typed_variants[global_idx].pos) / 2 + typed_variants[global_idx].pos;
      if (global_idx + 1 < typed_variants.size() && typed_variants[global_idx + 1].pos == typed_variants[global_idx].pos)
      {
        // TODO: handle multiallelic interpolation
      }

      if (reg.from() <= mid_pos) // TODO: how to deal with multiallelics
      {
        // TODO: write typed only sites


        while (ref_file && (ref_var.position() < mid_pos || ref_matches_typed))
        {
          if (ref_var.position() == 10000210)
          {
            auto a = 0;
          }
          ref_var.get_format("GT", genos);
          if (dense_genos.empty())
            dense_genos.resize(genos.size());
          assert(dense_genos.size() == genos.size());

          allele_counts.resize(ref_block.unique_haplotype_size());
          std::fill(allele_counts.begin(), allele_counts.end(), 0);
          for (auto it = genos.begin(); it != genos.end(); ++it)
          {
            ++allele_counts[ref_block.unique_map()[it.offset()]]; // reference panel is expected to be separated into biallelic variants
            dense_genos[it.offset()] = *it;
          }
          assert(std::accumulate(allele_counts.begin(), allele_counts.end(), 0) == std::accumulate(genos.begin(), genos.end(), 0));
          assert(std::accumulate(allele_counts.begin(), allele_counts.end(), 0) == genos.non_zero_size());
          if (!ref_matches_typed)
          {
            for (std::size_t j = 0; j < dims[1]; ++j)
            {
              float p = 0.f;
              std::size_t ac = 0;
              std::size_t an = 0;
              float prob_sum = 0.f;
              const auto& best_templates = hmm_results.best_templates(global_idx, j);
              const auto& best_probs = hmm_results.best_probs(global_idx, j);
              const auto& best_uniq_flags = hmm_results.best_uniq_flags(global_idx, j);
              assert(best_templates.size() == best_probs.size());

              for (std::size_t k = 0; k < best_templates.size(); ++k)
              {
                if (best_uniq_flags[k])
                {
                  // best_templates[k] is a unique haplotype
                  std::size_t uniq_idx = best_templates[k];
                  ac += allele_counts[uniq_idx];
                  prob_sum += best_probs[k];
                  an += ref_block.cardinalities()[uniq_idx];
                  p += best_probs[k] * (float(allele_counts[uniq_idx]) / float(ref_block.cardinalities()[uniq_idx]));
                }
                else
                {
                  std::int8_t g = dense_genos[best_templates[k]];
                  p += best_probs[k] * g;
                  prob_sum += best_probs[k];
                  ac += g;
                  ++an;
                }
              }

              p += (1.f - prob_sum) * (float(genos.non_zero_size() - ac) / float(genos.size() - an));
              p = std::max(0.f, std::min(1.f, p));
              p = float(std::int16_t(p * bin_scalar_ + 0.5f)) / bin_scalar_;
              // char s[256];
              // std::sprintf(s, "%0.3f", p);
              dosages[j] = p; // TODO: this may not correct
              assert(dosages[j] <= 1.f);
              assert(dosages[j] >= 0.f);
            }
          }
          else
          {
            for (std::size_t j = 0; j < dims[1]; ++j)
            {
              dosages[j] = float(std::int16_t(hmm_results.dosage(global_idx, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
              loo_dosages[j] = float(std::int16_t(hmm_results.loo_dosage(global_idx, j) * bin_scalar_ + 0.5f)) / bin_scalar_;
            }
          }

          for (auto it = genos.begin(); it != genos.end(); ++it)
            dense_genos[it.offset()] = 0; // make dense_genos all zero.

          out_var = savvy::site_info(ref_var.chromosome(), ref_var.position(), ref_var.ref(), ref_var.alts(), ref_var.id());
          prepare_output_variant(out_var, update_ctx, dosages, ref_matches_typed ? loo_dosages : empty_loo_vec, ref_matches_typed ? typed_variants[global_idx].gt : empty_gt_vec);
          out_file_ << out_var;

          ref_file >> ref_var;
          ref_matches_typed = sites_match(typed_variants[global_idx], ref_var);
        }
      }

      ++global_idx;
      ++local_idx;
      if (local_idx >= ref_block.variant_size())
      {
        ++block_idx;
        local_idx = 0;
      }
    }
    return out_file_.good();
  }
private:
  void prepare_output_variant(savvy::variant& out_var, variant_update_ctx& ctx, const std::vector<float>& dosages, const std::vector<float>& loo_dosages, const std::vector<std::int8_t>& observed)
  {
    ctx.sparse_dosages.assign(dosages.begin(), dosages.end());

    std::size_t n = ctx.sparse_dosages.size();
    float s_x = std::accumulate(ctx.sparse_dosages.begin(), ctx.sparse_dosages.end(), 0.f);
    float s_xx = std::inner_product(ctx.sparse_dosages.begin(), ctx.sparse_dosages.end(), ctx.sparse_dosages.begin(), 0.f);
    float af = s_x / n;

    float r2 = savvy::typed_value::missing_value<float>();
    if (af > 0.f && af < 1.f)
      r2 = ((s_xx - s_x * s_x / n) / n) / (af * (1.f - af));

    out_var.set_info("AF", af);
    out_var.set_info("MAF", af > 0.5f ? 1.f - af : af);
    out_var.set_info("R2", r2);

    bool hds_set = true;
    for (const auto& fmt : fmt_fields_)
    {
      if (fmt == "HDS")
      {
        if (!hds_set)
          ctx.sparse_dosages.assign(dosages.begin(), dosages.end());
        out_var.set_format("HDS", ctx.sparse_dosages);
        hds_set = true;
      }
      else if (fmt == "GT")
      {
        ctx.sparse_gt.assign(dosages.begin(), dosages.end(), [](float v){ return std::int8_t(v < 0.5f ? 0 : 1); });
        out_var.set_format("GT", ctx.sparse_gt);
      }
      else if (fmt == "DS")
      {
        if (!hds_set)
          ctx.sparse_dosages.assign(dosages.begin(), dosages.end());
        savvy::stride_reduce(ctx.sparse_dosages, ctx.sparse_dosages.size() / n_samples_);
        out_var.set_format("DS", ctx.sparse_dosages);
        hds_set = false;
      }
      else if (fmt == "GP")
      {
      }
      else if (fmt == "SD")
      {
      }
    }

    if (observed.size())
    {
      assert(observed.size() == loo_dosages.size());
      ctx.sparse_dosages.assign(loo_dosages.begin(), loo_dosages.end());
      ctx.sparse_gt.assign(observed.begin(), observed.end());
      s_x = std::accumulate(ctx.sparse_dosages.begin(), ctx.sparse_dosages.end(), 0.f);
      s_xx = std::inner_product(ctx.sparse_dosages.begin(), ctx.sparse_dosages.end(), ctx.sparse_dosages.begin(), 0.f);
      float s_y = std::accumulate(ctx.sparse_gt.begin(), ctx.sparse_gt.end(), 0.f);
      // since observed can only be 0 or 1, s_yy is the same as s_y
      float s_yy = s_y; //std::inner_product(sparse_gt.begin(), sparse_gt.end(), sparse_gt.begin(), 0.f); // TODO: allow for missing oberserved genotypes.
      float s_xy = 0.f;
      for (auto it = ctx.sparse_gt.begin(); it != ctx.sparse_gt.end(); ++it)
        s_xy += *it * loo_dosages[it.offset()];

      //                         n * Sum xy - Sum x * Sum y
      //  r = -------------------------------------------------------------------
      //      Sqrt(n * Sum xx - Sum x * Sum x) * Sqrt(n * Sum yy - Sum y * Sum y)
      float emp_r = (n * s_xy - s_x * s_y) / (std::sqrt(n * s_xx - s_x * s_x) * std::sqrt(n * s_yy - s_y * s_y));
      out_var.set_info("ER2", std::isnan(emp_r) ? savvy::typed_value::missing_value<float>() : emp_r * emp_r);
    }

    out_var.set_info(observed.size() ? "TYPED" : "IMPUTED", std::vector<std::int8_t>());
  }
};

bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path)
{
  std::vector<std::pair<std::string, std::string>> headers;
  std::vector<std::string> ids;

  shrinkwrap::gz::istream input_file(input_path);
  std::string line;

  std::uint8_t m3vcf_version = 0;
  const std::string m3vcf_version_line = "##fileformat=M3VCF";
  const std::string vcf_version_line = "##fileformat=VCF";
  while (std::getline(input_file, line))
  {
    std::size_t equal_pos = line.find('=');
    if (equal_pos != std::string::npos)
    {
      std::string key = line.substr(0, equal_pos);
      std::string val = line.substr(equal_pos + 1);
      key.erase(0, key.find_first_not_of('#'));

      if (key == "fileformat")
      {
        if (val.substr(0, 5) == "M3VCF")
        {
          if (val == "M3VCFv2.0")
            m3vcf_version = 2;
          else
            m3vcf_version = 1;
        }
      }
      else
      {
        headers.emplace_back(std::move(key), std::move(val));
      }
    }
    else
    {

      break;
    }
  }

  if (line.size() < 1 || line[0] != '#')
  {
    std::cerr << "Error: invalid reference file" << std::endl;
    return false;
  }

  headers.insert(headers.begin(), {"subfileformat","M3VCFv3.0"});
  headers.insert(headers.begin(), {"fileformat","VCFv4.2"});
  headers.emplace_back("INFO","<ID=REPS,Number=1,Type=Integer,Description=\"Number of distinct haplotypes in block\">");
  headers.emplace_back("INFO","<ID=VARIANTS,Number=1,Type=Integer,Description=\"Number of variants in block\">");
  headers.emplace_back("INFO","<ID=END,Number=1,Type=Integer,Description=\"End position of record\">");
  headers.emplace_back("INFO","<ID=UHA,Number=.,Type=Integer,Description=\"Unique haplotype allele\">");
  headers.emplace_back("FORMAT","<ID=UHM,Number=.,Type=Integer,Description=\"Unique haplotype mapping\">");

  //headers.emplace_back("ALT","<ID=DUP,Description=\"Duplication\">");


  std::size_t tab_cnt = 0;
  std::size_t last_pos = 0;
  std::size_t tab_pos;
  while ((tab_pos = line.find('\t', tab_pos)) != std::string::npos)
  {
    if (tab_cnt >= 9)
    {
      ids.emplace_back(line.substr(last_pos, tab_pos - last_pos));
    }
    last_pos = ++tab_pos;
    ++tab_cnt;
  }

  ids.emplace_back(line.substr(last_pos, tab_pos - last_pos));
  std::size_t n_samples = ids.size();

  savvy::writer output_file(output_path, savvy::file::format::bcf, headers, ids, 6);


  std::size_t block_cnt = 0;
  unique_haplotype_block block;
  std::vector<std::int8_t> tmp_geno;
  while (block.deserialize(input_file, m3vcf_version, m3vcf_version == 1 ? n_samples : 2 * n_samples))
  {
    if (block.variants().empty())
      break;

    block.serialize(output_file);
    ++block_cnt;
  }

  return !input_file.bad() && output_file.good();
}

int main(int argc, char** argv)
{
  //return convert_old_m3vcf("../1000g-test/topmed.chr20.gtonly.filtered.rehdr.pass_only.phased.ligated.shuffled.m3vcf.gz", "../1000g-test/topmed.chr20.gtonly.filtered.rehdr.pass_only.phased.ligated.shuffled.m3bcf") ? 0 : 1;
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);

    return EXIT_SUCCESS;
  }

  savvy::region extended_region =
    {
      args.region().chromosome(),
      std::uint64_t(std::max(std::int64_t(1), std::int64_t(args.region().from()) - args.overlap())),
      args.region().to() + args.overlap()
    };

  std::time_t start_time = std::time(nullptr);
  std::vector<std::string> sample_ids;
  std::vector<target_variant> target_sites;
  load_target_haplotypes(args.tar_path(), extended_region, target_sites, sample_ids);
  std::cerr << ("Loading target haplotypes took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  start_time = std::time(nullptr);

  reduced_haplotypes typed_only_reference_data(16, 512);
  reduced_haplotypes full_reference_data;
  load_reference_haplotypes(args.ref_path(), extended_region, args.region(), target_sites, typed_only_reference_data, args.deferred_interpolation() ? nullptr : &full_reference_data);
  std::cerr << ("Loading reference haplotypes took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  start_time = std::time(nullptr);
  std::vector<target_variant> target_only_sites = separate_target_only_variants(target_sites);
  std::cerr << ("Separating typed only variants took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  if (target_sites.empty())
    return std::cerr << "Error: no target variants\n", EXIT_FAILURE;

  target_sites.back().recom = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
  if (args.map_path().size())
  {
//    start_time = std::time(nullptr);
//    std::vector<reference_site_info> ref_sites;
//    ref_sites.reserve(full_reference_data.variant_size());
//    for (auto it = full_reference_data.begin(); it != full_reference_data.end(); ++it)
//      ref_sites.emplace_back(*it);

    if (!recombination::parse_map_file(args.map_path(), target_sites.begin(), target_sites.end()))
    //if (!recombination::parse_map_file(args.map_path(), ref_sites.begin(), ref_sites.end()))
      return std::cerr << "Error: parsing map file failed\n", EXIT_FAILURE;

#if 0
    auto tar_it = target_sites.begin();
    for (auto it = target_sites.begin(); it != target_sites.end(); ++it)
    {
      auto prev_it = it++;
      if (it != target_sites.end())
      {
        assert(prev_it->ref_cnt > 0);
        float recom = std::max(recombination::recom_min, prev_it->recom / prev_it->ref_cnt);
        float temp = (1.f - recom);
        for (std::size_t i = 0; i < prev_it->ref_cnt; ++i)
          temp *= (1.f - recom);
        prev_it->recom = 1.f - temp;
        assert(prev_it->recom <= 1.f && prev_it->recom >= 0.f);
      }
    }
#endif
    std::cerr << ("Loading switch probabilities took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }

  std::vector<std::vector<std::vector<std::size_t>>> reverse_maps;
  reverse_maps.reserve(typed_only_reference_data.blocks().size());
  for (auto it = typed_only_reference_data.blocks().begin(); it != typed_only_reference_data.blocks().end(); ++it)
  {
    reverse_maps.emplace_back();
    auto& map = reverse_maps.back();
    for (std::size_t i = 0; i < it->cardinalities().size(); ++i)
    {
      map.emplace_back();
      map.back().reserve(it->cardinalities()[i]);
    }

    for (std::size_t i = 0; i < it->unique_map().size(); ++i)
    { assert(it->unique_map()[i] < map.size());
      map[it->unique_map()[i]].push_back(i);
    }
  }

  std::cerr << "Running HMM with " << args.threads() << " threads ..." << std::endl;
  start_time = std::time(nullptr);
  omp::internal::thread_pool2 tpool(args.threads());
  std::vector<hidden_markov_model> hmms(args.threads());

  if (args.deferred_interpolation())
  {
    best_templates_results hmm_results;
    hmm_results.resize(target_sites.size(), target_sites[0].gt.size());

    omp::parallel_for_exp(omp::static_schedule(), omp::sequence_iterator(0), omp::sequence_iterator(target_sites[0].gt.size()), [&](int& i, const omp::iteration_context& ctx)
      {
        hmms[ctx.thread_index].traverse_forward(typed_only_reference_data.blocks(), target_sites, i);
        hmms[ctx.thread_index].traverse_backward(typed_only_reference_data.blocks(), target_sites, i, i, reverse_maps, hmm_results);
      },
      tpool);
    std::cerr << ("Running HMM took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
    start_time = std::time(nullptr);
    haplotype_interpolator interpolator(args.out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom);
    interpolator.piecewise_constant(hmm_results, target_sites, typed_only_reference_data, args.ref_path(), args.region(), args.out_path());
    std::cerr << ("Interpolating and writing output took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }
  else
  {
    std::size_t haplotype_buffer_size = 400; // TODO
    std::vector<std::string> temp_file_paths;
    temp_file_paths.reserve(target_sites[0].gt.size() / haplotype_buffer_size + 1);
    full_dosages_results hmm_results;
    hmm_results.resize(full_reference_data.variant_size(), target_sites.size(), std::min(haplotype_buffer_size, target_sites[0].gt.size()));

    double impute_time = 0.;
    double temp_write_time = 0.;
    bool use_temp_files = target_sites[0].gt.size() > haplotype_buffer_size;
    for (std::size_t i = 0; i < target_sites[0].gt.size(); i += haplotype_buffer_size)
    {
      std::size_t group_size = target_sites[0].gt.size() - i;
      if (group_size < haplotype_buffer_size)
      {
        for (std::size_t j = 0; j < hmm_results.dosages_.size(); ++j)
          hmm_results.dosages_[j].resize(group_size);
        for (std::size_t j = 0; j < hmm_results.loo_dosages_.size(); ++j)
          hmm_results.loo_dosages_[j].resize(group_size);
      }

      start_time = std::time(nullptr);
      omp::parallel_for_exp(omp::static_schedule(), omp::sequence_iterator(i), omp::sequence_iterator(std::min(i + haplotype_buffer_size, target_sites[0].gt.size())), [&](int& i, const omp::iteration_context& ctx)
        {
          hmms[ctx.thread_index].traverse_forward(typed_only_reference_data.blocks(), target_sites, i);
          auto reverse_ref_itr = --full_reference_data.end();
          std::size_t block_idx(-1);
          hmms[ctx.thread_index].traverse_backward(typed_only_reference_data.blocks(), target_sites, i, i % haplotype_buffer_size, reverse_maps, hmm_results, reverse_ref_itr, --full_reference_data.begin(),
            block_idx);
        }, tpool);
      impute_time += std::difftime(std::time(nullptr), start_time);

      start_time = std::time(nullptr);
      temp_file_paths.push_back(use_temp_files ? "/tmp/m4_g" + std::to_string(i / haplotype_buffer_size) + ".sav" : args.out_path()); //TODO: use real temp files
      haplotype_interpolator interpolator(temp_file_paths.back(),
        use_temp_files ? savvy::file::format::sav : args.out_format(),
        use_temp_files ? std::min<std::uint8_t>(3, args.out_compression()) : args.out_compression(),
        sample_ids,
        use_temp_files ? std::vector<std::string>{"HDS"} : args.fmt_fields(),
        target_sites.front().chrom);
      interpolator.write_dosages(hmm_results, target_sites, {i, std::min(i + haplotype_buffer_size, target_sites[0].gt.size())}, full_reference_data, args.out_path(), tpool);
      temp_write_time += std::difftime(std::time(nullptr), start_time);
    }

    std::cerr << ("Running HMM took " + std::to_string(impute_time) + " seconds") << std::endl;
    if (use_temp_files)
    {
      std::cerr << ("Writing temp files took " + std::to_string(temp_write_time) + " seconds") << std::endl;

      // TODO:
      std::cerr << "Merging temp files ... " << std::endl;
      start_time = std::time(nullptr);
      haplotype_interpolator interpolator(args.out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom);
      interpolator.merge_temp_files(temp_file_paths);
      std::cerr << ("Merging temp files took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
    }
    else
    {
      std::cerr << ("Writing output took " + std::to_string(temp_write_time) + " seconds") << std::endl;
    }
  }



  return EXIT_SUCCESS;
}
