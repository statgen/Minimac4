
#include "getopt_wrapper.hpp"
#include "input_prep.hpp"
#include "hidden_markov_model.hpp"
#include "recombination.hpp"
#include "dosage_writer.hpp"

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
  std::size_t temp_buffer_ = 200;
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
  std::size_t temp_buffer() const { return temp_buffer_ ; }

  prog_args() :
    getopt_wrapper(
      "Usage: minimac4 [opts ...] <reference.m3vcf.gz> <target.vcf.gz>",
      {
        {"temp-buffer", required_argument, 0, 'b', "Number of samples to impute before writing to temporary files (default: 200)"},
        {"help", no_argument, 0, 'h', "Print usage"},
        {"format", required_argument, 0, 'f', "Comma-separated list of format fields to generate (GT, HDS, DS, GP, or SD; default: HDS)"},
        {"map", required_argument, 0, 'm', "Genetic map file"},
        {"output", required_argument, 0, 'o', "Output path (default: /dev/stdout)"},
        {"output-format", required_argument, 0, 'O', "Output file format (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: bcf)"},
        {"region", required_argument, 0, 'r', "Genomic region to impute"},
        {"threads", required_argument, 0, 't', "Number of threads (default: 1)"},
        {"overlap", required_argument, 0, 'w', "Size (in basepairs) of overlap before and after region to use as input to HMM (default: 1000000)"}
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
      case 'b':
        temp_buffer_ = std::size_t(std::atoll(optarg ? optarg : ""));
        break;
      case 'h':
        help_ = true;
        return true;
      case 'f':
        fmt_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
        break;
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

  std::vector<std::string> split_string_to_vector(const char* in, char delim)
  {
    std::vector<std::string> ret;
    const char* d = nullptr;
    std::string token;
    const char* s = in;
    const char*const e = in + strlen(in);
    while ((d = std::find(s, e,  delim)) != e)
    {
      ret.emplace_back(std::string(s, d));
      s = d ? d + 1 : d;
    }
    ret.emplace_back(std::string(s,d));
    return ret;
  }
};

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
  load_reference_haplotypes(args.ref_path(), extended_region, args.region(), target_sites, typed_only_reference_data, full_reference_data);
  std::cerr << ("Loading reference haplotypes took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  start_time = std::time(nullptr);
  std::vector<target_variant> target_only_sites = separate_target_only_variants(target_sites);
  std::cerr << ("Separating typed only variants took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  if (target_sites.empty())
    return std::cerr << "Error: no target variants\n", EXIT_FAILURE;

  target_sites.back().recom = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
  if (args.map_path().size())
  {
    start_time = std::time(nullptr);
    if (!recombination::parse_map_file(args.map_path(), target_sites.begin(), target_sites.end()))
      return std::cerr << "Error: parsing map file failed\n", EXIT_FAILURE;
    std::cerr << ("Loading switch probabilities took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }

  auto reverse_maps = generate_reverse_maps(typed_only_reference_data);

  std::cerr << "Running HMM with " << args.threads() << " threads ..." << std::endl;
  start_time = std::time(nullptr);
  omp::internal::thread_pool2 tpool(args.threads());
  std::vector<hidden_markov_model> hmms(args.threads());

  std::size_t ploidy = target_sites[0].gt.size() / sample_ids.size();
  std::size_t haplotype_buffer_size = args.temp_buffer() * ploidy;
  assert(ploidy && target_sites[0].gt.size() % sample_ids.size() == 0);
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
    dosage_writer output(temp_file_paths.back(),
      use_temp_files ? savvy::file::format::sav : args.out_format(),
      use_temp_files ? std::min<std::uint8_t>(3, args.out_compression()) : args.out_compression(),
      {sample_ids.begin() + (i / ploidy), sample_ids.begin() + (i + group_size) / ploidy},
      use_temp_files ? std::vector<std::string>{"HDS"} : args.fmt_fields(),
      target_sites.front().chrom);
    output.write_dosages(hmm_results, target_sites, {i, std::min(i + haplotype_buffer_size, target_sites[0].gt.size())}, full_reference_data);
    temp_write_time += std::difftime(std::time(nullptr), start_time);
  }

  std::cerr << ("Running HMM took " + std::to_string(impute_time) + " seconds") << std::endl;
  if (use_temp_files)
  {
    std::cerr << ("Writing temp files took " + std::to_string(temp_write_time) + " seconds") << std::endl;

    // TODO:
    std::cerr << "Merging temp files ... " << std::endl;
    start_time = std::time(nullptr);
    dosage_writer output(args.out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom);
    output.merge_temp_files(temp_file_paths);
    std::cerr << ("Merging temp files took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }
  else
  {
    std::cerr << ("Writing output took " + std::to_string(temp_write_time) + " seconds") << std::endl;
  }




  return EXIT_SUCCESS;
}
