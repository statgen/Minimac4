
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

#include <cstdio>
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
  std::string prefix_; // deprecated
  std::string emp_out_path_;
  std::string sites_out_path_;
  savvy::file::format out_format_ = savvy::file::format::bcf;
  std::uint8_t out_compression_ = 6;
  std::vector<std::string> fmt_fields_ = {"GT","HDS","DS"};
  savvy::genomic_region reg_ = {""};
  std::size_t temp_buffer_ = 200;
  std::int64_t chunk_size_ = 20000000;
  std::int64_t overlap_ = 3000000;
  std::int16_t threads_ = 1;
  bool all_typed_sites_ = false;
  bool pass_only_ = false;
  bool meta_ = false; // deprecated
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
  const std::string& emp_out_path() const { return emp_out_path_; }
  const std::string& sites_out_path() const { return sites_out_path_; }
  savvy::file::format out_format() const { return out_format_; }
  std::uint8_t out_compression() const { return out_compression_; }
  const std::vector<std::string>& fmt_fields() const { return fmt_fields_; }
  const savvy::genomic_region& region() const { return reg_; }
  std::int64_t chunk_size() const { return chunk_size_; }
  std::int64_t overlap() const { return overlap_; }
  std::int16_t threads() const { return threads_; }
  std::size_t temp_buffer() const { return temp_buffer_ ; }
  bool all_typed_sites() const { return all_typed_sites_; }
  bool pass_only() const { return pass_only_; }

  prog_args() :
    getopt_wrapper(
      "Usage: minimac4 [opts ...] <reference.m3vcf.gz> <target.vcf.gz>",
      {
        {"all-typed-sites", no_argument, 0, 'a', "Include in the output sites that exist only in target VCF"},
        {"temp-buffer", required_argument, 0, 'b', "Number of samples to impute before writing to temporary files (default: 200)"},
        {"chunk", required_argument, 0, 'c', "Maximum chunk length in base pairs to impute at once"}, // TODO
        {"empirical-output", required_argument, 0, 'e', "Output path for empirical dosages"},
        {"help", no_argument, 0, 'h', "Print usage"},
        {"format", required_argument, 0, 'f', "Comma-separated list of format fields to generate (GT, HDS, DS, GP, or SD; default: HDS)"},
        {"map", required_argument, 0, 'm', "Genetic map file"},
        {"output", required_argument, 0, 'o', "Output path (default: /dev/stdout)"},
        {"output-format", required_argument, 0, 'O', "Output file format (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: bcf)"},
        {"pass-only", no_argument, 0, 'p', "Only imports variants with FILTER column set to PASS"},
        {"region", required_argument, 0, 'r', "Genomic region to impute"},
        {"sites", required_argument, 0, 's', "Output path for sites-only file"},
        {"threads", required_argument, 0, 't', "Number of threads (default: 1)"},
        {"overlap", required_argument, 0, 'w', "Size (in base pairs) of overlap before and after impute region to use as input to HMM (default: 3000000)"},
        // vvvv deprecated vvvv //
        {"allTypedSites", no_argument, 0, '\x01', ""},
        {"rsid", no_argument, 0, '\x01', ""},
        {"passOnly", no_argument, 0, '\x01', ""},
        {"meta", no_argument, 0, '\x01', ""},
        {"haps", required_argument, 0, '\x02', ""},
        {"refHaps", required_argument, 0, '\x02', ""},
        {"prefix", required_argument, 0, '\x02', ""},
        {"mapFile", required_argument, 0, '\x02', ""},
        {"chr", required_argument, 0, '\x02', ""},
        {"start", required_argument, 0, '\x02', ""},
        {"end", required_argument, 0, '\x02', ""},
        {"window", required_argument, 0, '\x02', ""},
        {"ChunkOverlapMb", required_argument, 0, '\x02', ""},
        {"cpus", required_argument, 0, '\x02', ""}
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
      case 'a':
        all_typed_sites_ = true;
        break;
      case 'b':
        temp_buffer_ = std::size_t(std::atoll(optarg ? optarg : ""));
        break;
      case 'c':
        chunk_size_ = std::atoll(optarg ? optarg : "");
        break;
      case 'e':
        emp_out_path_ = optarg ? optarg : "";
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
      case 'p':
        pass_only_ = true;
        break;
      case 'r':
        reg_ = string_to_region(optarg ? optarg : "");
        break;
      case 's':
        sites_out_path_ = optarg ? optarg : "";
        break;
      case 't':
        threads_ = atoi(optarg ? optarg : "");
        break;
      case 'w':
        overlap_ = std::atoll(optarg ? optarg : "");
        break;
      case '\x01':
        if (std::string(long_options_[long_index].name) == "allTypedSites")
        {
          std::cerr << "Warning: --allTypedSites is deprecated in favor of --all-typed-sites\n";
          all_typed_sites_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "rsid")
        {
          std::cerr << "Warning: --rsid is deprecated (on by default)\n";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "passOnly")
        {
          std::cerr << "Warning: --passOnly is deprecated in favor of --pass-only\n";
          pass_only_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "meta")
        {
          std::cerr << "Warning: --meta is deprecated in favor of --empirical-output\n";
          meta_ = true;
          break;
        }
        // else pass through to default
      case '\x02':
        if (std::string(long_options_[long_index].name) == "haps")
        {
          std::cerr << "Warning: --haps is deprecated\n";
          tar_path_ = optarg ? optarg : "";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "refHaps")
        {
          std::cerr << "Warning: --refHaps is deprecated\n";
          ref_path_ = optarg ? optarg : "";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "chr")
        {
          std::cerr << "Warning: --chr is deprecated in favor of --region\n";
          reg_ = savvy::genomic_region(optarg ? optarg : "", reg_.from(), reg_.to());
          break;
        }
        else if (std::string(long_options_[long_index].name) == "start")
        {
          std::cerr << "Warning: --start is deprecated in favor of --region\n";
          reg_ = savvy::genomic_region(reg_.chromosome(), std::atoll(optarg ? optarg : ""), reg_.to());
          break;
        }
        else if (std::string(long_options_[long_index].name) == "end")
        {
          std::cerr << "Warning: --end is deprecated in favor of --region\n";
          reg_ = savvy::genomic_region(reg_.chromosome(), reg_.from(), std::atoll(optarg ? optarg : ""));
          break;
        }
        else if (std::string(long_options_[long_index].name) == "prefix")
        {
          std::cerr << "Warning: --prefix is deprecated in favor of --output, --empirical-output, and --sites\n";
          prefix_ = optarg ? optarg : "";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "mapFile")
        {
          std::cerr << "Warning: --mapFile is deprecated in favor of --map\n";
          map_path_ = optarg ? optarg : "";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "window")
        {
          std::cerr << "Warning: --window is deprecated in favor of --overlap\n";
          overlap_ = std::atoll(optarg ? optarg : "");
          break;
        }
        else if (std::string(long_options_[long_index].name) == "ChunkOverlapMb")
        {
          std::cerr << "Warning: --ChunkOverlapMb is deprecated in favor of --overlap\n";
          overlap_ = std::atoll(optarg ? optarg : "") * 1000000;
          break;
        }
        // else pass through to default
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
    else if (!prefix_.empty())
    {
      std::string suffix = "sav";
      if (out_format_ == savvy::file::format::bcf)
        suffix = "bcf";
      else if (out_format_ == savvy::file::format::vcf)
      {
        suffix = "vcf";
        if (out_compression_)
          suffix += ".gz";
      }

      out_path_ = prefix_ + ".dose." + suffix;
      out_path_ = prefix_ + ".sites." + suffix;
      if (meta_)
        emp_out_path_ = prefix_ + ".empiricalDose." + suffix;
    }
    else if (remaining_arg_count < 2)
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

bool impute_chunk(const savvy::region& impute_region, const prog_args& args, dosage_writer& output)
{
  savvy::region extended_region =
    {
      impute_region.chromosome(),
      std::uint64_t(std::max(std::int64_t(1), std::int64_t(impute_region.from()) - args.overlap())),
      impute_region.to() + args.overlap()
    };

  std::cerr << "Imputing " << impute_region.chromosome() << ":" << impute_region.from() << "-" << impute_region.to() << std::endl;

  std::time_t start_time = std::time(nullptr);
  std::vector<std::string> sample_ids;
  std::vector<target_variant> target_sites;
  load_target_haplotypes(args.tar_path(), extended_region, target_sites, sample_ids);
  std::cerr << ("Loading target haplotypes took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  start_time = std::time(nullptr);

  reduced_haplotypes typed_only_reference_data(16, 512);
  reduced_haplotypes full_reference_data;
  load_reference_haplotypes(args.ref_path(), extended_region, impute_region, target_sites, typed_only_reference_data, full_reference_data);
  std::cerr << ("Loading reference haplotypes took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  start_time = std::time(nullptr);
  std::vector<target_variant> target_only_sites = separate_target_only_variants(target_sites);
  if (!args.all_typed_sites())
    target_only_sites.clear();
  std::cerr << ("Separating typed only variants took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;

  if (target_sites.empty())
    return std::cerr << "Error: no target variants\n", false;

  target_sites.back().recom = 0.f; // Last recom prob must be zero so that the first step of backward traversal will have no recombination.
  if (args.map_path().size())
  {
    start_time = std::time(nullptr);
    if (!recombination::parse_map_file(args.map_path(), target_sites.begin(), target_sites.end()))
      return std::cerr << "Error: parsing map file failed\n", false;
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

  //  std::list<savvy::reader> temp_files;
  //  std::list<savvy::reader> temp_emp_files;
  std::list<std::string> temp_files;
  std::list<std::string> temp_emp_files;

  full_dosages_results hmm_results;
  hmm_results.resize(full_reference_data.variant_size(), target_sites.size(), std::min(haplotype_buffer_size, target_sites[0].gt.size()));

  double impute_time = 0.;
  double temp_write_time = 0.;
  bool use_temp_files = target_sites[0].gt.size() > haplotype_buffer_size;
  for (std::size_t i = 0; i < target_sites[0].gt.size(); i += haplotype_buffer_size)
  {
    std::size_t group_size = std::min(target_sites[0].gt.size() - i, haplotype_buffer_size);
    if (group_size < haplotype_buffer_size)
    {
      for (std::size_t j = 0; j < hmm_results.dosages_.size(); ++j)
        hmm_results.dosages_[j].resize(group_size);
      for (std::size_t j = 0; j < hmm_results.loo_dosages_.size(); ++j)
        hmm_results.loo_dosages_[j].resize(group_size);
    }

    start_time = std::time(nullptr);
    omp::parallel_for_exp(omp::static_schedule(), omp::sequence_iterator(i), omp::sequence_iterator(i + group_size), [&](int& i, const omp::iteration_context& ctx)
      {
        hmms[ctx.thread_index].traverse_forward(typed_only_reference_data.blocks(), target_sites, i);
        hmms[ctx.thread_index].traverse_backward(typed_only_reference_data.blocks(), target_sites, i, i % haplotype_buffer_size, reverse_maps, hmm_results, full_reference_data);
      },
      tpool);
    impute_time += std::difftime(std::time(nullptr), start_time);

    start_time = std::time(nullptr);

    int tmp_fd = -1;
    int tmp_emp_fd = -1;

    if (use_temp_files)
    {
      std::string out_emp_path;
      std::string out_path = "/tmp/m4_" + std::to_string(i / haplotype_buffer_size) + "_XXXXXX";
      tmp_fd = mkstemp(&out_path[0]);
      if (tmp_fd < 0)
        return std::cerr << "Error: could not open temp file (" << out_path << ")" << std::endl, false;

      if (args.emp_out_path().size())
      {
        out_emp_path = "/tmp/m4_" + std::to_string(i / haplotype_buffer_size) + "_emp_XXXXXX";
        tmp_emp_fd = mkstemp(&out_emp_path[0]);
        if (tmp_emp_fd < 0)
          return std::cerr << "Error: could not open temp file (" << out_emp_path << ")" << std::endl, false;
      }

      dosage_writer temp_output(out_path, out_emp_path,
        "", // sites path
        savvy::file::format::sav,
        std::min<std::uint8_t>(3, args.out_compression()),
        {sample_ids.begin() + (i / ploidy), sample_ids.begin() + (i + group_size) / ploidy},
        {"HDS"},
        impute_region.chromosome(),
        use_temp_files);

      temp_files.emplace_back(out_path);
      ::close(tmp_fd);
      // std::remove(out_path.c_str()); // TODO: update savvy to indices to use FILE*
      assert(tmp_fd > 0);

      if (out_emp_path.size())
      {
        temp_emp_files.emplace_back(out_emp_path);
        ::close(tmp_emp_fd);
        // std::remove(out_emp_path.c_str()); // TODO: update savvy to indices to use FILE*
        assert(tmp_emp_fd > 0);
      }

      if (!temp_output.write_dosages(hmm_results, target_sites, target_only_sites, {i, i + group_size}, full_reference_data, impute_region))
        return std::cerr << "Error: failed writing output\n", false;
      temp_write_time += std::difftime(std::time(nullptr), start_time);
    }
  }

  std::cerr << ("Running HMM took " + std::to_string(impute_time) + " seconds") << std::endl;
  if (use_temp_files)
  {
    std::cerr << ("Writing temp files took " + std::to_string(temp_write_time) + " seconds") << std::endl;

    std::cerr << "Merging temp files ... " << std::endl;
    start_time = std::time(nullptr);
    //dosage_writer output(args.out_path(), args.emp_out_path(), args.sites_out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom, false);
    if (!output.merge_temp_files(temp_files, temp_emp_files))
      return std::cerr << "Error: failed merging temp files\n", false;
    std::cerr << ("Merging temp files took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }
  else
  {
    std::cerr << "Writing output ... " << std::endl;
    start_time = std::time(nullptr);
    if (!output.write_dosages(hmm_results, target_sites, target_only_sites, {0, target_sites[0].gt.size()}, full_reference_data, impute_region))
      return std::cerr << "Error: failed writing output\n", false;
    std::cerr << ("Writing output took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }

  return true;
}

int main(int argc, char** argv)
{
//  test_class t(test_class::file_ptr_t(fopen("","")));
//  test_class t2(fopen("",""));

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

  std::string chrom = args.region().chromosome();
  std::vector<std::string> sample_ids;

  {
    savvy::reader temp_rdr(args.tar_path());
    savvy::variant first_var;
    if (!temp_rdr || !temp_rdr.read(first_var))
      return std::cerr << "Error: could not open target file\n", EXIT_FAILURE;

    sample_ids = temp_rdr.samples();
    if (chrom.empty())
      chrom = first_var.chromosome();
  }

  dosage_writer output(args.out_path(),
    args.emp_out_path(),
    args.sites_out_path(),
    args.out_format(),
    args.out_compression(),
    sample_ids,
    args.fmt_fields(),
    chrom,
    false);

  for (std::uint64_t chunk_start_pos = std::max(std::uint64_t(1), args.region().from()); chunk_start_pos < args.region().to(); chunk_start_pos += args.chunk_size())
  {
    std::uint64_t chunk_end_pos = std::min(args.region().to(), chunk_start_pos + args.chunk_size() - 1ul);
    savvy::region impute_region =
      {
        chrom,
        chunk_start_pos,
        chunk_end_pos
      };

    if (!impute_chunk(impute_region, args, output))
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
