
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
  input.reset_bounds(reg);
  savvy::variant var;
  std::vector<std::int8_t> tmp_geno;
  while (input >> var)
  {
    var.get_format("GT", tmp_geno);
    for (std::size_t i = 0; i < var.alts().size(); ++i)
    {
      std::size_t allele_idx = i + 1;
      target_sites.push_back({var.chromosome(), var.position(), var.ref(), var.alts()[i], true, false, std::numeric_limits<float>::quiet_NaN(), 0.01, recombination::recom_min, {}});
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

bool load_reference_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, reduced_haplotypes& reference_data)
{
  savvy::reader input(file_path);
  input.reset_bounds(reg);
  savvy::variant var;
  std::vector<std::int8_t> tmp_geno;

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
        reference_data.compress_variant({it->chrom, it->pos, it->ref, it->alt}, tmp_geno);
        //freq.push_back(std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size());
        //it->af = std::accumulate(tmp_geno.begin(), tmp_geno.end(), 0.f) / tmp_geno.size(); // TODO; remove
        it->af = float((--reference_data.end())->ac) / tmp_geno.size();
        it->in_ref = true;
        if (it != tar_it)
          std::swap(*it, *tar_it);
        ++tar_it;
        break;
      }
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
}

class haplotype_interpolator
{
private:
  savvy::writer out_file_;
  savvy::file::format file_format_;
  std::vector<std::string> fmt_fields_;
  std::size_t n_samples_ = 0;
  const std::int16_t bin_scalar_ = 100; //256;
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

  static bool sites_match(const target_variant& t, const savvy::variant& r)
  {
    bool ret = t.pos == r.position() && t.alt == r.alts()[0] && t.ref == r.ref();
    return ret;
    //return typed_variants[global_idx].pos == ref_var.position() && typed_variants[global_idx].alt == ref_var.alts()[0] && typed_variants[global_idx].ref == ref_var.ref();
  }

  bool piecewise_constant(const std::vector<target_variant>& typed_variants, const reduced_haplotypes& reduced_reference, const best_templates_results& hmm_results, const std::string& reference_path, const savvy::genomic_region& reg, const std::string& output_path)
  {
    auto dims = hmm_results.dimensions();
    savvy::reader ref_file(reference_path);
    savvy::variant ref_var, out_var;
    ref_file.reset_bounds(reg);

    //savvy::writer out_file(output_path, savvy::file::format::bcf, {}, {});

    std::vector<float> dosages(dims[1]);
    std::vector<float> loo_dosages(dims[1]);
    const std::vector<float> empty_loo_vec;
    const std::vector<std::int8_t> empty_gt_vec;

    savvy::compressed_vector<std::int8_t> genos;
    std::vector<std::size_t> allele_counts;
    std::vector<std::int8_t> dense_genos;

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

          prepare_output_variant(out_var, ref_var, dosages, ref_matches_typed ? loo_dosages : empty_loo_vec, ref_matches_typed ? typed_variants[global_idx].gt : empty_gt_vec);
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
  }
private:
  void prepare_output_variant(savvy::variant& out_var, const savvy::variant& ref_var, const std::vector<float>& dosages, const std::vector<float>& loo_dosages, const std::vector<std::int8_t>& observed)
  {
    savvy::compressed_vector<std::int8_t> sparse_gt;
    savvy::compressed_vector<float> sparse_dosages;
    std::vector<float> gp_vec;
    sparse_dosages.assign(dosages.begin(), dosages.end());
    out_var = savvy::site_info(ref_var.chromosome(), ref_var.position(), ref_var.ref(), ref_var.alts(), ref_var.id());

    std::size_t n = sparse_dosages.size();
    float s_x = std::accumulate(sparse_dosages.begin(), sparse_dosages.end(), 0.f);
    float s_xx = std::inner_product(sparse_dosages.begin(), sparse_dosages.end(), sparse_dosages.begin(), 0.f);
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
          sparse_dosages.assign(dosages.begin(), dosages.end());
        out_var.set_format("HDS", sparse_dosages);
        hds_set = true;
      }
      else if (fmt == "GT")
      {
        sparse_gt.assign(dosages.begin(), dosages.end(), [](float v){ return std::int8_t(v < 0.5f ? 0 : 1); });
        out_var.set_format("GT", sparse_gt);
      }
      else if (fmt == "DS")
      {
        if (!hds_set)
          sparse_dosages.assign(dosages.begin(), dosages.end());
        savvy::stride_reduce(sparse_dosages, sparse_dosages.size() / n_samples_);
        out_var.set_format("DS", sparse_dosages);
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
      sparse_dosages.assign(loo_dosages.begin(), loo_dosages.end());
      sparse_gt.assign(observed.begin(), observed.end());
      s_x = std::accumulate(sparse_dosages.begin(), sparse_dosages.end(), 0.f);
      s_xx = std::inner_product(sparse_dosages.begin(), sparse_dosages.end(), sparse_dosages.begin(), 0.f);
      float s_y = std::accumulate(sparse_gt.begin(), sparse_gt.end(), 0.f);
      // since observed can only be 0 or 1, s_yy is the same as s_y
      float s_yy = s_y; //std::inner_product(sparse_gt.begin(), sparse_gt.end(), sparse_gt.begin(), 0.f); // TODO: allow for missing oberserved genotypes.
      float s_xy = 0.f;
      for (auto it = sparse_gt.begin(); it != sparse_gt.end(); ++it)
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

int main(int argc, char** argv)
{
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

  reduced_haplotypes reduced_reference_data(16, 512);
  load_reference_haplotypes(args.ref_path(), extended_region, target_sites, reduced_reference_data);
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
    if (!recombination::parse_map_file(args.map_path(), target_sites))
      return std::cerr << "Error: parsing map file failed\n", EXIT_FAILURE;
    std::cerr << ("Loading switch probabilities took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }

  std::vector<std::vector<std::vector<std::size_t>>> reverse_maps;
  reverse_maps.reserve(reduced_reference_data.blocks().size());
  for (auto it = reduced_reference_data.blocks().begin(); it != reduced_reference_data.blocks().end(); ++it)
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

  bool deferred_interoplation = true;
  if (deferred_interoplation)
  {
    best_templates_results hmm_results;
    hmm_results.resize(target_sites.size(), target_sites[0].gt.size());

    omp::parallel_for_exp(omp::static_schedule(), omp::sequence_iterator(0), omp::sequence_iterator(target_sites[0].gt.size()), [&](int& i, const omp::iteration_context& ctx)
      {
        hmms[ctx.thread_index].traverse_forward(reduced_reference_data.blocks(), target_sites, i);
        hmms[ctx.thread_index].traverse_backward(reduced_reference_data.blocks(), target_sites, i, reverse_maps, hmm_results);
      },
      tpool);
    std::cerr << ("Running HMM took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
    start_time = std::time(nullptr);
    haplotype_interpolator interpolator(args.out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom);
    interpolator.piecewise_constant(target_sites, reduced_reference_data, hmm_results, args.ref_path(), args.region(), args.out_path());
    std::cerr << ("Interpolating and writing output took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }
  else
  {
    full_dosages_results hmm_results;
    reduced_haplotypes full_reference_data(16, 512); // TODO: Load from m3vcf file
    hmm_results.resize(full_reference_data.variant_size(), target_sites[0].gt.size());

    omp::parallel_for_exp(omp::static_schedule(), omp::sequence_iterator(0), omp::sequence_iterator(target_sites[0].gt.size()), [&](int& i, const omp::iteration_context& ctx)
      {
        hmms[ctx.thread_index].traverse_forward(reduced_reference_data.blocks(), target_sites, i);
        auto reverse_ref_itr = --full_reference_data.end();
        hmms[ctx.thread_index].traverse_backward(reduced_reference_data.blocks(), target_sites, i, reverse_maps, hmm_results, reverse_ref_itr, --full_reference_data.begin());
      },
      tpool);
    std::cerr << ("Running HMM took " + std::to_string(std::difftime(std::time(nullptr), start_time)) + " seconds") << std::endl;
  }



  return EXIT_SUCCESS;
}
