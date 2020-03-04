#ifndef KMER_H
#define KMER_H

#include <iostream>
#include <string>
#include <zlib.h>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>

#include "stdc++.h" 
#include "kseq.h"


#include <sdsl/suffix_arrays.hpp>
using namespace sdsl;

char complement(char n)
{   
    switch(n)
    {   
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'N':
      return 'N';  
    default:
      throw std::domain_error("Invalid nucleotide.");
    }
}   



struct Container {
    
  int kmerlength;
  boost::filesystem::path outfile;
  boost::filesystem::path infile;
  boost::filesystem::path dumpfile;
};



KSEQ_INIT(gzFile,gzread)

int kmers(int argc, char **argv)
{


  Container c;
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
  ("help,?", "show help message")
  ("dump,u", boost::program_options::value<boost::filesystem::path>(&c.dumpfile)->default_value("kmers.txt"), "output file of k-mers")
  ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("kmers.fm"), "output FM-index")
  ("kmer,k", boost::program_options::value<int>(&c.kmerlength)->default_value(27), "k-mer length")
  ;
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
  ("input,i", boost::program_options::value<boost::filesystem::path>(&c.infile), "input FASTQ file")
  ;
  boost::program_options::positional_options_description pos_args;
  pos_args.add("input", -1);
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  if ((vm.count("help")) || (!vm.count("input"))) {

    std::cout << std::endl;
    std::cout << "Usage: paths " << argv[0] << " [OPTIONS] <input.fq/input.fq.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  
  }

  csa_wt<> fmi;
  kseq_t *seq;
  int n = 0;
  int kmer = c.kmerlength;


  std::unordered_set<std::string> set;
  gzFile fp;
  fp = gzopen(c.infile.string().c_str(), "r");
  seq = kseq_init(fp);
  while (kseq_read(seq) >= 0){
    
    ++n;
    std::string s = seq->seq.s;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    
    //do something with qualities?
    //std::string q = seq->qual.s;

    for (int i = 0; i <= s.length() - kmer; i++) {

      std::string forw = s.substr(i, kmer);
      set.insert(forw);
      std::string rev = forw;
      std::transform(rev.rbegin(), rev.rend(), rev.begin(), complement);
      set.insert(rev);
    
    }

  }

  kseq_destroy(seq);
  gzclose(fp);
  std::ofstream kmdump(c.dumpfile.string().c_str()); 

  for (auto s : set)
  
  {
   
    kmdump << s << std::endl;
  
  }

  kmdump.close();
  std::cout << "Processed " << n << " sequences" << '\n';
  construct(fmi,c.dumpfile.string().c_str(), 1);
  store_to_file(fmi,c.outfile.string().c_str());

  return 0;
}

#endif
