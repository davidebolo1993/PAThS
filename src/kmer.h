#ifndef KMER_H
#define KMER_H

#include <iostream>
#include <string>
#include <fstream>
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

struct Container {    
  int kmerlength;
  int minquality;
  boost::filesystem::path outfile;
  boost::filesystem::path infile;
  boost::filesystem::path dumpfile;
  boost::filesystem::path infilelist;
};

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

int avgq(std::string const& s) {

  int aq = 0;
  
  for (int i = 0; i < s.size(); ++i) {

    aq += (int) s[i]-33;
  
  }

  return aq/s.size();

}


KSEQ_INIT(gzFile,gzread)


int kmers(int argc, char **argv)
{

  // Declare container

  Container c;

  // Parse command line

  boost::program_options::options_description generic("Generic options");

  generic.add_options()
  ("help,?", "show help message")
  ("dump,u", boost::program_options::value<boost::filesystem::path>(&c.dumpfile)->default_value("kmers.txt"), "output unique k-mers list")
  ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("kmers.fm"), "output unique k-mers FM index")
  ("list,l", boost::program_options::value<boost::filesystem::path>(&c.infilelist), "multiple input FASTQ list")
  ("kmer,k", boost::program_options::value<int>(&c.kmerlength)->default_value(27), "k-mers length")
  ("quality,q", boost::program_options::value<int>(&c.minquality)->default_value(30), "minimum average qscore to retain k-mers")
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

  if (vm.count("help")) {

    std::cout << std::endl;
    std::cout << "Usage: paths " << argv[0] << " [OPTIONS] <input.fq/input.fq.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  
  }

  else if (!vm.count("input") && !vm.count("list"))

  {
    
    std::cout << std::endl;
    std::cout << "Missing input file/s" << std::endl;
    std::cout << "Usage: paths " << argv[0] << " [OPTIONS] <input.fq/input.fq.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // Generate FM index of k-mers from FASTQ/FASTA

  csa_wt<> fmi;
  std::unordered_set<std::string> set;
  kseq_t *seq;
  gzFile fp;
  int n = 0;

  if (vm.count("list")) {

    //open multiple input file (can be gzipped or not, can be FASTA or FASTQ), one per time, and store unique k-mers for all of them in the same FM index

    std::ifstream inlist(c.infilelist.string().c_str());

    if (inlist.is_open()) {

      std::string line;
      
      while (std::getline(inlist,line)) {

        fp = gzopen(line.c_str(), "rb");
        seq = kseq_init(fp);

        while (kseq_read(seq) >= 0){

          ++n; //count processed sequences
          std::string s = seq->seq.s;
          std::string q = seq->qual.s;

          std::transform(s.begin(), s.end(), s.begin(), ::toupper);

          for (int i = 0; i <= s.length() - c.kmerlength; i++) {

            std::string forw = s.substr(i, c.kmerlength);
            std::string qual = q.substr(i, c.kmerlength);

            int avgqual = avgq(qual);

            if (avgqual < c.minquality) {

              continue;

            }

            set.insert(forw);
            std::reverse(forw.begin(),forw.end());
            std::transform(forw.begin(), forw.end(), forw.begin(), complement);
            set.insert(forw);

          }

        }

        kseq_destroy(seq);
        gzclose(fp);
      }

      inlist.close();
    }
  }

  else {

    //open single input file (can be gzipped or not, can be FASTA or FASTQ)

    fp = gzopen(c.infile.string().c_str(), "rb");
    seq = kseq_init(fp);
    
    //quickly scan FASTQ/FASTQ

    while (kseq_read(seq) >= 0){
      
      ++n; //count processed sequences
      std::string s = seq->seq.s;
      std::string q = seq->qual.s;

      std::transform(s.begin(), s.end(), s.begin(), ::toupper);
      
      for (int i = 0; i <= s.length() - c.kmerlength; i++) {

        std::string forw = s.substr(i, c.kmerlength);
        std::string qual = q.substr(i, c.kmerlength);

        int avgqual = avgq(qual);

        if (avgqual < c.minquality) {

          continue;

        }

        set.insert(forw);
        std::reverse(forw.begin(),forw.end());
        std::transform(forw.begin(), forw.end(), forw.begin(), complement);
        set.insert(forw);
      
      }

    }

    kseq_destroy(seq);
    gzclose(fp);

  }

  std::cout << "Processed " << n << " sequences" << std::endl;

  // output compressed list of unique k-mers

  std::ofstream kmdump;
  kmdump.open(c.dumpfile.string().c_str());
  n=0;

  for (auto s : set)
  
  {
    
    kmdump << s << std::endl;
    n++;
  
  }

  kmdump.close();

  std::cout << "Found " << n << " unique k-mers" << std::endl;

  // construct FM index

  std::cout << "Constructing FM index of unique k-mers from file" << std::endl;

  construct(fmi,c.dumpfile.string().c_str(), 1);
  store_to_file(fmi,c.outfile.string().c_str());

  std::cout << "Done" << std::endl;

  return 0;

}

#endif
