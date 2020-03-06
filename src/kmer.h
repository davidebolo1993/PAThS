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
  boost::filesystem::path mapfile;
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

  boost::program_options::options_description generic("K-mer options");
  generic.add_options()
  ("help,?", "show help message")
  ("kmer,k", boost::program_options::value<int>(&c.kmerlength)->default_value(61), "k-mers length")
  ("quality,q", boost::program_options::value<int>(&c.minquality)->default_value(30), "minimum average qscore to retain k-mers in FASTQ")
  ;

  boost::program_options::options_description input("Input options");
  input.add_options()
  ("list,l", boost::program_options::value<boost::filesystem::path>(&c.infilelist), "multiple input FASTQ/FASTA list")
  ;

  boost::program_options::options_description output("Output options");
  output.add_options()
  ("dump,u", boost::program_options::value<boost::filesystem::path>(&c.dumpfile)->default_value("kmers.txt"), "output unique k-mers list")
  ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("kmers.fm9"), "output unique k-mers FM index")
  ("map,m", boost::program_options::value<boost::filesystem::path>(&c.mapfile)->default_value("kmers.json"), "output k-mers spectra")
  ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
  ("input,i", boost::program_options::value<boost::filesystem::path>(&c.infile), "input FASTQ file")
  ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input", -1);
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(input).add(output).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(input).add(output);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {

    std::cout << std::endl;
    std::cout << "Usage: paths " << argv[0] << " [OPTIONS] <input.fq/input.fa>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  
  }

  else if (!vm.count("input") && !vm.count("list"))

  {
    
    std::cout << std::endl;
    std::cout << "Missing input file/s" << std::endl;
    std::cout << "Usage: paths " << argv[0] << " [OPTIONS] <input.fq/input.fa>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // Generate FM index of k-mers from FASTQ

  csa_wt<> fmi;
  std::unordered_map<std::string,int> hashmap;
  kseq_t *seq;
  gzFile fp;
  int n = 0;
  int klen,avgqual;
  std::string s,q,forw,qual;

  if (vm.count("list")) {

    //open multiple input files (can be gzipped or not), one per time, and store unique k-mers for all of them in the same FM index

    std::ifstream inlist(c.infilelist.string().c_str());

    if (inlist.is_open()) {

      std::string line;
      
      while (std::getline(inlist,line)) {

        fp = gzopen(line.c_str(), "rb");
        seq = kseq_init(fp);

        while (kseq_read(seq) >= 0){

          ++n; //count processed sequences
          s = seq->seq.s; //extract sequence
          std::transform(s.begin(), s.end(), s.begin(), ::toupper);

          if (seq->is_fastq) q = seq->qual.s; //extract quality (if FASTQ)
          
          if (s.length() < c.kmerlength) { //FASTA can have different sequences length

            klen = s.length();

          }

          else {

            klen = c.kmerlength; 

          }

          for (int i = 0; i <= s.length() - klen; i++) {

            forw = s.substr(i,klen);
            
            if (seq->is_fastq) {

              qual = q.substr(i,klen);
              avgqual = avgq(qual);
              
              if (avgqual < c.minquality) {

                continue;

              }

            }

            hashmap[forw] ++;
            std::reverse(forw.begin(),forw.end());
            std::transform(forw.begin(), forw.end(), forw.begin(), complement);
            hashmap[forw] ++;

          }

        }

        kseq_destroy(seq);
        gzclose(fp);
      }

      inlist.close();
    }
  }

  else {

    //open single input file (can be gzipped or not)

    fp = gzopen(c.infile.string().c_str(), "rb");
    seq = kseq_init(fp);
    
    //quickly scan FASTQ/FASTA

    while (kseq_read(seq) >= 0){
      
      ++n; //count processed sequences
      s = seq->seq.s;
      std::transform(s.begin(), s.end(), s.begin(), ::toupper);

      if (seq->is_fastq) q = seq->qual.s; //extract quality (if FASTQ)

      if (s.length() < c.kmerlength) { //FASTA can have different sequences length

        klen = s.length();

      }

      else {

        klen = c.kmerlength; 

      }

      for (int i = 0; i <= s.length() - klen; i++) {

        forw = s.substr(i, klen);
        
        if (seq->is_fastq) {

          qual = q.substr(i,klen);
          avgqual = avgq(qual);

          if (avgqual < c.minquality) {

            continue;

          }

        }

        hashmap[forw] ++;
        std::reverse(forw.begin(),forw.end());
        std::transform(forw.begin(), forw.end(), forw.begin(), complement);
        hashmap[forw] ++;
      
      }

    }

    kseq_destroy(seq);
    gzclose(fp);

  }

  std::cout << "Processed " << n << " sequences" << std::endl;

  // output uncompressed list of unique k-mers

  std::ofstream kmdump,kmmap;

  kmdump.open(c.dumpfile.string().c_str());
  kmmap.open(c.mapfile.string().c_str());

  kmmap << "{" << std::endl;
  std::string delim = "";
  
  for (auto it = hashmap.cbegin(); it != hashmap.cend(); ++it) {

    kmdump << (*it).first << std::endl;
    kmmap << delim << "\"" << (*it).first << "\": " << (*it).second;
    delim = ",\n";

  }

  kmmap << std::endl << "}" << std::endl;

  kmdump.close();
  kmmap.close();

  std::cout << "Found " << hashmap.size() << " unique k-mers" << std::endl;

  // construct FM index from file

  std::cout << "Constructing FM index of unique k-mers from file" << std::endl;

  construct(fmi,c.dumpfile.string().c_str(), 1);
  store_to_file(fmi,c.outfile.string().c_str());

  std::cout << "Done" << std::endl;

  return 0;

}

#endif
