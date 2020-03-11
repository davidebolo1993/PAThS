#ifndef KMER_H
#define KMER_H

//system libraries

#include <iostream>
#include <string>
#include <fstream>
#include <zlib.h>
#include <thread>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

//header files in folder

#include "stdc++.h" 
#include "kseq.h"

//header-only sparsepp library

#include "sparsepp/sparsepp/spp.h"
using spp::sparse_hash_map;

struct Container {    
  int kmerlength;
  int minquality;
  boost::filesystem::path infile;
  boost::filesystem::path infilelist;
  boost::filesystem::path jsonfile;
  boost::filesystem::path mapfile;
};

class FileSerializer 
{
public:
    // serialize basic types to FILE
    // -----------------------------
    template <class T>
    bool operator()(FILE *fp, const T& value) 
    {
        return fwrite((const void *)&value, sizeof(value), 1, fp) == 1;
    }

    template <class T>
    bool operator()(FILE *fp, T* value) 
    {
        return fread((void *)value, sizeof(*value), 1, fp) == 1;
    }

    // serialize std::string to FILE
    // -----------------------------
    bool operator()(FILE *fp, const std::string& value) 
    {
        const size_t size = value.size();
        return (*this)(fp, size) && fwrite(value.c_str(), size, 1, fp) == 1;
    }

    bool operator()(FILE *fp, std::string* value) 
    {
        size_t size;
        if (!(*this)(fp, &size)) 
            return false;
        char* buf = new char[size];
        if (fread(buf, size, 1, fp) != 1) 
        {
            delete [] buf;
            return false;
        }
        new (value) std::string(buf, (size_t)size);
        delete[] buf;
        return true;
    }

    // serialize std::pair<const A, B> to FILE - needed for maps
    // ---------------------------------------------------------
    template <class A, class B>
    bool operator()(FILE *fp, const std::pair<const A, B>& value)
    {
        return (*this)(fp, value.first) && (*this)(fp, value.second);
    }

    template <class A, class B>
    bool operator()(FILE *fp, std::pair<const A, B> *value) 
    {
        return (*this)(fp, (A *)&value->first) && (*this)(fp, &value->second);
    }
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
    //case 'N': //excluding kmers with N makes this un-necessary
      //return 'N';  
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

unsigned djb2_hash(const char *str) { //excellent hash string function
  
  unsigned hash = 5381;
  int c;

  while (c = *str++)
        
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

  return hash;
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
  ("quality,q", boost::program_options::value<int>(&c.minquality)->default_value(30), "minimum average qscore to retain k-mers (only works for FASTQ)")
  ;

  boost::program_options::options_description input("Input options");
  input.add_options()
  ("list,l", boost::program_options::value<boost::filesystem::path>(&c.infilelist), "multiple FASTQ/FASTA (one per line)")
  ;

  boost::program_options::options_description output("Output options");
  output.add_options()
  ("json,j", boost::program_options::value<boost::filesystem::path>(&c.jsonfile)->default_value("kmers.json.gz"), "output k-mers spectra in JSON")
  ("map,m", boost::program_options::value<boost::filesystem::path>(&c.mapfile)->default_value("kmers.map"), "output binary k-mers hash map")
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

  sparse_hash_map<unsigned, int> hashmap;
  sparse_hash_map<int,int> khash; //this stores kmers pectra
  kseq_t *seq;
  gzFile fp;
  int n = 0;
  int klen,avgqual;
  std::string s,q,forw,qual;
  unsigned keyfw,keyrc;
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

  if (vm.count("list")) {

    //open multiple input files (can be gzipped or not), one per time, and store unique k-mers for all of them in the same hash map

    std::ifstream inlist(c.infilelist.string().c_str());

    if (inlist.is_open()) {

      std::string line;
      
      while (std::getline(inlist,line)) {

        fp = gzopen(line.c_str(), "rb");
        now = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing \"" <<  line.c_str() << "\"" << std::endl;
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

            if (forw.find('N') != std::string::npos) {

              continue;

            }

            if (seq->is_fastq) {

              qual = q.substr(i,klen);
              avgqual = avgq(qual);
              
              if (avgqual < c.minquality) {

                continue;

              }

            }

            keyfw=djb2_hash(forw.c_str());
            hashmap[keyfw] ++;
            //std::cout << forw << ":"<< keyfw << std::endl;
            std::reverse(forw.begin(),forw.end());
            std::transform(forw.begin(), forw.end(), forw.begin(), complement);
            keyrc=djb2_hash(forw.c_str());
            hashmap[keyrc] ++;
            //std::cout << forw << ":"<< keyrc << std::endl;

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
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing \"" <<  c.infile.string().c_str() << "\"" << std::endl;
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

        if (forw.find('N') != std::string::npos) {

          continue;

        }
        
        if (seq->is_fastq) {

          qual = q.substr(i,klen);
          avgqual = avgq(qual);

          if (avgqual < c.minquality) {

            continue;

          }

        }

        keyfw=djb2_hash(forw.c_str());
        hashmap[keyfw] ++;
        std::cout << forw << ":"<< keyfw << std::endl;
        std::reverse(forw.begin(),forw.end());
        std::transform(forw.begin(), forw.end(), forw.begin(), complement);
        keyrc=djb2_hash(forw.c_str());
        hashmap[keyrc] ++;
        std::cout << forw << ":"<< keyrc << std::endl;
      }

    }

    kseq_destroy(seq);
    gzclose(fp);

  }

  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done" << std::endl;  
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Found " << hashmap.size() << " unique k-mers in " << n << " sequences" << std::endl;
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Storing k-mers spectra \"" << c.jsonfile.string().c_str() << "\"" << std::endl;

  for (auto it = hashmap.cbegin(); it != hashmap.cend(); ++it) {

    khash[(*it).second] ++;

  }

  boost::iostreams::filtering_ostream kmjson;
  kmjson.push(boost::iostreams::gzip_compressor());
  kmjson.push(boost::iostreams::file_sink(c.jsonfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
  kmjson << "{" << std::endl;
  std::string delim = "";
  
  for (auto it = khash.cbegin(); it != khash.cend(); ++it) {

    kmjson << delim << "\"" << (*it).first << "\": " << (*it).second;
    delim = ",\n";

  }

  kmjson << std::endl << "}" << std::endl;
  kmjson.pop();

  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done" << std::endl;
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Storing binary k-mers hash map \"" << c.mapfile.string().c_str() << "\"" << std::endl;

  FILE *out = fopen(c.mapfile.string().c_str(), "wb");
  hashmap.serialize(FileSerializer(), out);
  fclose(out);

  //recheck

  //sparse_hash_map<unsigned, int> hashtest;
  //FILE *in = fopen(c.mapfile.string().c_str(), "rb");
  //hashtest.unserialize(FileSerializer(), in);
  //fclose(in);

  //for (auto it = hashtest.cbegin(); it != hashtest.cend(); ++it) {

    //std::cout << (*it).first << ":" << (*it).second << std::endl;

  //}

  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done" << std::endl;
  
  return 0;

}

#endif
