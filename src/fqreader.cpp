#include <iostream>
#include <string>
#include <zlib.h>
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


KSEQ_INIT(gzFile,gzread)

int main(int argc, char **argv)
{

  csa_wt<> fmi;
  kseq_t *seq;
  int n = 0;
  int kmer = atoi(argv[2]);
  std::unordered_set<std::string> set;
  gzFile fp;
  fp = gzopen(argv[1], "r");
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
  std::ofstream kmdump("kmers.txt"); 

  for (auto s : set)
  
  {
   
    kmdump << s << std::endl;
  
  }

  std::cout << "Processed " << n << " sequences" << '\n';
  construct(fmi,"kmers.txt", 1);
  store_to_file(fmi,"kmers.fm");

  return 0;
}
