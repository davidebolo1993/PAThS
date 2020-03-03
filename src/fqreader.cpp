// compile with g++ -o fqreader fqreader.cpp -lz -std=c++11

#include <iostream>
#include <string>
#include <zlib.h>

#include "stdc++.h" 
#include "kseq.h"

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
      std::cout << forw << std::endl;
      std::string rev = forw;
      std::transform(rev.rbegin(), rev.rend(), rev.begin(), complement);
      set.insert(rev);
      std::cout << rev << std::endl;

    }

  }

  for (auto s : set)
  
  {
   
    std::cout << s << " "; 
  
  }

  std::cout << std::endl << "Processed " << n << " sequences" << '\n';

  kseq_destroy(seq);
  gzclose(fp);
  

  return 0;
}