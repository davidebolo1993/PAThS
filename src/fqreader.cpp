// compile with g++ -o fqreader fqreader.cpp -lz -std=c++11

#include <iostream>
#include <string>
#include <zlib.h>
#include "stdc++.h" 

#include "kseq.h"


KSEQ_INIT(gzFile,gzread)

int main(int argc, char **argv)
{

  kseq_t *seq;
  int n = 0;
  int kmer = 3;
  std::unordered_set<std::string> set;

  gzFile fp;
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);
  while (kseq_read(seq) >= 0){
    ++n;

    std::string s = seq->seq.s;
    std::string q = seq->qual.s;

    for (int i = 0; i <= s.length() - kmer; i++) {

      std::string subs = s.substr(i, kmer);
      set.insert(subs);

    }

    std::cout << s << std::endl << q << std::endl;

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