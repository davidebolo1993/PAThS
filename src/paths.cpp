#include <iostream>
#include <fstream>

#include "kmer.h"

inline void displayUsage() {
  std::cout << "PAThS (https://github.com/davidebolo1993/PAThS)" << std::endl;
  std::cout << std::endl;
  std::cout << "Usage: paths <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Commands:" << std::endl;
  std::cout << std::endl;
  std::cout << "    kmers       generate a FM index for k-mers in FASTQ/FASTA" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

inline void asciiArt() {

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout <<" ________  ________  _________  ___       ________  " << std::endl;
  std::cout <<"|\\   __  \\|\\   __  \\|\\___   ___\\\\  \\     |\\   ____\\   "  << std::endl;
  std::cout <<"\\ \\  \\|\\  \\ \\  \\|\\  \\|___ \\  \\_\\ \\  \\____\\ \\  \\___|_   " << std::endl;
  std::cout <<" \\ \\   ____\\ \\   __  \\   \\ \\  \\ \\ \\   __  \\ \\_____  \\   " << std::endl;
  std::cout <<"  \\ \\  \\___|\\ \\  \\ \\  \\   \\ \\  \\ \\ \\  \\ \\  \\|____|\\  \\  " << std::endl;
  std::cout <<"   \\ \\__\\    \\ \\__\\ \\__\\   \\ \\__\\ \\ \\__\\ \\__\\____\\_\\  \\ "<< std::endl;
  std::cout <<"    \\|__|     \\|__|\\|__|    \\|__|  \\|__|\\|__|\\_________\\" << std::endl;
  std::cout <<"                                            \\|_________|" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}


int main(int argc, char **argv) {
  
  std::string pathsversion = "0.1 (dirty)";
  if (argc < 2) { 
    asciiArt();
    displayUsage();
    return 0;
  }  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    asciiArt();
    std::cout << "paths version: v" << pathsversion << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    asciiArt();
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "kmers")) {
    return kmers(argc-1,argv+1);
  } else {
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
  }
}