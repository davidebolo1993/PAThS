#include <iostream>
#include "kseq_cpp.h"

int main(const int argc, const char *argv[])
{
	kseq_cpp::kseq_parser kseq("test.fq");
	while (kseq.read_seq() >= 0)
	{
		std::cout << kseq.name << std::endl;
		std::cout << kseq.seq << std::endl;
		if (kseq.qual.size())
			std::cout << "+" << std::endl << kseq.qual << std::endl;
	}
	return 0;
}
