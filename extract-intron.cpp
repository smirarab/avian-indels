#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include <cmath>
using namespace std;

#define MAXSEQ 100000

int main(int argc, char* argv[]) {

	if (argc != 5) {
		cerr << "usage: ./extract-intron <infile> <start> <end> <length>\n";
		cerr << "  infile = aligned fasta file\n";
		cerr << "  start  = start of partition to extract\n";
		cerr << "  end    = end of partition to extract\n";
		cerr << "  length = length of the sequence alignment\n";
		return 0;
	}
	
	unsigned start = atoi(argv[2]) - 1;
	unsigned end = atoi(argv[3]) + 1;
	unsigned length = atoi(argv[4]);
	
	unsigned i, pos, nseq, allmissing;
	char nt;
	char name[100];
	char seq[MAXSEQ];
	
	ifstream inf(argv[1]); // argv[1] is the infile name

	nt='x';
	while ( nt != '\0' ) {
		nt='\0';
		inf >> nt;
		if ( nt == '>' ) {
			inf >> name;
		//	cout << '>' << name << endl;
			seq[0]='\0';
			pos=1; nseq=0; allmissing=1;
			for (i=0; i<length; i++) {
				inf >> nt;
				if ( pos > start && pos < end ) { 
				//	cout << nt; 
					seq[nseq]=nt;
					// convert missing data to N's
					if ( nt == '?' ) { seq[nseq]='N'; }
					nseq++;
					seq[nseq]='\0';
					if ( nt != '?' ) { allmissing=0; }
				}
				pos++;
			}
			// only output sequence if there is some defined data
			if ( allmissing == 0 ) {
				cout << '>' << name << endl;
				cout << seq << endl;
			}
		}
	}
	
	inf.close();
	
	return 0;
}
