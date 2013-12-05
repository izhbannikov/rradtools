#include <stdlib.h>     /* malloc, free, rand */
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

//Структура содержащая правильный нуклеотид и позицию в чтении:
struct Corr {
 char base; // base that replaces 'wrong' nucleotide
 int pos; //position where that 'wrong' nucleotide occurs
};



std::vector<char> GetBases(char b) {
	std::vector<char> bases;
	bases.push_back('G'); bases.push_back('C'); bases.push_back('T'); 
	if(b == 'A') {
		return bases;
	} else if(b == 'G') {
		bases.clear();
		bases.push_back('A'); bases.push_back('C'); bases.push_back('T'); 
	} else if(b == 'C') {
		bases.clear();
		bases.push_back('A'); bases.push_back('G'); bases.push_back('T'); 
	} else if(b == 'T') {
		bases.clear();
		bases.push_back('A'); bases.push_back('G'); bases.push_back('C'); 
	} 

	return bases;
	
}

int main(int argc, char *argv[]) {

	
	std::map<std::string, std::map<int, char> > reads;
	std::map< std::string, std::vector<Read> > kmers;
	std::map< std::string, int > kmer_freqs;
	std::map< std::string, int >::iterator it_kmer_freqs;
	std::map< std::string, std::vector<Read> >::iterator it_kmers;
		
	int ii = 0;
        int k = 30;
	int rad_length = 50;
	std::string line;
       	std::ifstream in("/Users/ilya/bio/app/RAD/first4k.fastq");
       	std::vector<std::string> record_block;
		
        //Первый прогон: разбиваем чтения на кмеры:
       	while ( getline(in,line) )
       	{
                	//Read ID
                	if(ii==0) 
                	{
                	    record_block.push_back(line); 
			    std::map<int, char> corr; corr[-1] = "";
			    reads[line] = corr;
                            ii++;
                	    continue;
                	}
                	//DNA string
                	if(ii==1) 
                	{
                	    record_block.push_back(line.substr(12,50));
                	    ii++;
                	    continue;
                	}
                	//a symbol "+"
                	if(ii==2) 
                	{
                	     record_block.push_back(line);
                	     ii++;
                	     continue;
                	}
                	if(ii==3) 
                	{
                	     ii=0;
           
                	     record_block.push_back(line); //Сохраняем качества чтения нуклеотидов
                	     
			     for(int i = 0; i<record_block[1].length() - k + 1; ++i) {
				Read r;
				r.r_id = record_block[0]; r.pos = i;
				kmer_freqs[record_block[1].substr(i,k)].push_back(r);
				kmers[record_block[1].substr(i,k)] += 1;
			     }

			     
			     record_block.clear();
		     
                	     
                	}	
	}


	std::cout << "Before correction: " << kmers.size() << std::endl;

	//Второй прогон: корректировка чтений:
	//Рассматриваем только те кмеры, у которых частота = 1:
	for(it_kmer_freqs=kmer_freqs.begin(); it_kmer_freqs != kmer_freqs.end(); ++it_kmer_freqs) {
		std::string kmer = (*it_kmer_freqs).first;
		if( (*it_kmer_freqs).second == 1 ) {
			//Начинаем корректировать:
			for(int i=0; i<k; ++i) {
				//Меняем каждый нуклеотид в k-mere:
				if(kmer[i] != 'N') {
					char base = kmer[i];
					std::vector<char> bases = GetBases(base);

					kmer[i] = bases[0];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second > 1) ) {
						kmers[kmer] += 1;
						//Удаляем предположительно ошибочный кмер:
						kmer[i] = base;
						kmers.erase(kmer);
						break;
					}

					kmer[i] = bases[1];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second > 1) ) {
						kmers[kmer] += 1;
						//Удаляем предположительно ошибочный кмер:
						kmer[i] = base;
						kmers.erase(kmer);
						break;
					}

					kmer[i] = bases[2];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second > 1) ) {
						kmers[kmer] += 1;
						//Удаляем предположительно ошибочный кмер:
						kmer[i] = base;
						kmers.erase(kmer);
						break;
					}
				}	
			}			
			//it_kmer = kmers.find(kmer);
			
		}
	}

	//Теперь корректируем нуклеотиды и выводим скорректированный файл:
	

	std::cout << "After correction:" << kmers.size() << std::endl;
	

}
