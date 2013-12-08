#include <stdlib.h>     /* malloc, free, rand */
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

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
	std::map<std::string, std::map<int, char> >::iterator it_reads;
	
	std::map< std::string, std::map< std::string, std::map<int,int> > > kmers;
	std::map< std::string, std::map< std::string, std::map<int,int> > >::iterator it_kmers;
	
	std::map< std::string, int > kmer_freqs;
	std::map< std::string, int >::iterator it_kmer_freqs;
	
	
		
	int ii = 0;
        int k = 15;
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
			    reads[line][-1] = 'N';
                            ii++;
                	    continue;
                	}
                	//DNA string
                	if(ii==1) 
                	{
                	    //record_block.push_back(line.substr(12,line.length()-12));
			    record_block.push_back(line.substr(12,60));
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
                	     
			     for(int i = 12; i<record_block[1].length() - k + 1; ++i) {

				kmers[ record_block[1].substr(i,k) ][ record_block[0] ][ i ] += 1;
				
				kmer_freqs[ record_block[1].substr(i,k) ] += 1;
				
			     }

			     
			     record_block.clear();
		     
                	     
                	}	
	}

	in.close();

	std::cout << "Before correction: " << kmers.size() << std::endl;

	//Второй прогон: корректировка чтений:
	//Рассматриваем только те кмеры, у которых частота = 1:
	for(it_kmer_freqs=kmer_freqs.begin(); it_kmer_freqs != kmer_freqs.end(); ++it_kmer_freqs) {
		std::string kmer = (*it_kmer_freqs).first;
		if( (*it_kmer_freqs).second == 1 ) {
			bool stopflag = false;
			std::string rid = kmers[kmer].begin()->first;
			int pos = kmers[kmer].begin()->second.begin()->first;

			//Начинаем корректировать:
			for(int i=0; i<k; ++i) {
				//Меняем каждый нуклеотид в k-mere:
				if(kmer[i] != 'N') {
					char base = kmer[i];
					std::vector<char> bases = GetBases(base);
					
					kmer[i] = bases[0];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 10) ) {
						//Удаляем предположительно ошибочный кмер:
						reads[rid][12+pos+i] = bases[0];
						std::map<std::string,std::map<int,int> >::iterator _it = it_kmers->second.begin();
						for(_it = it_kmers->second.begin(); _it != it_kmers->second.end(); ++_it) {
							if(	_it->second.find( pos ) != _it->second.end() ) {
								kmers.erase( (*it_kmer_freqs).first ); 
								stopflag = true;break;
							}
						}
						if(stopflag) 
							break;
						//std::cout << rid << ", pos in read " << 12+pos+i << ", base to replace " << base << ", kmer " << (*it_kmer_freqs).first << ", pos in kmer " << i << ", kmer pos: " << (*it_kmers).second[0].pos << '\n';
						
					}
					
					kmer[i] = bases[1];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 10) ) {
						reads[rid][12+pos+i] = bases[1];
						std::map<std::string,std::map<int,int> >::iterator _it = it_kmers->second.begin();
						for(_it = it_kmers->second.begin(); _it != it_kmers->second.end(); ++_it) {
							if(	_it->second.find( pos ) != _it->second.end() ) {
								kmers.erase( (*it_kmer_freqs).first ); 
								stopflag = true;break;
							}
						}
						if(stopflag) 
							break;
					}
					
					kmer[i] = bases[2];
					it_kmers = kmers.find(kmer);
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 10) ) {
						reads[rid][12+pos+i] = bases[2];
						std::map<std::string,std::map<int,int> >::iterator _it = it_kmers->second.begin();
						for(_it = it_kmers->second.begin(); _it != it_kmers->second.end(); ++_it) {
							if(	_it->second.find( pos ) != _it->second.end() ) {
								kmers.erase( (*it_kmer_freqs).first ); 
								stopflag = true;break;
							}
						}
						if(stopflag) 
							break;
					}
				}	
			}			
			
			
		}
	}

	
	//Теперь корректируем нуклеотиды и выводим скорректированный файл:
	ii = 0;
        in.open("/Users/ilya/bio/app/RAD/first4k.fastq");
	std::ofstream out("/Users/ilya/bio/app/RAD/first4k_corrected.fastq");

        while ( getline(in,line) )
       	{
                	//Read ID
                	if(ii==0) 
                	{
                	    record_block.push_back(line); 
			    ii++;
                	    continue;
                	}
                	//DNA string
                	if(ii==1) 
                	{
                	    record_block.push_back(line);
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
                	     
			     //Now physically correct the bases:
			     it_reads = reads.find(record_block[0]);
			     if(it_reads != reads.end()) {
				std::map<int, char>::iterator it_corr;
				for(it_corr=(*it_reads).second.begin(); it_corr != (*it_reads).second.end(); ++it_corr) {
					record_block[1][(*it_corr).first] = (*it_corr).second;
				}
			     }
			     			     
			     //Write corrected reads into a file:
			     out << record_block[0] << "\n" << 	record_block[1] << "\n" << record_block[2] << "\n" << record_block[3] << "\n";					     

			     record_block.clear();
		     
                	     
                	}	
	}

	in.close();
	out.close();
	
	
	std::cout << "After correction:" << kmers.size() << std::endl;
	

}
