#include <stdlib.h>     /* malloc, free, rand */
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;

int Hamming(std::string s1, std::string s2, int rad_site_length) {
	int dist = 0;
	for(int i=0; i < rad_site_length; ++i) {
		if(s1[i] != s2[i]) 
			dist += 1;
	}
	
	return dist;
}

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

std::vector<std::string> CorrectReads(std::vector<char*> files, int rad_site_length) {
	std::vector<std::string> corrected_reads;
	
	std::map<std::string, std::map<int, char> > reads;
	std::map<std::string, std::map<int, char> >::iterator it_reads;
	std::map<std::string, std::string > read_data;
	std::map<std::string, std::string >::iterator it_read_data;
	
	std::map< std::string, std::map< std::string, std::map<int,int> > > kmers;
	std::map< std::string, std::map< std::string, std::map<int,int> > >::iterator it_kmers;
	
	std::map< std::string, int > kmer_freqs;
	std::map< std::string, int >::iterator it_kmer_freqs;
	
	int k = 15;
	
	for(int f=0; f < files.size(); ++f) {	
		int ii = 0;
        	std::string line;
        	std::ifstream in(files[f]);
        	std::vector<std::string> record_block;
		
	        std::cout << "Processing file: " << files[f] << "\n";
        	
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
			    record_block.push_back(line.substr(12,rad_site_length));
			    read_data[record_block[0]] = line;
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
                	     
			     for(int i = 12; i< 12+rad_site_length- k + 1; ++i) {

				kmers[ record_block[1].substr(i,k) ][ record_block[0] ][ i ] += 1;
				
				kmer_freqs[ record_block[1].substr(i,k) ] += 1;
				
			     }

			
			     record_block.clear();
		     
                	     
                	}	
		}

		in.close();
	}

	
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
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 5) ) {
						//Удаляем предположительно ошибочный кмер:
						reads[rid][12+pos+i] = bases[0];
						std::map<std::string,std::map<int,int> >::iterator _it = it_kmers->second.begin();
						for(_it = it_kmers->second.begin(); _it != it_kmers->second.end(); ++_it) {
							//std::cout << _it->first << "\n";
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
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 5) ) {
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
					if( (it_kmers != kmers.end()) && ((*it_kmers).second.size() > 5) ) {
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

	
	for( it_read_data = read_data.begin(); it_read_data != read_data.end(); it_read_data++ ) {
                //Now physically correct the bases:
		it_reads = reads.find(it_read_data->first);
		if(it_reads != reads.end()) {
			std::map<int, char>::iterator it_corr;
			for(it_corr=(*it_reads).second.begin(); it_corr != (*it_reads).second.end(); ++it_corr) {
				it_read_data->second[(*it_corr).first] = (*it_corr).second;
			}
		}
		//std::cout << it_read_data->second << "\n";	     			     
		corrected_reads.push_back(it_read_data->second);	     
		
	}

	std::cout << "After correction:" << kmers.size() << std::endl;
	
	return corrected_reads;
	
}

RcppExport SEXP BuildRadSites(SEXP fnames, SEXP ks, SEXP merge_sites, SEXP correction, SEXP rad_site_length, SEXP max_distance) 
{
	
	std::map<std::string, int> key_list;
	std::map<std::string,int>::iterator it_key_list;
	
	Rcpp::CharacterVector cx(fnames); 
	std::vector<char*> files;
	for (int i=0; i<cx.size(); i++) 
        {
		files.push_back(cx[i]);  
    	} 
	
	cx = Rcpp::CharacterVector(ks);  
	std::string *keys = new std::string[cx.size()];
    	for (int i=0; i<cx.size(); i++) 
    	{  
      		key_list[std::string(cx[i])] = 0;  
		keys[i] = cx[i];
    	} 

	std::map<std::string, std::map<std::string, int> > RS_Counts;
	std::map<std::string, std::map<std::string, int> >::iterator it_RS_Counts;

	unsigned long read_counter = 0;

	if(as<bool>(correction)) {
		std::vector<std::string> reads;		
		reads = CorrectReads(files,as<int>(rad_site_length));
		for(unsigned long i=0; i < reads.size(); ++i) {
			std::string bc = reads[i].substr(0,6);
			std::string rs = reads[i].substr(6,6);
			     
			std::string rad = reads[i].substr(12,as<int>(rad_site_length));
			     
			it_key_list = key_list.find(bc);
			if( (it_key_list != key_list.end() ) && (rs == "TGCAGG") ) 
			{
				it_RS_Counts = RS_Counts.find(rad);
				if( it_RS_Counts == RS_Counts.end() ) 
				{
					//Создаем новый словарь:
					std::map<std::string, int> RS_Counts_keys;
					RS_Counts[rad] = RS_Counts_keys;
					for(int j=0; j<cx.size(); ++j) 
					{
						RS_Counts[rad][keys[j]] = 0;
					}
				}
					
				RS_Counts[rad][bc] += 1;
					
			}

			if( (read_counter > 0) && (read_counter % 100000 == 0 ) )
			{
        			std::cout << "Record: " << read_counter << std::endl;
			}
			read_counter++;
		}

	} else {
	
		for(int f=0; f < files.size(); ++f) {
			int ii = 0;
        		std::string line;
        		std::ifstream in(files[f]);
        		std::vector<std::string> record_block;
		
		        std::cout << "Processing files: " << files[f] << "\n";
	        
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
                	     
				     std::string bc = record_block[1].substr(0,6);
				     std::string rs = record_block[1].substr(6,6);
			     
				     std::string rad = record_block[1].substr(12,as<int>(rad_site_length));
			     
				     it_key_list = key_list.find(bc);
				     if( (it_key_list != key_list.end() ) && (rs == "TGCAGG") ) 
				     {
					it_RS_Counts = RS_Counts.find(rad);
					if( it_RS_Counts == RS_Counts.end() ) 
					{
						//Создаем новый словарь:
						std::map<std::string, int> RS_Counts_keys;
						RS_Counts[rad] = RS_Counts_keys;
						for(int j=0; j<cx.size(); ++j) 
						{
							RS_Counts[rad][keys[j]] = 0;
						}
					}
					
					RS_Counts[rad][bc] += 1;
					
				     }

				     record_block.clear();

				     if( (read_counter > 0) && (read_counter % 100000 == 0 ) )
				     {
        			    	std::cout << "Record: " << read_counter << std::endl;
				     }
				     read_counter++;
                	     
        	        	}	
			}
		}
	}
	
	/*
	std::vector<std::string> rad_sites;
	//Выводим на печать и в файл:
	for(it_RS_Counts = RS_Counts.begin(); it_RS_Counts != RS_Counts.end(); it_RS_Counts++) {
    		// iterator->first = key
    		// iterator->second = value
    		//std::cout << (*it_RS_Counts).first << std::endl;
		rad_sites.push_back((*it_RS_Counts).first);
		std::map<std::string,int>::iterator it_RS_Counts_keys;
		for(it_RS_Counts_keys = (*it_RS_Counts).second.begin(); it_RS_Counts_keys != (*it_RS_Counts).second.end(); it_RS_Counts_keys++) {
			std::cout << (*it_RS_Counts_keys).second << std::endl;
		}
	}
	*/
	std::cout << "Records processed: " << read_counter << std::endl;
	
	if(as<bool>(merge_sites)) {
		std::vector< std::map<std::string, std::map<std::string, int> >  > rad_sites;
		for(it_RS_Counts = RS_Counts.begin(); it_RS_Counts != RS_Counts.end(); it_RS_Counts++) {
			std::map<std::string, std::map<std::string, int> > site;
			site[it_RS_Counts->first] = it_RS_Counts->second;
		
			std::map<std::string, std::map<std::string, int> >::iterator it_RS_Counts_2;

			for(it_RS_Counts_2 = it_RS_Counts; it_RS_Counts_2 != RS_Counts.end(); it_RS_Counts_2++) {		
				if( Hamming( it_RS_Counts->first, it_RS_Counts_2->first, as<int>(rad_site_length) ) <= as<int>(max_distance) ) {
					site[it_RS_Counts_2->first] = it_RS_Counts_2->second;
				}
			}
			rad_sites.push_back(site);
		}
		return Rcpp::wrap( rad_sites );
	}
	
	return Rcpp::wrap( RS_Counts );
}

