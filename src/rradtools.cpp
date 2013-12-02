#include <stdlib.h>     /* malloc, free, rand */
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP BuildRadSites(SEXP fnames, SEXP ks) 
{
	std::map<std::string, int> key_list;
	std::map<std::string,int>::iterator it_key_list;
	
	Rcpp::CharacterVector cx(fnames); 
	std::string *files = new std::string[cx.size()];
	int num_files = cx.size();
	for (int i=0; i<cx.size(); i++) 
        {  
      		files[i] = cx[i];  
    	} 
	
	cx = Rcpp::CharacterVector(ks);  
	std::string *keys = new std::string[cx.size()];
    	for (int i=0; i<cx.size(); i++) 
    	{  
      		key_list[std::string(cx[i])] = 0;  
		keys[i] = cx[i];
    	} 
	
	/*
	key_list["CCTGTG"] = 0; key_list["CTGCGA"] = 0; key_list["GAGGGA"] = 0; key_list["CATAGA"] = 0; key_list["GATCCA"] = 0;
	key_list["CGTTAA"] = 0; key_list["ACTGAT"] = 0; key_list["AACAAC"] = 0; key_list["GTACCG"] = 0; key_list["ATCGGG"] = 0;
	key_list["CATCAG"] = 0; key_list["AGTCAC"] = 0; key_list["GCACAC"] = 0; key_list["TCGTCA"] = 0; key_list["TCCACG"] = 0;
	key_list["TTCGAC"] = 0; key_list["ATAGTT"] = 0; key_list["GGTAAG"] = 0; key_list["GGGATT"] = 0; key_list["ACCTAA"] = 0;
	key_list["ACCAGT"] = 0; key_list["CTAGGC"] = 0; key_list["TCACGG"] = 0; key_list["TTCCCA"] = 0;
	
	std::string keys[] = {"CCTGTG", "CTGCGA", "GAGGGA", "CATAGA", "GATCCA", "CGTTAA", "ACTGAT", "AACAAC", "GTACCG", "ATCGGG", "CATCAG", "AGTCAC", "GCACAC", "TCGTCA", "TCCACG", "TTCGAC", "ATAGTT", "GGTAAG", "GGGATT", "ACCTAA", "ACCAGT", "CTAGGC", "TCACGG", "TTCCCA"};
	*/
	//char *files[] = {"first4k.fastq"};

	std::map<std::string, std::map<std::string, int> > RS_Counts;
	std::map<std::string, std::map<std::string, int> >::iterator it_RS_Counts;

	
	
	unsigned long read_counter = 0;
	
	for(int f=0; f < num_files; ++f) {
		int ii = 0;
        	std::string line;
        	std::fstream in(files[f].c_str());
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
			     
			     std::string rad = record_block[1].substr(12,20);
			     
			     it_key_list = key_list.find(bc);
			     if( (it_key_list != key_list.end() ) && (rs == "TGCAGG") ) 
			     {
				it_RS_Counts = RS_Counts.find(rad);
				if( it_RS_Counts == RS_Counts.end() ) 
				{
					//Создаем новый словарь:
					std::map<std::string, int> RS_Counts_keys;
					RS_Counts[rad] = RS_Counts_keys;
					for(int j=0; j<24; ++j) 
					{
						RS_Counts[rad][keys[j]] = 0;
					}
				}
				
				RS_Counts[rad][bc] += 1;
				
			     }

			     record_block.clear();

			     if( read_counter % 100000 == 0 ) 
			     {
        		    	std::cout << "Record: " << read_counter << std::endl;
			     }
			     read_counter++;
                	     
                	}	
		}
	}
	
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
	
	std::cout << "Records processed: " << read_counter << std::endl;

	//free(files);
	
	return Rcpp::wrap( rad_sites );
}
