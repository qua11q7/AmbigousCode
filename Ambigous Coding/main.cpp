//
//  main.cpp
//  enumnpf
//
//  Copyright © 2017 muhammed oguzhan kulekci. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_helper.hpp>

#include "arithmetic_codec.h"
#include "binary_codec.h"

//#include <string.h>
//#include <sstream>
//#include <ctime>

using namespace std;

using namespace sdsl;


struct codewordList{
	unsigned int codewordLength;
	unsigned char* codeword;
};

struct alphabetSymbolType{
    unsigned int symbolID = 0;
    unsigned long int frequency  = 0;
    unsigned int rank ;
	int codewordID = -1;
};

bit_vector rawBits_EliasW[64];//there can be at most 64 different labels 0..63
wt_huff_int<rrr_vector<63>> wt_EliasW;

double staticHFbitspersymbol;
double staticACbitspersymbol;
double adaptiveHFbitspersymbol;
double adaptiveACbitspersymbol;

size_t EliasW_encode(int_vector<> x, unsigned int* WTsizeInBytes)
{
    int_vector<> labels;
    size_t i,j,bitvectorSize;
    uint64_t one=1;
    
    for(i=0; i<64; i++) rawBits_EliasW[i].resize(0);
    
    size_t seqLength = x.size();
    labels.resize(seqLength);
    
    for(i=0; i<seqLength; i++)
    {
        labels[i]=bits::hi(x[i]); // compute the label
        bitvectorSize = rawBits_EliasW[labels[i]].size();//current size of the corresponding rawbits
        if (0==labels[i])
        {
            rawBits_EliasW[0].resize(bitvectorSize+1); //since label 0 requires special processing
            (x[i]==0) ? rawBits_EliasW[0][bitvectorSize]=0 :  rawBits_EliasW[0][bitvectorSize]=1;
        }
        else
        {
            rawBits_EliasW[labels[i]].resize(bitvectorSize+labels[i]); //compute the payload, which is the least significant labels[i] bits. Append these bits to the rawBits[label]
            for(j=0; j<labels[i]; j++) rawBits_EliasW[labels[i]][bitvectorSize+j] = (x[i]) & (one<<j);//insert from least significant to most significant to be able to retrieve in one 64-bit read operation
        }
    }
    
    construct_im(wt_EliasW, labels);//create the wavelet tree of the labels
    
    j = 0;
    for(i=0; i<64; i++)   j += size_in_bytes(rawBits_EliasW[i]);
    
    *WTsizeInBytes =  (unsigned int) size_in_bytes(wt_EliasW);

    return (size_in_bytes(wt_EliasW)+j)*8;
}


size_t ReadAllBytes(char const* filename, unsigned char** data)
{
    ifstream ifs(filename, std::ios::binary|std::ios::ate);
    size_t pos = ifs.tellg();
    
    *data = new unsigned char[pos + 4];
    
    ifs.seekg(0, ios::beg);
    ifs.read((char*)data[0], pos);
    (*data)[pos]=0;
    (*data)[pos+1]=0;
    (*data)[pos+2]=0;
    (*data)[pos+3]=0;
    return pos;
}


unsigned long psi(unsigned int k, unsigned int d, unsigned int v){//returns number of vectors with 'd' dimensions; each dimension can take values '1..k'; the sum of values in d dimensions is 'v'
    
    if ( (v>(k*d)) || (v<d) )
        return 0;
    if (d==1)   return 1;
    if (v==d)   return 1;
    //if (v==k*d) return 1;
    if (v==d+1) return d;
    
    unsigned int lastalpha;
    (k<v-d+1)?(lastalpha=k):(lastalpha=v-d+1);
    int firstalpha;
    ((int)(v+k-k*d)>1)?(firstalpha=v+k-k*d):(firstalpha=1);
    
    unsigned long int sum=0;
    for(unsigned int i=firstalpha;i<=lastalpha;i++) sum += psi(k,d-1,v-i);
    return sum;
}

void index2vector(unsigned long int index, unsigned int* vec, unsigned int k, unsigned int d,unsigned int sum)
{
	for(unsigned int j=0;j<d-1;j++){
		vec[j]=1;
		unsigned long int z;
		while( (z = psi(k,d-j-1,sum-vec[j])) <= index ){
			index -= z;
			vec[j]++;
		}
		sum-=vec[j];
	}
	vec[d-1]=sum;
}
	
	
unsigned long int vector2index(unsigned int* dim, unsigned int d, unsigned int k)
{
    if (d==1) return 0;
    unsigned int vecsum=0;
    for(unsigned int i=0;i<d;i++) vecsum += dim[i];
    unsigned long int index=0;
    for(unsigned int i=1; i < dim[0]; i++) index += psi(k,d-1,vecsum-i);
    index += vector2index(dim+1,d-1,k);
    return index;
}

void testEnumeration(unsigned int k, unsigned int d){// for all possible v values of range 'k .. k*d', find the ith  d dimensional vector, where each dimension can take values from '1..k', and test if this vector is the ith one
		
	unsigned int* vec = new unsigned int[d];
		
	for(unsigned int sum=d; sum<=k*d; sum++){
		unsigned long int q = psi(k,d,sum);
		cout << "Number of vectors\t" << "(k= " << k  << ", d= " << d  << ", v= " << sum << ")=\t" << q << endl;
		for(unsigned int index=0; index<q; index++){
			index2vector(index,vec,k,d,sum);
			cout << index << '\t';
			for(unsigned int a=0;a<d;a++) cout << vec[a] << '\t' ;
			unsigned long int t = vector2index(vec,d,k);
			cout << t << endl;
			if (t!=index){
				cout << "ERROR\t" << "k: " << k  << "\td: " << d  << "\tv: " << sum << endl;
				return;
			}
		}
	}
}

void createNPFCodeWordList(unsigned int numberofcodewords, struct codewordList* codewordDictionary){// codewordList = <'0',1>,<'1',1>,<'00',2>,<'01',2>,<'10',2>,......,<minimalBinaryRespresentation(numberofcodewords+1),log2(numberofcodewords+1)>
	
	for(unsigned long int i=0; i< numberofcodewords; i++){
        codewordDictionary[i].codewordLength = (unsigned char) log2((double)(i+2));
        codewordDictionary[i].codeword       = new unsigned char[codewordDictionary[i].codewordLength+1];
        unsigned long long one = 1;
        unsigned long long t = one << (codewordDictionary[i].codewordLength-1);
        for(unsigned int j=0;j<codewordDictionary[i].codewordLength;j++){
            codewordDictionary[i].codeword[j] = (t&(i+2))?'1':'0';
            t=t>>1;
        }
        codewordDictionary[i].codeword[codewordDictionary[i].codewordLength]=0;
    }
}

int cmpSymbolID(const void* ptr1, const void* ptr2){
    alphabetSymbolType* a1 = (alphabetSymbolType*) ptr1;
    alphabetSymbolType* a2 = (alphabetSymbolType*) ptr2;
    if (a1->symbolID > a2->symbolID) return 1;
    else if (a1->symbolID < a2->symbolID) return -1;
    else return 0;
}

int cmpFreq(const void* ptr1, const void* ptr2){
    alphabetSymbolType* a1 = (alphabetSymbolType*) ptr1;
    alphabetSymbolType* a2 = (alphabetSymbolType*) ptr2;
    if (a1->frequency < a2->frequency) return 1;
    else if (a1->frequency > a2->frequency) return -1;
    else return 0;
}

unsigned int computeStatistics(unsigned char* data, size_t dataSizeInBytes, unsigned int ngram, struct alphabetSymbolType* symbolList,size_t* totalNPFcodeWordLength, double* entropy){
    unsigned long int one = 1;
    unsigned long int alphabetSize = one << (ngram*8);
    unsigned long int symbol,k=0,symbolcount=0;
	unsigned int      numberofobservedSymbols=0;
    while(k<dataSizeInBytes){
        symbol = 0;
        for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + data[k+j];
        symbolList[symbol].frequency++;
        k+=ngram;
        symbolcount++;
    }
    
    // sort the according to frequencies
    qsort(&symbolList[0],alphabetSize,sizeof(struct alphabetSymbolType),cmpFreq);
    //*****************************************
    
 
    
    *entropy = 0.0;
	*totalNPFcodeWordLength = 0;
    double probabilities[alphabetSize];
    
	for(unsigned int i=0; i<alphabetSize;i++){
		if (symbolList[i].frequency>0){
            symbolList[i].rank = i;
            probabilities[i] = (double) symbolList[i].frequency / (double)(dataSizeInBytes/ngram) ;
            
            if (probabilities[i]< 0.0001){
                probabilities[0] -= (0.0001 - probabilities[i]);
                probabilities[i] = 0.0001;
            }
            *totalNPFcodeWordLength += (symbolList[i].frequency * (unsigned char)log2((double)(i+2)));
			symbolList[i].codewordID = (int) i;
            numberofobservedSymbols++;
            /* EMPIRICAL ENTROPY of the FILE as bits per symbol *************************************/
			(*entropy) += (double)symbolList[i].frequency  * log2((double)symbolcount/(double)symbolList[i].frequency);
            /*****************************************************************************************/
		}else{
			break;//since the list is sorted in descending order, when we see a zero, rest is 0 for sure
		}
	}
    (*entropy) /= symbolcount;
   

    //******************************************
    // sort the according to symbolID
    qsort(&symbolList[0],alphabetSize,sizeof(struct alphabetSymbolType),cmpSymbolID);
    //*****************************************

    /* 0 order static Huffman compression ratio as bits/symbol*/
    Static_Huffman_Code staticHuffman;
    staticHuffman.set_distribution(numberofobservedSymbols,probabilities);
    unsigned char* tmpbuffer = new unsigned char[dataSizeInBytes];
    Binary_Codec staticHFencoder(dataSizeInBytes,tmpbuffer);
    staticHFencoder.start_encoder();
    
    for(unsigned int i = 0; i<dataSizeInBytes; i+=ngram){
        symbol = 0;
        for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + data[i+j];
        staticHFencoder.encode(symbolList[symbol].rank, staticHuffman);
    }
    staticHFbitspersymbol = (double) (staticHFencoder.stop_encoder())*8.0 / (double)(dataSizeInBytes/ngram);
   /************************************************************************************************************/
    
    /* 0 order adaptive Huffman compression ratio as bits/symbol*/
    Adaptive_Huffman_Code adaptiveHuffman(numberofobservedSymbols);
    Binary_Codec adaptiveHFencoder(dataSizeInBytes,tmpbuffer);
    adaptiveHFencoder.start_encoder();
    
    for(unsigned int i = 0; i<dataSizeInBytes; i+=ngram){
        symbol = 0;
        for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + data[i+j];
        adaptiveHFencoder.encode(symbolList[symbol].rank, adaptiveHuffman);
    }
    adaptiveHFbitspersymbol = (double) (adaptiveHFencoder.stop_encoder())*8.0 / (double)(dataSizeInBytes/ngram);
    /************************************************************************************************************/

    /* 0 order static Arithmetic Coding compression ratio as bits/symbol*/
    Static_Data_Model staticDataModel;
    staticDataModel.set_distribution(numberofobservedSymbols,probabilities);
    Arithmetic_Codec staticACencoder(dataSizeInBytes);
    staticACencoder.start_encoder();
    
    for(unsigned int i = 0; i<dataSizeInBytes; i+=ngram){
        symbol = 0;
        for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + data[i+j];
        staticACencoder.encode(symbolList[symbol].rank, staticDataModel);
    }
    staticACbitspersymbol = (double) (staticACencoder.stop_encoder())*8.0 / (double)(dataSizeInBytes/ngram);
    /************************************************************************************************************/
    
    /* 0 order adaptive Arithmetic Coding compression ratio as bits/symbol*/
    Adaptive_Data_Model adaptiveDataModel(numberofobservedSymbols);
    Arithmetic_Codec adaptiveACencoder(dataSizeInBytes);
    adaptiveACencoder.start_encoder();
    
    for(unsigned int i = 0; i<dataSizeInBytes; i+=ngram){
        symbol = 0;
        for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + data[i+j];
        adaptiveACencoder.encode(symbolList[symbol].rank, adaptiveDataModel);
    }
   adaptiveACbitspersymbol = (double) adaptiveACencoder.stop_encoder()*8.0 / (double)(dataSizeInBytes/ngram);
    /************************************************************************************************************/
    return numberofobservedSymbols;
}




int main(int argc, const char * argv[]) {
    
	for(int v=20; v>0;v--){
		for(int d=7; d>0;d--){
			cout << psi(6,d,v) << '\t';
		}
		cout << endl;
	}
		
		
	
    
    //for(unsigned i=6;i<=42;i++) cout << psi(7,6,i) << endl;
    
    //for(unsigned i=4;i<=28;i++) cout << psi(7,4,i) << endl;
    
    //for(unsigned i=2;i<=14;i++) cout << psi(7,2,i) << endl;
    
    //testEnumeration(7,8);
	//for(unsigned int k=7; k<=21;k+=7)
		//for(unsigned int d=4; d<=8;d++)
			//cout << "k = " << k << " d = "<< d << " v = " << k + (k*d-d)/2 << " psi (k,d,v) = " << psi(k,d,(k*d-d)/2) << endl;
	
	// usage filename ngram sparseFactor 
	
    unsigned char* inputData;
    unsigned int   ngram        = atoi(argv[2]);
    unsigned int   sparseFactor = atoi(argv[3]);
    
	size_t inputLen             = ReadAllBytes(argv[1], &inputData);
	size_t alphabetSize         = 1<<(ngram*8);
	
    struct alphabetSymbolType symbols[alphabetSize];
	for(unsigned int i=0;i<alphabetSize;i++) symbols[i].symbolID=i;
    
	size_t totalNPFlen               = 0;
	double entropyTheoretical        = 0.0;
    unsigned int observedSymbolCount = computeStatistics(inputData, inputLen, ngram, &symbols[0],&totalNPFlen,&entropyTheoretical);
	unsigned int maxCodeWordLength   = (unsigned int) log2((double) (observedSymbolCount+1));
    
	struct codewordList NPFcodewordDict[alphabetSize];
	createNPFCodeWordList(observedSymbolCount, &NPFcodewordDict[0]);
	 
    size_t symbol_count=0,codewordBitArrayPosition=0;
    unsigned int blockBitLength=0;
	unsigned int numberOfDistinctBlockLengths = (maxCodeWordLength-1)*sparseFactor+1; 
    unsigned int block_codewordlens[sparseFactor];
        
    // statistics for the block_length and enumID 
    unsigned int statBlockLength[numberOfDistinctBlockLengths];
    for(unsigned int i=0;i<numberOfDistinctBlockLengths;i++) statBlockLength[i] = 0;
	
	unsigned int** statEnumID = new unsigned int*[numberOfDistinctBlockLengths];
    for(unsigned int i=0;i<numberOfDistinctBlockLengths;i++){
        unsigned int numberofitems = psi(maxCodeWordLength,sparseFactor,i+sparseFactor);
        statEnumID[i] = new unsigned int[numberofitems];
        for(unsigned int h=0;h<numberofitems;h++) statEnumID[i][h]=0;
    }
    //*************************************************************
	
    
    //testEnumeration(maxCodeWordLength,sparseFactor);
    
	bit_vector codewordBitArray(totalNPFlen);

	Arithmetic_Codec blockLengthCoder(0xF0000000U);
	Adaptive_Data_Model blockLengthDataModel(sparseFactor*(maxCodeWordLength-1)+1);
	blockLengthCoder.start_encoder();
	
	
		unsigned int middleIndex = ((sparseFactor*(maxCodeWordLength-1))-1)/2;
		
		if (psi(maxCodeWordLength, sparseFactor, middleIndex+sparseFactor) < (1<<11)){
			
			/*adaptive AC codec to encode enumID data*/
			Arithmetic_Codec     enumIDEncoder(0xF0000000U);
			Adaptive_Data_Model* enumIDDataModel[(sparseFactor*(maxCodeWordLength-1))-1];
			unsigned int numberOfModels = (sparseFactor*(maxCodeWordLength-1))-1; 
			for(unsigned int i=0; i<numberOfModels; i++){
				enumIDDataModel[i]  = new Adaptive_Data_Model(psi(maxCodeWordLength, sparseFactor, sparseFactor+i+1));
			}
			enumIDEncoder.start_encoder();		
		
			for(size_t i=0;i<inputLen;i+=ngram){
				
				unsigned long int symbol = 0;
				for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + inputData[i+j];
				for(unsigned int j=0;j< NPFcodewordDict[symbols[symbol].codewordID].codewordLength ;j++){
					if (NPFcodewordDict[symbols[symbol].codewordID].codeword[j]=='0')
						codewordBitArray[codewordBitArrayPosition] = 0;
					else
						codewordBitArray[codewordBitArrayPosition] = 1;
					codewordBitArrayPosition++;
				}
				
				block_codewordlens[symbol_count % sparseFactor] = NPFcodewordDict[symbols[symbol].codewordID].codewordLength;
				symbol_count++;
				blockBitLength += NPFcodewordDict[symbols[symbol].codewordID].codewordLength;
				
				if (0 == symbol_count % sparseFactor){
					
					statBlockLength[blockBitLength-sparseFactor]++;
					blockLengthCoder.encode(blockBitLength - sparseFactor, blockLengthDataModel);
					
					if ((blockBitLength!=sparseFactor) && (blockBitLength != (sparseFactor*maxCodeWordLength)) ){
						unsigned int id = vector2index(block_codewordlens,sparseFactor,maxCodeWordLength);
						//unsigned int modelID = min(sparseFactor*maxCodeWordLength - blockBitLength,blockBitLength - sparseFactor-1);  
						enumIDEncoder.encode(id, *enumIDDataModel[blockBitLength - sparseFactor-1]);
						statEnumID[blockBitLength-sparseFactor][id]++;
					}
					blockBitLength=0;
				}
			}
			
			unsigned long int blockLengthsCompressedSizeInBits = 8 * blockLengthCoder.stop_encoder();
			unsigned long int enumIDCompressedSizeInBits = 8 * enumIDEncoder.stop_encoder();
			/*
			cout << argv[1] << " File size in bytes: " << inputLen << endl;
			cout << "Alphabet size:\t" << observedSymbolCount << endl;
			cout << "Block Length        in bits/symbol " << (double)blockLengthsCompressedSizeInBits/(double)inputLen << endl;
			cout << "CodeWordStream      in bits/symbol " << (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count << endl;
			cout << "EnumID Arithmetic   in bits/symbol " << (double)enumIDCompressedSizeInBits / (double) symbol_count << endl << endl;

			cout << "TOTAL space usage   in bits/symbol " << (double)blockLengthsCompressedSizeInBits/(double)inputLen + (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count + (double)enumIDCompressedSizeInBits / (double) symbol_count << endl;
			cout << "Theoretical entropy in bits/symbol "   << entropyTheoretical << endl;
			cout << "Static Huffman      in bits/symbol "   << staticHFbitspersymbol << endl;
			cout << "Static Arithmetic   in bits/symbol "   << staticACbitspersymbol << endl;
			cout << "Adaptive Huffman    in bits/symbol " << adaptiveHFbitspersymbol << endl;
			cout << "Adaptive Arithmetic in bits/symbol " << adaptiveACbitspersymbol << endl;
			*/
			cout << argv[1] << '\t' << inputLen << '\t';
			cout << observedSymbolCount<< '\t';
			cout << (double)blockLengthsCompressedSizeInBits/(double)inputLen << '\t';
			cout << (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count << "\tA\t";
			cout << (double)enumIDCompressedSizeInBits / (double) symbol_count << '\t';

			cout << (double)blockLengthsCompressedSizeInBits/(double)inputLen + (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count + (double)enumIDCompressedSizeInBits / (double) symbol_count << '\t';
			cout << entropyTheoretical      << '\t' ;
			cout << staticHFbitspersymbol   << '\t' ;
			cout << staticACbitspersymbol   << '\t' ;
			cout << adaptiveHFbitspersymbol << '\t' ;
			cout << adaptiveACbitspersymbol << '\t' << endl;
			
		}else{	
			/*adaptive Huffman codec to encode enumID data*/
			Adaptive_Huffman_Code* adaptiveHFcode[(sparseFactor*(maxCodeWordLength-1))-1];
			for(unsigned int i=0; i<(sparseFactor*(maxCodeWordLength-1))-1; i++){
				adaptiveHFcode[i] = new Adaptive_Huffman_Code(psi(maxCodeWordLength, sparseFactor, sparseFactor+i+1));
			}
			Binary_Codec enumIDEncoder_HF(0xF0000000U);
			enumIDEncoder_HF.start_encoder();

		   for(size_t i=0;i<inputLen;i+=ngram){
				
				unsigned long int symbol = 0;
				for(unsigned int j=0;j<ngram;j++) symbol = symbol*256 + inputData[i+j];
				for(unsigned int j=0;j< NPFcodewordDict[symbols[symbol].codewordID].codewordLength ;j++){
					if (NPFcodewordDict[symbols[symbol].codewordID].codeword[j]=='0')
						codewordBitArray[codewordBitArrayPosition] = 0;
					else
						codewordBitArray[codewordBitArrayPosition] = 1;
					codewordBitArrayPosition++;
				}
				
				block_codewordlens[symbol_count % sparseFactor] = NPFcodewordDict[symbols[symbol].codewordID].codewordLength;
				symbol_count++;
				blockBitLength += NPFcodewordDict[symbols[symbol].codewordID].codewordLength;
				
				if (0 == symbol_count % sparseFactor){
					
					statBlockLength[blockBitLength-sparseFactor]++;
					blockLengthCoder.encode(blockBitLength - sparseFactor, blockLengthDataModel);
					
					if ((blockBitLength!=sparseFactor) && (blockBitLength != (sparseFactor*maxCodeWordLength)) ){
						unsigned int id = vector2index(block_codewordlens,sparseFactor,maxCodeWordLength);
						enumIDEncoder_HF.encode(id, *adaptiveHFcode[blockBitLength - sparseFactor-1]);
                        statEnumID[blockBitLength-sparseFactor][id]++;
					}else if (blockBitLength==sparseFactor) statEnumID[blockBitLength-sparseFactor][0]++;
                    else if (blockBitLength != (sparseFactor*maxCodeWordLength)) statEnumID[blockBitLength-sparseFactor][0]++;
                    
					blockBitLength=0;
				}
			}
			
			unsigned long int blockLengthsCompressedSizeInBits = 8 * blockLengthCoder.stop_encoder();
			unsigned long int HFenumIDCompressedSizeInBits = 8 * enumIDEncoder_HF.stop_encoder();
			/*
			cout << argv[1] << " File size in bytes: " << inputLen << endl;
			cout << "Alphabet size:\t" << observedSymbolCount << endl;
			cout << "Block Length        in bits/symbol " << (double)blockLengthsCompressedSizeInBits/(double)inputLen << endl;
			cout << "CodeWordStream      in bits/symbol " << (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count << endl;
			cout << "EnumID Huffman      in bits/symbol " << (double)HFenumIDCompressedSizeInBits / (double) symbol_count << endl << endl;

			cout << "TOTAL space usage   in bits/symbol " << (double)blockLengthsCompressedSizeInBits/(double)inputLen + (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count + (double)HFenumIDCompressedSizeInBits / (double) symbol_count << endl;
			cout << "Theoretical entropy in bits/symbol "   << entropyTheoretical << endl;
			cout << "Static Huffman      in bits/symbol "   << staticHFbitspersymbol << endl;
			cout << "Static Arithmetic   in bits/symbol "   << staticACbitspersymbol << endl;
			cout << "Adaptive Huffman    in bits/symbol " << adaptiveHFbitspersymbol << endl;
			cout << "Adaptive Arithmetic in bits/symbol " << adaptiveACbitspersymbol << endl;
			*/
			cout << argv[1] << '\t' << inputLen << '\t';
			cout << observedSymbolCount<< '\t';
			cout << (double)blockLengthsCompressedSizeInBits/(double)inputLen << '\t';
			cout << (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count << "\tH\t";
			cout << (double)HFenumIDCompressedSizeInBits / (double) symbol_count << '\t';

			cout << (double)blockLengthsCompressedSizeInBits/(double)inputLen + (double)size_in_bytes(codewordBitArray)*8.0 / (double) symbol_count + (double)HFenumIDCompressedSizeInBits / (double) symbol_count << '\t';
			cout << entropyTheoretical      << '\t' ;
			cout << staticHFbitspersymbol   << '\t' ;
			cout << staticACbitspersymbol   << '\t' ;
			cout << adaptiveHFbitspersymbol << '\t' ;
			cout << adaptiveACbitspersymbol << '\t' << endl;
		}	 
   
    
 
   //print statistics of the block_length and enumID
    for(unsigned int i=0;i<numberOfDistinctBlockLengths;i++) 
		cout << "statBlockLength[" << i+sparseFactor <<"]=" << statBlockLength[i] << endl;
	
    for(unsigned int i=0;i<numberOfDistinctBlockLengths;i++){
        unsigned int numberofitems = psi(maxCodeWordLength,sparseFactor,i+sparseFactor);
        unsigned int qaz=0;
        if ((i+sparseFactor)==15)
            for(unsigned int h=0;h<numberofitems;h++){
                cout << i+sparseFactor << " " << h << "  " << statEnumID[i][h] << endl;
                qaz+=statEnumID[i][h];
            }
        cout << i+sparseFactor << " " << qaz << endl;
    }
    /*************************************************************/
 
  return 0;
}
