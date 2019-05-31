#pragma once

#include "HuffmanTableDataStructures.h"

void HandleCommandArguments(int, char**, unsigned char*);
void InspectFile();
void InitializeTableMaps();
void Encode();
void Decode();
void PermutateAlphabet();
void CreateDecodeTableMap();
void CreateLengthTableMap();
void CreateReverseDisInfoTableMap();
void CreateDisInfoTableMap();
void CreateNPFTableMap();
void MBR(unsigned char, unsigned char*, unsigned char*);
unsigned char ReverseMBR(const unsigned char*, const unsigned char*);
