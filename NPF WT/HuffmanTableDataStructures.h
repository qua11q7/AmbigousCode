#pragma once

#include "CodingParameters.h"

// NPF Coding Structs
// Each row in a npf huffman table is represented by this data structure.
// Holds most of the important data.
typedef struct {
	unsigned char mbrLength;
	unsigned char canOutput;
	unsigned char output;
	unsigned char nextTableId;
} huffmanRow;

// Each table has its own id and rows.
typedef struct {
	huffmanRow rows[SYMBOL_SIZE];
} npfTable;

// Holds all npf huffman tables for ease of access.
typedef struct {
	npfTable tables[256];
} npfTableMap;

// DisInfo Coding Structs
// Huffman Table structure for storing DisInfo rows.
// Created a new data structure because DisInfo uses less table and rows.
typedef struct {
	huffmanRow rows[D_SIZE + 1];
} disInfoTable;

// Holds all DisInfo Huffman Table structure for ease of access.
typedef struct {
	disInfoTable tables[256];
} disInfoTableMap;

// NPF Decoding Structs
typedef struct {
	unsigned char output;
	unsigned char remaining;
	unsigned char remainingLength;
} decodeNPFRow;

typedef struct {
	decodeNPFRow rows[D_SIZE];
} decodeNPFTable;

typedef struct {
	decodeNPFTable tables[256];
} decodeTableMap;

typedef struct {
	unsigned char getNextByte;
	unsigned char remainingLength;
} lengthRow;

typedef struct {
	lengthRow rows[D_SIZE];
} lengthTable;

typedef struct {
	lengthTable tables[D_SIZE];
} lengthTableMap;

// DisInfo Decoding Structs
// This data structure is used in decoding the encoded DisInfo.
// This row can have more than one output and doesn't use input length.
typedef struct {	// Reverse DisInfo Huffman Row
	unsigned char outputs[8];	// a byte can contain at most 8 outputs.
	unsigned char outputCount;
	unsigned char nextTableId;
} reverseDIHRow;

typedef struct {	// Reverse DisInfo Table
	reverseDIHRow rows[256];
} reverseDITable;

typedef struct {		// Reverse DisInfo Table Map
	reverseDITable tables[D_SIZE];
} reverseDITableMap;
