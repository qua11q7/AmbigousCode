#include "CodingParameters.h"
#include "HuffmanTableDataStructures.h"
#include "FunctionDeclerations.h"
#include "Utilities.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace sdsl;

unsigned char encode;
unsigned char permutatedAlphabet[SYMBOL_SIZE];
unsigned char reversePermutatedAlphabet[SYMBOL_SIZE];
char* originalFileName;
char* payloadFileName;
char* disInfoFileName;
char* seedFileName;

npfTableMap npfTM;       // 256 kb
disInfoTableMap diTM;    // 9 kb
decodeTableMap decodeTM; // 6 kb
lengthTableMap lengthTM; // 128 byte
reverseDITableMap rdiTM; // 20 kb

// USAGE
// Encode : fileName                -> Generates payload, disInfo, RRR vector and seed
// Decode : payload disInfo seed    -> Generates the original file
int main(int argc, char **argv) {
    HandleCommandArguments(argc, argv, &encode);
	if (encode == -1) {
        return -1;
	}

	InitializeTableMaps();

	if (encode) {
		Encode();

		free(payloadFileName);
		free(disInfoFileName);
        free(seedFileName);
	}
	else {
		Decode();

		free(originalFileName);
	}
}

void HandleCommandArguments(int argc, char** argv, unsigned char* encode) {
	*encode = -1;
	if (argc == 2) {
        *encode = 1;
        originalFileName = argv[1];

        // We want to add "_payload.bin" and "_disInfo.bin" to end of the encoded file names.
        // The length of these strings are 12 charachter long.
        char* fileNameWithoutExtension = remove_ext(originalFileName, '.', '/');
        int payloadFileNameLength = strlen(fileNameWithoutExtension) + 12 + 1;

        payloadFileName = (char*)malloc(payloadFileNameLength + 1);
        snprintf(payloadFileName, payloadFileNameLength, "%s%s", fileNameWithoutExtension, "_payload.bin");
        disInfoFileName = (char*)malloc(payloadFileNameLength + 1);
        snprintf(disInfoFileName, payloadFileNameLength, "%s%s", fileNameWithoutExtension, "_disInfo.bin");
        seedFileName = (char*)malloc(payloadFileNameLength + 1);
        snprintf(seedFileName, payloadFileNameLength - 3, "%s%s", fileNameWithoutExtension, "_seed.csv");

        InspectFile();

        free(fileNameWithoutExtension);
	}
	else if(argc == 4) {
        *encode = 0;
        payloadFileName = argv[1];
        disInfoFileName = argv[2];
        seedFileName = argv[3];

        // We want to add "_decoded.bin" to the end of the decoded file name.
        char* fileNameWithoutExtension = remove_ext(payloadFileName, '.', '/');
        int decodedFileNameLength = strlen(fileNameWithoutExtension) + 12 + 1;

        originalFileName = (char*)malloc(decodedFileNameLength + 1);
        snprintf(originalFileName, decodedFileNameLength, "%s%s", fileNameWithoutExtension, "_decoded.bin");

        free(fileNameWithoutExtension);
	}
}

void InspectFile() {
    FILE *fptr = fopen(originalFileName, "rb");
    fseek(fptr, 0, SEEK_END);
    unsigned long long int fileLength = ftell(fptr);
    fseek(fptr, 0, SEEK_SET);
    unsigned char* data = (unsigned char*)malloc(fileLength * sizeof(char));
    memset(data, 0, fileLength * sizeof(char));
    fread(data, 1, fileLength, fptr);
    fclose(fptr);

    unsigned long long int frequencies[256];
    unsigned char indexes[256];

    memset(frequencies, 0, 256 * sizeof(unsigned long long int));

    for(int i = 0; i < 256; i++) {
        indexes[i] = i;
    }
    for(int i = 0; i < fileLength; i++) {
        int index = data[i];
        frequencies[index] = frequencies[index] + 1;
    }

    // sort array
    for(int i = 0; i < 256; i++) {
        for(int j = i+1; j < 256; j++) {
            if(frequencies[i] < frequencies[j]) {
                unsigned long long int temp = frequencies[i];
                frequencies[i] = frequencies[j];
                frequencies[j] = temp;

                unsigned char tempIndex = indexes[i];
                indexes[i] = indexes[j];
                indexes[j] = tempIndex;
            }
        }
    }

    // calculate entropy
    double entropy = 0;
    for(int i = 0; i < 256; i++) {
        unsigned long long int f = frequencies[i];
        if(f > 0) {
            double p = (double)f / fileLength;
            entropy += p * log2(1/p);
        } else {
            break;
        }
    }
    printf("Entropy %f \r\n", entropy);

    FILE *seedFile = fopen(seedFileName, "w");
    for(int i = 0; i < 256; i++) {
        fprintf(seedFile, "%d;%d\n", indexes[i], i);
    }

    fclose(seedFile);
}

void InitializeTableMaps() {
    PermutateAlphabet();
    if(encode) {
        CreateNPFTableMap();
        CreateDisInfoTableMap();
    } else {
        CreateDecodeTableMap();
        CreateLengthTableMap();
        CreateReverseDisInfoTableMap();
    }
}

void Encode() {
    FILE *fptr = fopen(originalFileName, "rb");
    clock_t start, stop;

    start = clock();
    fseek(fptr, 0, SEEK_END);
    unsigned long long int fileLength = ftell(fptr);
    fseek(fptr, 0, SEEK_SET);
    unsigned char* data = (unsigned char*)malloc((fileLength) * sizeof(char));
    memset(data, 0, (fileLength) * sizeof(char));
    fread(data, 1, fileLength, fptr);
    fclose(fptr);

    stop = clock();

    double elapsed_secs_read = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("\r\nReading file took %f seconds.", elapsed_secs_read);

	clock_t startEncode, stopEncode;
	startEncode = clock();

	// We need to figure out how much space we need to allocate before hand.
	// approximately, disInfo takes 2/d size of the original data
	// payload takes d-2/d size of the original data
	unsigned char* encodedData = (unsigned char*)malloc(fileLength);
	unsigned char* disInfoData = (unsigned char*)malloc(fileLength);

    memset(encodedData, 0, fileLength);
    memset(disInfoData, 0, fileLength);

	long long int currentPos = 0;
	long long int currentIndex = 0;
	long long int disInfoIndex = 4;
	unsigned char disInfoTableIndex = 1;
	unsigned char tableIndex = 1;

	unsigned char* valuePtr = &data[0];
	huffmanRow* encodeRow;
	huffmanRow* disInfoRow;

	while (currentPos < fileLength) {
		valuePtr = &data[currentPos];
		++currentPos;

		// We need to split the data according to D_SIZE.
		encodeRow = &(npfTM.tables[tableIndex].rows[*valuePtr]);
		encodedData[currentIndex] = encodeRow->output;
		currentIndex += encodeRow->canOutput;
		tableIndex = encodeRow->nextTableId;


		disInfoRow = &(diTM.tables[disInfoTableIndex].rows[encodeRow->mbrLength]);
		disInfoData[disInfoIndex] = disInfoRow->output;
		disInfoIndex += disInfoRow->canOutput;
		disInfoTableIndex = disInfoRow->nextTableId;
	}

	unsigned char skipAmount = 0;
	// 0    -> Don't skip any DisInfo length while decoding.
	// 8    -> Skip last DisInfo length while decoding.
	// 16   -> Skip last two DisInfo lengths while decoding.
	// 24   -> Skip last three DisInfo lengths while decoding.

	if(disInfoTableIndex != 1) {
        unsigned char lastDIOutput, lastDIRemaining;
        MBR(disInfoTableIndex, &lastDIOutput, &lastDIRemaining);
        disInfoData[disInfoIndex] = lastDIOutput << (8 - lastDIRemaining);
        ++disInfoIndex;

        if(lastDIRemaining == 1)
            skipAmount = 8;
	}
	if(tableIndex != 1) {
        unsigned char lastPLOutput, lastPLRemaining;
        MBR(tableIndex, &lastPLOutput, &lastPLRemaining);
        encodedData[currentIndex] = lastPLOutput << (8 - lastPLRemaining);
        ++currentIndex;
	}

	// Create Header
	unsigned char dType = 32;
	// 0    -> D4
	// 32   -> D8
	// 64   -> D12
	// 96   -> D16
	// 128  -> D20

	unsigned char remainingDataCount = 0;
	// 0    -> No remaining data left
	// 2    -> 1 byte of remaining data left
	// 4    -> 2 bytes of remaining data left
	// 6    -> 3 bytes of remaining data left

	disInfoData[0] = dType | skipAmount | remainingDataCount | 1;
	disInfoData[1] = 0; // First byte of the remaining data
	disInfoData[2] = 0; // Second byte of the remaining data
	disInfoData[3] = 0; // Third byte of the remaining data

	bit_vector b(disInfoIndex * 8, 0);
	for(unsigned long long int i = 0; i < disInfoIndex; i++) {
        unsigned char value = disInfoData[i];
        for(int j = 0; j < 8; j++) {
            unsigned char lsb = (value >> j) & 1;
            if(lsb == 1) {
                b[i * 8 + j] = 1;
            }
        }
	}

	rrr_vector<63> rrrb(b);

	printf("\r\nBit vector size in MB: %f", size_in_mega_bytes(b));
	printf("\r\nRRR vector size in MB: %f", size_in_mega_bytes(rrrb));

	std::ofstream outRRR("rrr_vector_data.sdsl");
	rrrb.serialize(outRRR);

	FILE* fDisInfo = fopen(disInfoFileName, "wb");
	fwrite(disInfoData, 1, disInfoIndex, fDisInfo);
	fclose(fDisInfo);

	FILE* fEncoded = fopen(payloadFileName, "wb");
	fwrite(encodedData, 1, currentIndex, fEncoded);
	fclose(fEncoded);

	stopEncode = clock();
	double elapsed_secs = (double)(stopEncode - startEncode) / CLOCKS_PER_SEC;
	printf("\r\nEncoding time: %f seconds.", elapsed_secs);
	printf("\r\nEncoding speed: %f MiB/sec\r\n", (fileLength / 1024.0f / 1024.0f) / elapsed_secs);

    long long int disInfoSize = disInfoIndex;
    long long int payloadSize = currentIndex;
    long long int totalSize = disInfoSize + payloadSize;
    long long int overheadSize = totalSize - fileLength;
    float overheadPercentage = (overheadSize * 100.0f) / totalSize;

    printf("\r\nOriginal Size: %lld", fileLength);
    printf("\r\nPayload Size: %lld", payloadSize);
    printf("\r\nPayload Percentage: %f", payloadSize * 100.0f / fileLength);
    printf("\r\nDisInfo Size: %lld", disInfoSize);
    printf("\r\nDisInfo Percentage: %f", disInfoSize * 100.0f / fileLength);
    printf("\r\nTotal Size: %lld", totalSize);
    printf("\r\nOverhead Size: %lld", overheadSize);
    printf("\r\nOverhead Percentage: %f\r\n\r\n", overheadPercentage);

    free(data);
	free(encodedData);
	free(disInfoData);
}

void Decode() {
	FILE *fPayload = fopen(payloadFileName, "rb");
	FILE *fDisInfo = fopen(disInfoFileName, "rb");
	clock_t start, stop;

	start = clock();
	fseek(fDisInfo, 0, SEEK_END);
	unsigned long long int fileLenDisInfo = ftell(fDisInfo);
	fseek(fDisInfo, 0, SEEK_SET);
	unsigned char* disInfoData = (unsigned char*)malloc((fileLenDisInfo) * sizeof(char));
	fread(disInfoData, 1, fileLenDisInfo, fDisInfo);
	fclose(fDisInfo);

	fseek(fPayload, 0, SEEK_END);
	unsigned long long int fileLenPayload = ftell(fPayload);
	fseek(fPayload, 0, SEEK_SET);
	unsigned char* payloadData = (unsigned char*)malloc((fileLenPayload) * sizeof(char));
	fread(payloadData, 1, fileLenPayload, fPayload);
	fclose(fPayload);
	stop = clock();

	double elapsed_secs = (double)(stop - start) / CLOCKS_PER_SEC;
	printf("Reading payload and disInfo took %f seconds.", elapsed_secs);

	clock_t decodeStart, decodeStop;
	decodeStart = clock();

	// we know for sure that decodedData size must be smaller than the total size of DisInfo and Payload.
	unsigned char* decodedData = (unsigned char*)malloc((fileLenDisInfo + fileLenPayload) * 5);
	long long int decodedDataIndex = 0;

	long long int currentPayloadIndex = 1;
	unsigned char currentPayloadByte = payloadData[currentPayloadIndex];
	unsigned char prevPayloadByte = payloadData[currentPayloadIndex - 1];
	unsigned short payloadShort = prevPayloadByte << 8 | currentPayloadByte;
	unsigned char currentPayloadShortBitIndex = 0;

	unsigned char decodeTableIndex = payloadData[0];
	unsigned char disInfoTableIndex = 0;
	unsigned char lengthTableIndex = 0;

	unsigned char headerData = disInfoData[0];
	unsigned char dTypeData = headerData >> 5;      // Get first 3 most significant bits
	unsigned char skipData = (headerData >> 3) & 3; // Get 4th and 5th most significant bits
	unsigned char remainingCountData = (headerData >> 1) & 3; // Get 2nd and 3rd least significant bits
	unsigned char isSeededData = headerData & 1;    // Get the least significant bit.

	long long int currentPos = 4;   // Skip Header
	while (currentPos < fileLenDisInfo) {
		unsigned char data = disInfoData[currentPos];
		++currentPos;

		reverseDIHRow* row = &(rdiTM.tables[disInfoTableIndex].rows[data]);
		for (unsigned char i = 0; i < row->outputCount; i++)
		{
            if(currentPos == fileLenDisInfo && i == (row->outputCount - skipData)) {
                break;
            }

			int length = row->outputs[i];

			lengthRow* lengthRow = &(lengthTM.tables[lengthTableIndex].rows[length - 1]);
			lengthTableIndex = lengthRow->remainingLength;
			currentPayloadIndex += lengthRow->getNextByte;
			prevPayloadByte = currentPayloadByte;
			currentPayloadByte = payloadData[currentPayloadIndex];
			payloadShort = prevPayloadByte << 8 | currentPayloadByte;

			decodeNPFRow* currentDecodeRow = &(decodeTM.tables[decodeTableIndex].rows[length - 1]);
			decodedData[decodedDataIndex] = currentDecodeRow->output;

			++decodedDataIndex;

			unsigned char shiftedData = (payloadShort << currentPayloadShortBitIndex) >> (8 + currentDecodeRow->remainingLength);

			decodeTableIndex = currentDecodeRow->remaining | shiftedData;

			currentPayloadShortBitIndex = lengthRow->remainingLength;
		}

		disInfoTableIndex = row->nextTableId;
	}

	FILE* fDecoded = fopen(originalFileName, "wb");
	fwrite(decodedData, 1, decodedDataIndex, fDecoded);
	fclose(fDecoded);

	decodeStop = clock();
	elapsed_secs = (double)(decodeStop - decodeStart) / CLOCKS_PER_SEC;
	printf("\r\nDecoding Time: %f seconds.", elapsed_secs);
	printf("\r\nDecoding Speed: %f MiB/sec\r\n\r\n", (decodedDataIndex) / 1024.0f / 1024.0f / elapsed_secs);

	free(decodedData);
	free(payloadData);
	free(disInfoData);
}

void PermutateAlphabet() {
	memset(permutatedAlphabet, 0, SYMBOL_SIZE * sizeof(char));
	memset(reversePermutatedAlphabet, 0, SYMBOL_SIZE * sizeof(char));

    FILE* seedFile = fopen(seedFileName, "r");
    char line[SYMBOL_SIZE];
    while(fgets(line, SYMBOL_SIZE, seedFile)) {
        char* temp1 = strdup(line);
        char* temp2 = strdup(line);

        const char* value1Str = get_field(temp1, 1);
        const char* value2Str = get_field(temp2, 2);
        int value1 = atoi(value1Str);
        int value2 = atoi(value2Str);

        free(temp1);
        free(temp2);

        permutatedAlphabet[value1] = value2;
        reversePermutatedAlphabet[value2] = value1;
    }
}

void CreateNPFTableMap() {
	memset(&npfTM, 0, sizeof(npfTableMap));

	for (unsigned char t = 1; t <= 255 && t > 0; t++) {
		npfTable table;
		memset(&table, 0, sizeof(npfTable));

		unsigned char initialValue = 0;
		unsigned char initialLength = 0;
		if (t != 1) {
			MBR(t, &initialValue, &initialLength);
		}

		// for each table, create rows. SYMBOL_SIZE - 1 many rows.
		for (unsigned char r = 0; r <= (SYMBOL_SIZE - 1); r++)
		{
			huffmanRow row;
			memset(&row, 0, sizeof(huffmanRow));

            unsigned char totalLength = 0;
            unsigned char remainingLength = 0;
            unsigned char remaining = 0;

			unsigned char value = permutatedAlphabet[r];
			if (value < (SYMBOL_SIZE - 2)) {
				value += 2;
			} else {
                // Since D_SIZE is 8, we will always output!
				totalLength = initialLength + D_SIZE;
				row.canOutput = 1;

				if (value == (SYMBOL_SIZE - 2)) {	            // Consider this row as always D_SIZE many zeroes : 00000000
                    // If we can output, it means we have more than or equal to 8 bits stored.
                    // we have to shift output to have 8 bits, put remaining bits to remaining variable.
                    row.output = initialValue << (8 - initialLength);
                    remainingLength = initialLength;
                    remaining = 0;                           // Remaining will always be zero since this row always means D_SIZE many zeroes.
				}
				else if (value == (SYMBOL_SIZE - 1)) {			// Consider this row as always (D_SIZE - 1) many zeroes + 1 : 00000001
                    remainingLength = initialLength;
                    row.output = (initialValue << (8 - initialLength)) | (1 >> (remainingLength));
                    remaining = remainingLength == 0 ? 0 : 1;
				}
				row.nextTableId = ReverseMBR(&remaining, &remainingLength);
                row.mbrLength = D_SIZE;

				table.rows[r] = row;

				if(r == (SYMBOL_SIZE - 1))
                    break;
                continue;
			}

            unsigned char mbrValue = 0;
            unsigned char mbrLength = 0;
			MBR(value, &mbrValue, &mbrLength);

			row.mbrLength = mbrLength;
			totalLength = initialLength + mbrLength;
			row.canOutput = totalLength >= 8 ? 1 : 0;

			if (row.canOutput) {
				// initialValue comes from the table's own starting value.
				// If we add another MBR value to it and the length surpasses D_SIZE (8 bits),
				// we have to output the first 8 bits.
				// Firstly, we shift right the initialValue to make space for the new MBR value.
				// This shifting amount will be determined by how many bits we can shift to left without
				// losing any bit of information. If the initialValue's length (initialLength) is 5 bits,
				// we have to shift it 8 - 5 = 3 bits.
				// After making space for the upcoming MBR bits, we have to calculate how many of the
				// MBR bits we can concat to the end of the shiftedInitialValue. Lets say that MBR length
				// is 5 bits. We have only 3 bits of space to concat. So we have to shift 5 - 3 = 2 times right to
				// get rid of unnecessary bits (bits that won't be written quite yet).
				// Finally we apply or bitwise operation to shiftedInitialValue and shiftedMBRValue to get the output.
				// The right shifted bits of MBR is the remaining bits.

				remainingLength = mbrLength - (8 - initialLength);
				unsigned char shiftedInitialValue = initialValue << (8 - initialLength);
				unsigned char shiftedMBRValue = mbrValue >> remainingLength;
				row.output = shiftedInitialValue | shiftedMBRValue;

				// same logic in MBR is applied here. We need to get the last remaningLength bits of MBRValue.
				// So we are aplying a bitMask to remove first 8 - remaningLength of bits.
				unsigned char bitMask = (1 << remainingLength) - 1;
				remaining = mbrValue & bitMask;

				row.nextTableId = ReverseMBR(&remaining, &remainingLength);
			}
			else {
				row.output = 0;
				unsigned char shiftedInitialValue = initialValue << mbrLength;
				remaining = shiftedInitialValue | mbrValue;
				remainingLength = totalLength;	// the total length is the actual length since we didn't output anything.

				row.nextTableId = ReverseMBR(&remaining, &remainingLength);
			}

			table.rows[r] = row;

			if (r == (SYMBOL_SIZE - 1))
				break;
		}

		npfTM.tables[t] = table;
	}
}

void CreateDisInfoTableMap() {
	memset(&diTM, 0, sizeof(disInfoTableMap));

	for (unsigned char t = 1; t <= 255 && t > 0; t++)
	{
		disInfoTable table;
		memset(&table, 0, sizeof(disInfoTable));

		unsigned char initialValue = 0;
		unsigned char initialLength = 0;
		if (t != 1) {
			MBR(t, &initialValue, &initialLength);
		}

		for (unsigned char r = 1; r <= D_SIZE; r++)
		{
			unsigned char value = 1;
			unsigned char length = r;

			huffmanRow row;
			memset(&row, 0, sizeof(huffmanRow));

            unsigned char totalLength = 0;
            unsigned char remainingLength = 0;
            unsigned char remaining = 0;

			totalLength = initialLength + length;
			row.canOutput = totalLength >= 8 ? 1 : 0;

			if (row.canOutput) {
                remaining = value;
				remainingLength = totalLength - 8;
				row.output = initialValue << (8 - initialLength);
				if(remainingLength == 0) {
                    row.output |= value;
                    remaining = 0;
				}
			}
			else
			{
				row.output = 0;
				unsigned char shiftedInitialValue = initialValue << length;
				remaining = shiftedInitialValue | value;
				remainingLength = totalLength;	// the total length is the actual length since we didn't output anything.
			}

			row.nextTableId = ReverseMBR(&remaining, &remainingLength);

			table.rows[r] = row;
		}

		diTM.tables[t] = table;
	}
}

void CreateDecodeTableMap() {
	memset(&decodeTM, 0, sizeof(decodeTableMap));

	for (unsigned char t = 0; t <= 255; t++) {
		decodeNPFTable decodeTable;
		memset(&decodeTable, 0, sizeof(decodeNPFTable));

		for (unsigned char r = 0; r < D_SIZE; r++) {
			decodeNPFRow row;
			memset(&row, 0, sizeof(decodeNPFRow));

			unsigned char length = r + 1;
			unsigned char codeword = t >> (8 - length);
			if (length != D_SIZE) {
				unsigned char codewordValue = ReverseMBR(&codeword, &length);
				row.output = reversePermutatedAlphabet[codewordValue - 2];
			}
			else {
				if (codeword == 0)	        // if length is D_SIZE, this can only be (D_SIZE-1)0's => SYMBOL_SIZE - 2
					row.output = reversePermutatedAlphabet[SYMBOL_SIZE - 2];
				else if (codeword == 1)    // if length is D_SIZE, this can only be (D_SIZE-1)0's1 => SYMBOL_SIZE - 1
					row.output = reversePermutatedAlphabet[SYMBOL_SIZE - 1];
			}
			row.remaining = t << (length);
			row.remainingLength = 8 - length;

			decodeTable.rows[r] = row;
		}

		decodeTM.tables[t] = decodeTable;
		if (t == 255)
			break;
	}
}

void CreateLengthTableMap() {
	memset(&lengthTM, 0, sizeof(lengthTableMap));

	for (unsigned char t = 0; t < D_SIZE; t++)
	{
		lengthTable table;
		memset(&table, 0, sizeof(lengthTable));

		for (unsigned char r = 0; r < D_SIZE; r++)
		{
			lengthRow row;
            memset(&row, 0, sizeof(lengthRow));

			unsigned char totalLength = t + (r + 1);
			if (totalLength >= 8) {
				row.getNextByte = 1;
				row.remainingLength = totalLength - 8;
			}
			else {
				row.getNextByte = 0;
				row.remainingLength = totalLength;
			}

			table.rows[r] = row;
		}

		lengthTM.tables[t] = table;
	}
}

void CreateReverseDisInfoTableMap() {
    memset(&rdiTM, 0, sizeof(reverseDITableMap));

	// Remaining value can only be zero, only its length can change. If there is a 1 bit in the remaining,
	// we can already calculate the corresponding symbol length. So the remaining value cannot contain any 1 bits.
	// Actually, t can be at most 6, when t is 7, it means there are 7 zeroes. We know that 7 zeroes mean codeword length of 8 and cannot be a remaining data.
	for (unsigned char t = 0; t < D_SIZE -1; t++) {
		reverseDITable table;
		memset(&table, 0, sizeof(reverseDITable));

		unsigned char initialLength = t;

		for (unsigned char r = 1; r <= 255; r++) {
			reverseDIHRow row;
			memset(&row, 0, sizeof(reverseDIHRow));

			row.outputCount = 0;
			unsigned char totalLength = initialLength + 8;
			unsigned short currentValue = r;

			// left align the sequence.
			unsigned short shiftedValue = currentValue << (16 - totalLength);

			unsigned short leftMostBitMask = 32768;		// 10000000 00000000 => if shiftedValue & 32768 == 0, left most bit 0; else == 32768, left most bit 1.
			unsigned char zeroCount = 0;
			for (unsigned char i = 0; i < totalLength; i++) {
				unsigned char leftMostBit = (shiftedValue & leftMostBitMask) == 0 ? 0 : 1;
				shiftedValue <<= 1;
				if (leftMostBit == 0) {
					zeroCount++;
				}
				if(leftMostBit != 0 || zeroCount == 7) {
					row.outputs[row.outputCount] = zeroCount + 1;
					row.outputCount++;
					zeroCount = 0;
				}
			}

			row.nextTableId = zeroCount;

			table.rows[r] = row;
			if (r == 255) {
				break;
			}
		}

		rdiTM.tables[t] = table;
	}
}

void MBR(unsigned char value, unsigned char* mbrValue, unsigned char* mbrLength) {
	// MBR => Minimum Binary Representation, removes the most significant bit and returns the remaining value with its length.
	// Cannot calculate MBR of 0 and 1.
	if (value == 0 || value == 1) {
        unsigned char zero = 0;
        *mbrValue = zero;
        *mbrLength = zero;
        return;
	}

	unsigned char orgValue = value;		// cache the original value
	unsigned char mbrIndex = 0;			// variable to store where the most significant bit (MSB) is
	while (value != 1) {				// shift the value to the right till it becomes 1, the shift amount will give us the MSB.
		value >>= 1;
		++mbrIndex;
	}

	unsigned char bitMask = (1 << mbrIndex) - 1;		// To get rid of the MSB, we have to and the value with 2^(MSBindex) - 1
	unsigned char mbr = orgValue & bitMask;				// Obtain MBR.

	// Ex. value = 23 : 00010111
	// MSB index is 4
	// bitMask : 2^4 - 1 = 15 : 00001111
	// 00010111 & 00001111 = 00000111 (MBR)
	// MBR length is actually 4 bits (0111)
	// so we have to return the length of the MBR as well.
	// MBR length is actually the MBR index.

    *mbrValue = mbr;
    *mbrLength = mbrIndex;
}

unsigned char ReverseMBR(const unsigned char* value, const unsigned char* length) {
	// Ex. value = 7, length = 5 : 00111
	// Result should be 100111 : 39
	// (1 << length) : 100000 | 00111 = 100111
	return (1 << (*length)) | (*value);
}
