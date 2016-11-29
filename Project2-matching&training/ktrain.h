#ifndef KTRAIN_H
#define KTRAIN_H

#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>

using namespace std;

#define SEGNUM (5)
#define SAMPLENUM (2)
#define INF (100000)

void readSampleDCT(const char* filename, int flag);
void initSeg();
void testPrintSample();
void calSegmentMeansCovar();
void calSampDistance (int num);
void kmeansDistance();
void kBackTrack(int num);
void segKMeansProcess();
#endif