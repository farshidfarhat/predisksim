/*
 * predisksim.h
 *
 *  Created on: Sep 19, 2014
 *      Author: user
 */

#ifndef PREDISKSIM_H_
#define PREDISKSIM_H_

#define ull unsigned long long
#define ldbl long double

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;
struct stat st = {0};

//#define folderpath "new_trace_files"
#define WARMUPNUM 0
#define TraceIndex 0
#define MaxTraceNum 2000000
#define MaxTraceFile 2
#define PageSize 8
#define UtilFactor 1.0 // Utilization Factor may decrease/increase the App's IOPS proportionally
#define EstimateFactor 1.0 // Estimate Factor may decrease/increase the HDD & SSD 's IOPS
#define AT_SqueezeFactor 5.0 // Squeeze Factor may squeeze the traces' arrival times
#define BN_SqueezeFactor 1.0 // Squeeze Factor may squeeze the traces' block numbers
#define MAX_SPACE 10000000 // Maximum Address Space of a Trace

#define TraceFileNum 2
string tracefilepath[TraceFileNum] =
{
	//tracefilepath[0] =  "/home/user/workspace/disksim-4.0/frashtest/fin1";//argv[1];
	//tracefilepath[1] =  "/home/user/workspace/disksim-4.0/frashtest/fin2";//argv[1];
		"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin1",// IOPS = 122
		"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin2"// IOPS = 90
		//"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin1",// IOPS = 122
		//"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin2"// IOPS = 90
	//tracefilepath[0] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/hm0" ;// IOPS =
	//tracefilepath[1] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/hm1" ;// IOPS =
	//tracefilepath[1] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/proj1";// IOPS = 40
	//tracefilepath[0] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/proj2";// IOPS = 48
};
//ldbl SpaceRatio[2] = {0.1, 0.9};
//ldbl SpaceRatio[2] = {0.2, 0.8};
//ldbl SpaceRatio[2] = {0.3, 0.7};
//ldbl SpaceRatio[2] = {0.4, 0.6};
ldbl SpaceRatio[2] = {0.5, 0.5};
//ldbl SpaceRatio[2] = {0.6, 0.4};
//ldbl SpaceRatio[2] = {0.7, 0.3};
//ldbl SpaceRatio[2] = {0.8, 0.2};
//ldbl SpaceRatio[2] = {0.9, 0.1};



struct TRACE
{
	ull id; //identifier
	ldbl at; //arrival time
	int dn; //device number
	int fd; //file descriptor
	int md; //read/write
	ull lba; //logical block address
	ull nbs; //number of blocks
	ldbl dt; //departure time
};
TRACE traces [MaxTraceFile][MaxTraceNum+1]; // +1 to have the stats for traces


struct TraceStat
{
	// Primary Requirements
	ull space; //required size == MaxLBA
	ldbl iops; //required IOPS
	ldbl iops_t; //overall IOPS
	ldbl bw; //overall bw
	ldbl f_time; //final arrival time
	ldbl spaceratio; // static ratio assigned to each trace like 0.1 (the other 0.9) or 0.5/0.5
	unsigned int SpaceDist [MAX_SPACE];

	// Secondary Requirements
	ldbl rt; //required response time
	ldbl th; //required throughput
	ldbl ut; //required utilization
	ull HDDspace; // size in HDD
	ull SSDspace; // size in SSD
	ldbl HDDiops; // IOPS in HDD
	ldbl SSDiops; // IOPS in SSD
	ull TraceSize; // number of traces
	ull SpaceFreqThreshold; // Space-Frequency Threshold
	ull MaxNBS;

	TRACE* traces;
};
TraceStat tracestats[MaxTraceFile];

struct HDD
{
	ull space;
	ldbl bw;
	ldbl iops;
};
HDD hdd;

struct SSD
{
	ull space;
	ldbl bw;
	ldbl iops;
};
SSD ssd;

string line;
//int TraceFileNum;
ull MaxTraceSize;
ldbl MaxTime;


ldbl GetNextNumber(string* Line, int* iL);
ldbl GetMaxThroughput(string tracefile, int iTF, int warmupnum);
void GetTraceStat(string tracefile, int iTF, int warmupnum, ldbl spaceratio);
void SetTrace(string tracefile, string newtracefile);
void ArrangeTraces(string hddtracefilepath, string ssdtracefilepath, int nTF, int Algorithm); // Space Partitioning Algorithm: 0 = Basic partitioning, 1 = Smart
void BAAv1(string folderpath, int nTF, int Algorithm);
void RAMTSv1(string folderpath, int nTF, int Algorithm);
void DRFv1(string folderpath, int nTF, int Algorithm);
void DRFv2(string folderpath, int nTF, int Algorithm);
void BAAv2(string folderpath, int nTF, int Algorithm);
void BAAeq(string folderpath, int nTF, int Algorithm);
void RMSv2(string folderpath, int nTF, int Algorithm);
void EQAL_SPAREv1(string folderpath, int nTF, int Algorithm);

#endif /* PREDISKSIM_H_ */
