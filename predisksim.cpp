/*
 * predisksim.cpp
 *
 *  Created on: Sep 19, 2014
 *      Author: user
 */

#include "predisksim.h"

int main ( int arc, char **argv )
{
	//string TraceFilePath, NewTraceFilePath;
	//TraceFilePath    =  argv[1]; //"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/CAMRESHMSA01-lvm1.csv";//argv[1];
	//NewTraceFilePath =  argv[2]; //"/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/hm1";//argv[1];
	//SetTrace(TraceFilePath, NewTraceFilePath); // warmup can be 1000, 1000000, ...
	//exit(0);
	// Init the space distribution of all traces
	for(int itf = 0; itf < MaxTraceFile; itf++)
		for(int sc = 0; sc < MAX_SPACE; sc++)
			tracestats[itf].SpaceDist[sc] = 0;
	// Make a dir for new trace files
	string folderpath = "new_trace_files";
	if(stat(&folderpath[0], &st) == -1)
		mkdir(&folderpath[0], 0700);
	///////////////////////// Test an output (.dat file) to get the tracestat /////////////////////
	ldbl HDDMaxThr = GetMaxThroughput("/home/user/workspace/disksim-4.0/frashtest/hdd.dat", 0, WARMUPNUM);
	ldbl SSDMaxThr = GetMaxThroughput("/home/user/workspace/disksim-4.0/frashtest/ssd.dat", 1, WARMUPNUM);
	///////////////////////// Get stats of I/O traces of all apps /////////////////////////////////
	//TraceFileNum = 2;
	//string tracefilepath[TraceFileNum];
	//tracefilepath[0] =  "/home/user/workspace/disksim-4.0/frashtest/fin1";//argv[1];
	//tracefilepath[1] =  "/home/user/workspace/disksim-4.0/frashtest/fin2";//argv[1];
	//tracefilepath[0] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin1";// IOPS = 122
	//tracefilepath[1] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/fin2";// IOPS = 90
	//tracefilepath[0] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/hm0" ;// IOPS =
	//tracefilepath[1] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/hm1" ;// IOPS =
	//tracefilepath[1] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/proj1";// IOPS = 40
	//tracefilepath[0] =  "/media/user/22584d95-e836-4c95-83ef-3cb279fa3d10/traces/proj2";// IOPS = 48
	//ldbl SpaceRatio[2] = {0.1, 0.9};
	//ldbl SpaceRatio[2] = {0.2, 0.8};
	//ldbl SpaceRatio[2] = {0.3, 0.7};
	//ldbl SpaceRatio[2] = {0.4, 0.6};
	//ldbl SpaceRatio[2] = {0.5, 0.5};
	//ldbl SpaceRatio[2] = {0.6, 0.4};
	//ldbl SpaceRatio[2] = {0.7, 0.3};
	//ldbl SpaceRatio[2] = {0.8, 0.2};
	//ldbl SpaceRatio[2] = {0.9, 0.1};

	MaxTraceSize = 1000000000; // set to min of all app
	MaxTime = 1000000000;
	for(int i = 0; i < TraceFileNum; i++) {
		GetTraceStat(tracefilepath[i], i, TraceIndex, SpaceRatio[i]); // warmup can be 1000, 1000000, ...
		if(MaxTime > tracestats[i].f_time)
			MaxTime = tracestats[i].f_time;
		}
	ldbl MaxDiskSize = 0;
	for(int i = 0; i < TraceFileNum; i++)
		if(MaxDiskSize < tracestats[i].space)
			MaxDiskSize = tracestats[i].space;

	// SPACE and IOPS Allocation
	// Set ssd and hdd parameters
	hdd.space = MaxDiskSize * TraceFileNum; // must be set in disksim config
	ssd.space = hdd.space / TraceFileNum; // must be set in disksim config
	hdd.iops = EstimateFactor * HDDMaxThr; // or must be set in disksim config by response time
	ssd.iops = EstimateFactor * SSDMaxThr;
	cout << "HDD Size = " << hdd.space << " should be at least in hdd.parv\n";
	cout << "SSD Size = " << ssd.space << " should be at least in ssd.parv\n";
	cout << "HDD IOPS = " << hdd.iops << "\n";
	cout << "SSD IOPS = " << ssd.iops << "\n";

	// Generate the new traces
	BAAeq(folderpath, 2, 0);
	DRFv2(folderpath, 2, 0);
	RMSv2(folderpath, 2, 0);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void RMSv2(string folderpath, int nTF, int Algorithm) // number of trace files
{
	//OverUtilFactor = 1.0; // This factor over utilize the resource to get any more advantage of the system

	/////////////////////////////////// SMART SPACE ALLOCATION ////////////////////////////////////
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SSDspace = ((tracestats[itf].space < ssd.space * tracestats[itf].spaceratio) ? tracestats[itf].space : ssd.space * tracestats[itf].spaceratio);
		tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
	}
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////

	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 0;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF;			 itf++) {
		if(traces[itf][itr].lba < tracestats[itf].SSDspace) {
			(hit[itf])++; }
		else
			(miss[itf])++; }

	ldbl HDDdemandediops = 0;
	ldbl SSDdemandediops = 0;
	for(int itf = 0; itf < nTF; itf++)
	{
		HDDdemandediops += tracestats[itf].iops * miss[itf] / MaxTraceNum;
		SSDdemandediops += tracestats[itf].iops * hit [itf] / MaxTraceNum;
	}
	ldbl hEQC = (hdd.iops - HDDdemandediops) / (nTF * hdd.iops - HDDdemandediops); // HDD Equilibrium Constant
	ldbl sEQC = (ssd.iops - SSDdemandediops) / (nTF * ssd.iops - SSDdemandediops); // SSD Equilibrium Constant
	//ldbl HDDrestofiopsperApp = (hdd.iops > HDDdemandediops)?(hdd.iops - HDDdemandediops) / nTF:0;
	//ldbl SSDrestofiopsperApp = (ssd.iops > SSDdemandediops)?(ssd.iops - SSDdemandediops) / nTF:0;

	/////////////////////////////////////////////
	for(int itf = 0; itf < nTF; itf++)
	{
		ldbl hResidualService = (hEQC * (hdd.iops - tracestats[itf].iops * miss[itf] / MaxTraceNum));
		ldbl sResidualService = (sEQC * (ssd.iops - tracestats[itf].iops * hit [itf] / MaxTraceNum));
		tracestats[itf].HDDiops = (hResidualService > 0) ? (tracestats[itf].iops * miss[itf] / MaxTraceNum + hResidualService) : (tracestats[itf].iops * miss[itf] / MaxTraceNum);
		tracestats[itf].SSDiops = (sResidualService > 0) ? (tracestats[itf].iops * hit [itf] / MaxTraceNum + sResidualService) : (tracestats[itf].iops * hit [itf] / MaxTraceNum);
	}
	/////////////////////////////////////////////

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;
		printf("RMS:	App[%d]: HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//rms_hybrid_hdd.trace", folderpath + "//rms_hybrid_ssd.trace", nTF, Algorithm);
}

void BAAeq(string folderpath, int nTF, int Algorithm) // number of trace files
{
	//OverUtilFactor = 1.0; // This factor over utilize the resource to get any more advantage of the system

	// Parameters of HDD and SSD
	// HDD/SSD size allocation
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SSDspace = ((tracestats[itf].space < ssd.space * tracestats[itf].spaceratio) ? tracestats[itf].space : ssd.space * tracestats[itf].spaceratio);
		tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
	}
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	// Get hit/miss ratios based on size allocation
	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 0;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF; 		itf++) {
		if(traces[itf][itr].lba < tracestats[itf].HDDspace) {
			(miss[itf])++; }
		else
			(hit[itf])++; }

	// IOPS allocation
	double iops_hitbal = ssd.iops / (ssd.iops + hdd.iops);
	for(int itf = 0; itf < nTF; itf++) {
		if( (hit[itf] / (hit[itf] + miss[itf])) <= iops_hitbal) {
			//HDD iops fit
			//cout << "BAA: App[" << itf << "] is bottlenecked on HDD IOPS!\n";
			tracestats[itf].HDDiops = ((tracestats[itf].iops * nTF < hdd.iops) ? tracestats[itf].iops : hdd.iops / nTF);
			tracestats[itf].SSDiops = tracestats[itf].HDDiops * hit[itf] / miss[itf]; }
		else {
			//SSD iops fit
			//cout << "BAA: App[" << itf << "] is bottlenecked on SSD IOPS!\n";
			tracestats[itf].SSDiops = ((tracestats[itf].iops * nTF < ssd.iops) ? tracestats[itf].iops : ssd.iops / nTF);
			tracestats[itf].HDDiops = tracestats[itf].SSDiops * miss[itf] / hit[itf]; }
		//cout << "BAA: App[" << itf << "]: SSDiops = " << tracestats[itf].SSDiops << " and HDDiops = " << tracestats[itf].HDDiops << "\n";
	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;

		//double alfa = 0.8; // 0.1 0.2 0.3 ... 1
		//if(itf == 0)
		//{
		//	tracestats[itf].HDDiops += (182 - 132.7) * alfa;
		//	tracestats[itf].SSDiops += (1009.2 - 65) * alfa;
		//}
		//else
		//{
		//	tracestats[itf].HDDiops += (197 - 132.7) * alfa;
		//	tracestats[itf].SSDiops += (991.4 - 76.5) * alfa;
		//}

		printf("BAA:	App[%d]: HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//baa_hybrid_hdd.trace", folderpath + "//baa_hybrid_ssd.trace", nTF, Algorithm);
}

void EQAL_SPARE(string folderpath, int nTF, int Algorithm) // number of trace files
{
	//OverUtilFactor = 1.0; // This factor over utilize the resource to get any more advantage of the system

	//	ull totalspace = 0;
//	for(int itf = 0; itf < nTF; itf++) { totalspace += tracestats[itf].space; }
//	for(int itf = 0; itf < nTF; itf++) {
//		tracestats[itf].SSDspace = ((totalspace < ssd.space) ? tracestats[itf].space : ssd.space * tracestats[itf].space / totalspace); //( (ssd.space < tracestats[itf].space) ? ssd.space : tracestats[itf].space );
//		tracestats[itf].HDDspace = tracestats[itf].space - tracestats[itf].SSDspace;	}
	for(int itf = 0; itf < nTF; itf++)
	{
		//if( (tracestats[itf].space/totalspace) <= space_hitbal ) { // (hit[itf] / (hit[itf] + miss[itf]))
		//if( itf%2 ) { //0.5 <= space_hitbal[itf] ) {
			tracestats[itf].SSDspace = ((tracestats[itf].space < ssd.space * tracestats[itf].spaceratio) ? tracestats[itf].space : ssd.space * tracestats[itf].spaceratio);
			tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
		//} else {
			//tracestats[itf].HDDspace = ((tracestats[itf].space * nTF < hdd.space) ? tracestats[itf].space : hdd.space / nTF);
			//tracestats[itf].SSDspace = ((tracestats[itf].space > tracestats[itf].HDDspace) ? tracestats[itf].space - tracestats[itf].HDDspace : 0); // * miss[itf] / hit[itf];
		//}
		//cout << "BAA: App[" << itf << "]: SSDspace = " << tracestats[itf].SSDspace << " and HDDspace = " << tracestats[itf].HDDspace << "\n";
	}
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 0;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF;			 itf++) {
		if(traces[itf][itr].lba < tracestats[itf].SSDspace) {
			(hit[itf])++; }
		else
			(miss[itf])++; }

	ldbl HDDdemandediops = 0;
	ldbl SSDdemandediops = 0;
	for(int itf = 0; itf < nTF; itf++)
	{
		HDDdemandediops += tracestats[itf].iops * miss[itf] / MaxTraceNum;
		SSDdemandediops += tracestats[itf].iops * hit [itf] / MaxTraceNum;
	}
	ldbl HDDrestofiopsperApp = (hdd.iops > HDDdemandediops)?(hdd.iops - HDDdemandediops) / nTF:0;
	ldbl SSDrestofiopsperApp = (ssd.iops > SSDdemandediops)?(ssd.iops - SSDdemandediops) / nTF:0;

	for(int itf = 0; itf < nTF; itf++) {
		tracestats[itf].HDDiops = (tracestats[itf].iops * miss[itf] / MaxTraceNum) + HDDrestofiopsperApp;
		tracestats[itf].SSDiops = (tracestats[itf].iops * hit [itf] / MaxTraceNum) + SSDrestofiopsperApp; }

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;
		printf("EQS:	App[%d]: HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//eqs_hybrid_hdd.trace", folderpath + "//eqs_hybrid_ssd.trace", nTF, Algorithm);
}

void DRFv2(string folderpath, int nTF, int Algorithm) // number of trace files
{// Parameters of HDD and SSD
	// Rearrange input I/O traces based on DRF
	// (1) SSDspace0 + SSDspace1 = ssd.space (2) SSDspace0 / ssd.space = HDDspace1 / hdd.space
	// Note: we have two valid allocations, but I chose the more fair one!
	//if(tracestats[1].nbs > tracestats[0].nbs) {
	//	tracestats[1].SSDspace = ((ull) (ssd.space * (hdd.space - tracestats[1].nbs) / (hdd.space - ssd.space) ) ) / PageSize * PageSize;
	//	tracestats[1].HDDspace = tracestats[1].nbs - tracestats[1].SSDspace;
	//	tracestats[0].SSDspace = ssd.space - tracestats[1].SSDspace;
	//	tracestats[0].HDDspace = tracestats[0].nbs - tracestats[0].SSDspace; }
	//else {
	//	tracestats[0].SSDspace = ((ull) (ssd.space * (hdd.space - tracestats[0].nbs) / (hdd.space - ssd.space) ) ) / PageSize * PageSize;
	//	tracestats[0].HDDspace = tracestats[0].nbs - tracestats[0].SSDspace;
	//	tracestats[1].SSDspace = ssd.space - tracestats[0].SSDspace;
	//	tracestats[1].HDDspace = tracestats[1].nbs - tracestats[1].SSDspace; }

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SSDspace = ((tracestats[itf].space * nTF < ssd.space) ? tracestats[itf].space : ssd.space / nTF);
		tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
	}
	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}


	if(tracestats[1].iops < tracestats[0].iops)
	{
	//if(tracestats[1].bw > tracestats[0].bw) {
		tracestats[1].HDDiops = hdd.iops * (ssd.iops - tracestats[1].iops) / (ssd.iops - hdd.iops);
		tracestats[1].SSDiops = tracestats[1].iops - tracestats[1].HDDiops;
		tracestats[0].HDDiops = hdd.iops - tracestats[1].HDDiops; //(hdd.iops > tracestats[1].HDDiops) ? hdd.iops - tracestats[1].HDDiops : 0;
		tracestats[0].SSDiops = tracestats[0].iops - tracestats[0].HDDiops;
	}
	else
	{
		tracestats[0].HDDiops = hdd.iops * (ssd.iops - tracestats[0].iops) / (ssd.iops - hdd.iops);
		tracestats[0].SSDiops = tracestats[0].iops - tracestats[0].HDDiops;
		tracestats[1].HDDiops = hdd.iops - tracestats[0].HDDiops; //(hdd.iops > tracestats[0].HDDiops) ? hdd.iops - tracestats[0].HDDiops : 0;
		tracestats[1].SSDiops = tracestats[1].iops - tracestats[1].HDDiops;
	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;
		printf("DRF:	App[%d]:	HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//drf_hybrid_hdd.trace", folderpath + "//drf_hybrid_ssd.trace", nTF, Algorithm);
}



void ArrangeTraces(string hddtracefilepath, string ssdtracefilepath, int nTF, int Algorithm)
{
	/////////////////////////////////// INITIALIZE THE ARRAYS /////////////////////////////////////
	//int hddtrc [nTF][MaxTraceNum];
	//int ssdtrc [nTF][MaxTraceNum];
	TRACE** temptrc = new TRACE* [nTF];
	int** hddtrc = new int * [nTF]; // trace index in hdd for each app
	int** ssdtrc = new int * [nTF]; // trace index in ssd for each app
	for(int itf = 0; itf < nTF; itf++) {
		temptrc[itf]  = new TRACE [MaxTraceNum];
		hddtrc [itf]  = new int   [MaxTraceNum];
		ssdtrc [itf]  = new int   [MaxTraceNum];
		for(ull itr = 0; itr < MaxTraceNum; itr++) {
			hddtrc[itf][itr] = -1;
			ssdtrc[itf][itr] = -1;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// ARRANGE THE TRACES ////////////////////////////////////////
	ldbl hddtimeoffset;
	ldbl ssdtimeoffset;
	ldbl Ltime;
	ull hddsizeoffset = 0;
	ull ssdsizeoffset = 0;
	//ldbl hddtimestep = 1000 / hdd.iops;
	//ldbl ssdtimestep = 1000 / ssd.iops;
	unsigned int hddi[nTF]; // index of hdd for each app from 0 to max trace number in hdd
	unsigned int ssdi[nTF]; // index of ssd for each app from 0 to max trace number in ssd
	//unsigned int SpaceFreqThreshold[nTF];
	for(int itf = 0; itf < nTF; itf++)
	{
		hddi[itf] = 0; ssdi[itf] = 0;
	}
	for(int itf = 0; itf < nTF; itf++)
	{
		Ltime = 0;
		ssdtimeoffset = 0;
		hddtimeoffset = 0;
		for(ull itr = 0; (itr < MaxTraceNum) && (Ltime < MaxTime); itr++) // May be MaxTime
		{
			if(Algorithm == 0) // Basic
			{
				if(traces[itf][itr].lba < tracestats[itf].SSDspace)//+ traces[itf][itr].nbs < tracestats[itf].SSDspace)
				{// SSD trace
					temptrc[itf][itr].at = ((traces[itf][itr].at < ssdtimeoffset) ? ssdtimeoffset : traces[itf][itr].at);
					temptrc[itf][itr].dn = traces[itf][itr].dn;
					temptrc[itf][itr].lba = traces[itf][itr].lba + ssdsizeoffset;
					temptrc[itf][itr].nbs = traces[itf][itr].nbs;
					temptrc[itf][itr].md = traces[itf][itr].md;
					//if(!(temptrc[itf][itr].lba || temptrc[itf][itr].nbs))
					//	printf("Both are zero!");
					//ssdtimeoffset += ssdtimestep;
					ssdtimeoffset += 1000 / tracestats[itf].SSDiops;
					Ltime = temptrc[itf][itr].at;
					ssdtrc[itf][ssdi[itf]] = itr;
					(ssdi[itf])++;
				}
				else if (traces[itf][itr].lba >= tracestats[itf].SSDspace)
				{// HDD trace
					temptrc[itf][itr].at = ((traces[itf][itr].at < hddtimeoffset) ? hddtimeoffset : traces[itf][itr].at);
					temptrc[itf][itr].dn = traces[itf][itr].dn;
					temptrc[itf][itr].lba = traces[itf][itr].lba - tracestats[itf].SSDspace + hddsizeoffset;
					temptrc[itf][itr].nbs = traces[itf][itr].nbs;
					temptrc[itf][itr].md = traces[itf][itr].md;
					//if(!(temptrc[itf][itr].lba || temptrc[itf][itr].nbs))
					//	printf("Both are zero!");
					//hddtimeoffset += hddtimestep;
					hddtimeoffset += 1000 / tracestats[itf].HDDiops;
					Ltime = temptrc[itf][itr].at;
					hddtrc[itf][hddi[itf]] = itr;
					(hddi[itf])++;
				}
			}
			else if(Algorithm == 1) // Smart
			{
				if( tracestats[itf].SpaceDist[ traces[itf][itr].lba ] >= tracestats[itf].SpaceFreqThreshold ) // test also tracestats[itf].SpaceFreqThreshold - 1
				//if(traces[itf][itr].lba < tracestats[itf].SSDspace)//+ traces[itf][itr].nbs < tracestats[itf].SSDspace)
				{// SSD trace
					temptrc[itf][itr].at = ((traces[itf][itr].at < ssdtimeoffset) ? ssdtimeoffset : traces[itf][itr].at);
					temptrc[itf][itr].dn = traces[itf][itr].dn;
					temptrc[itf][itr].lba = (traces[itf][itr].lba % tracestats[itf].SSDspace) + ssdsizeoffset;
					temptrc[itf][itr].nbs = traces[itf][itr].nbs;
					temptrc[itf][itr].md = traces[itf][itr].md;
					//if(!(temptrc[itf][itr].lba || temptrc[itf][itr].nbs))
					//	printf("Both are zero!");
					//ssdtimeoffset += ssdtimestep;
					ssdtimeoffset += 1000 / tracestats[itf].SSDiops;
					Ltime = temptrc[itf][itr].at;
					ssdtrc[itf][ssdi[itf]] = itr;
					(ssdi[itf])++;
				}
				else if ( tracestats[itf].SpaceDist[ traces[itf][itr].lba ] < tracestats[itf].SpaceFreqThreshold )
				//else if (traces[itf][itr].lba >= tracestats[itf].SSDspace)
				{// HDD trace
					temptrc[itf][itr].at = ((traces[itf][itr].at < hddtimeoffset) ? hddtimeoffset : traces[itf][itr].at);
					temptrc[itf][itr].dn = traces[itf][itr].dn;
					temptrc[itf][itr].lba = traces[itf][itr].lba + hddsizeoffset;
					temptrc[itf][itr].nbs = traces[itf][itr].nbs;
					temptrc[itf][itr].md = traces[itf][itr].md;
					//if(!(temptrc[itf][itr].lba || temptrc[itf][itr].nbs))
					//	printf("Both are zero!");
					//hddtimeoffset += hddtimestep;
					hddtimeoffset += 1000 / tracestats[itf].HDDiops;
					Ltime = temptrc[itf][itr].at;
					hddtrc[itf][hddi[itf]] = itr;
					(hddi[itf])++;
				}
			}
		}
		ssdsizeoffset += tracestats[itf].SSDspace;
		hddsizeoffset += tracestats[itf].HDDspace;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// WRITE TO HDD/SSD TRACE FILES //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// WRITE TO HDD TRACE FILES //////////////////////////////////
	// HDD
	//string hddtracepath = folderpath + "//hybrid_hdd.trace";
	ofstream hybridhdd (&hddtracefilepath[0]);
	hybridhdd.precision(20);

	unsigned int hi[nTF]; // count hddi
	unsigned int hsum = 0;
	unsigned int hsumi = 0;
	ldbl hddmaxat = -1;
	for(int itf = 0; itf < nTF; itf++)
	{
		if(hddi[itf])
			hi[itf] = 1;
		else // no hdd trace!
			hi[itf] = 0;
		hsum+=(hddi[itf]);
		if((hddi[itf] > 0) && (hddmaxat < temptrc[itf][hddtrc[itf][hddi[itf]-1]].at))
			hddmaxat = temptrc[itf][hddtrc[itf][hddi[itf]-1]].at;
	}
	int iTR;
	int iTF;
	int itf = 0;
	ldbl hddminat;
	//int whilecond = 1;
	//while(whilecond)
	while(hsumi < hsum)
	{
		hddminat = hddmaxat; //temptrc[iTF][hddtrc[iTF][hi[iTF]]].at;
		iTF = -1;
		for(itf = 0; itf < nTF; itf++)
			if( (hi[itf] > 0) && (hi[itf] <= hddi[itf]) && (temptrc[itf][hddtrc[itf][hi[itf]-1]].at <= hddminat) )
			{
				hddminat = temptrc[itf][hddtrc[itf][hi[itf]-1]].at;
				iTF = itf;
			}
		if(iTF >= 0)
		{
			iTR = hddtrc[iTF][hi[iTF]-1];
			hybridhdd << iTF
					<< "	" << temptrc[iTF][iTR].at
					<< "	" << 0 //temptrc[iTF][iTR].dn
					<< "	" << temptrc[iTF][iTR].lba
					<< "	" << temptrc[iTF][iTR].nbs
					<< "	" << (((char)(temptrc[iTF][iTR].md)=='R')?1:0)
					<< "\n";
			(hi[iTF])++;
			//if(hi[iTF] >= hddi[iTF])
			//	whilecond = 0;
			hsumi++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// WRITE TO SSD TRACE FILES //////////////////////////////////
	// SSD
	//string ssdtracepath = folderpath + "//hybrid_ssd.trace";
	ofstream hybridssd (&ssdtracefilepath[0]);
	hybridssd.precision(20);

	unsigned int ssumi = 0;
	unsigned int ssum = 0;
	unsigned int si[nTF]; // count ssdi
	ldbl ssdmaxat = -1;
	for(int itf = 0; itf < nTF; itf++)
	{
		if(ssdi[itf])
			si[itf] = 1;
		else // no ssd trace!
			si[itf] = 0;
		ssum+=(ssdi[itf]);
		if((ssdi[itf] > 0) && (ssdmaxat < temptrc[itf][ssdtrc[itf][ssdi[itf]-1]].at))
			ssdmaxat = temptrc[itf][ssdtrc[itf][ssdi[itf]-1]].at;
	}
	ldbl ssdminat;
	//whilecond = 1;
	//while(whilecond)
	while(ssumi < ssum)
	{
		ssdminat = ssdmaxat; //temptrc[iTF][ssdtrc[iTF][hi[iTF]]].at;
		iTF = -1;
		for(int itf = 0; itf < nTF; itf++)
			if( (si[itf] > 0) && (si[itf] <= ssdi[itf]) && (temptrc[itf][ssdtrc[itf][si[itf]-1]].at <= ssdminat) )
			{
				ssdminat = temptrc[itf][ssdtrc[itf][si[itf]-1]].at;
				iTF = itf;
			}
		if(iTF >= 0)
		{
			iTR = ssdtrc[iTF][si[iTF]-1];
			ull new_bn = (((temptrc[iTF][iTR].lba) / PageSize) * PageSize);
			ull new_space = (((temptrc[iTF][iTR].lba + temptrc[iTF][iTR].nbs + PageSize - 1) / PageSize * PageSize) - new_bn);
			hybridssd << iTF
					<< "	" << temptrc[iTF][iTR].at
					<< "	" << 0 //temptrc[iTF][iTR].dn
					<< "	" << new_bn
					<< "	" << new_space
					<< "	" << (((char)(temptrc[iTF][iTR].md)=='R')?1:0)
					<< "\n";
			(si[iTF]) ++;
			//if(si[iTF] >= ssdi[iTF])
			//	whilecond = 0;
			ssumi++;
		}
	}
	/////////////////////////////////// DELETE ALL TEMP ARRAYS ////////////////////////////////////
	for(int i = 0; i < nTF; i++) {
		delete [] (temptrc[i]);
		delete [] (hddtrc[i]);
		delete [] (ssdtrc[i]);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
}

void SetTrace(string tracefile, string newtracefile)
{
	ifstream logfile    (&   tracefile[0]);
	ofstream newlogfile (&newtracefile[0]);
	logfile.precision(20);
	newlogfile.precision(20);
	TRACE tr;
	if (logfile.is_open())
	{
		int i = 0, j = 0;
		getline (logfile,line);
		ldbl BaseIAT = GetNextNumber(&line,&i);
		do
		{
			i = 0;
			tr.at = (ull)(floor((GetNextNumber(&line,&i) - BaseIAT) / 10)); // 100ns -> us
			tr.dn = GetNextNumber(&line,&i);
			tr.fd = 0;
			tr.md = GetNextNumber(&line,&i); //((line[i+1] == 'R')?(unsigned)('R'):(unsigned)('W'));
			tr.lba = (ull)(GetNextNumber(&line,&i) / 4096);
			tr.nbs = (ull)(GetNextNumber(&line,&i) / 4096);
			newlogfile << j << "	" << tr.at << "	" << tr.dn << "	"
					<< tr.fd << "	" << (((char)(tr.md)=='R')?"Read":"Write") << "	"
					<< tr.lba << "	" << tr.nbs << "\n";
			j++;
		} while (getline (logfile,line));
	}
	logfile.close();
	newlogfile.close();
}

void GetTraceStat(string tracefile, int iTF, int warmupnum, ldbl spaceratio)
{
	ifstream logfile;
	logfile.open(&tracefile[0]); // command line program gets log filename
	ull tByte = 0;
	ull tIO = 0;
	ull minLBA = 1000000000; //traces[0].lba;
	ull maxLBA = 0; //traces[0].lba + traces[0].nbs;
	ull minNBS = 1000000000;
	ull maxNBS = 0;
	//int wunum = WARMUPNUM; // warm up
	ull ti = 0;
	if (logfile.is_open())
	{
		int i = 0;
		while(warmupnum) {
			getline (logfile,line); // warmup
			warmupnum--;
		}
		getline (logfile,line);
		GetNextNumber(&line,&i);
		ldbl BaseIAT = GetNextNumber(&line,&i) / 1000; // us -> ms or micro-sec to mili-sec
		do
		{
			i = 0;
			GetNextNumber(&line,&i); // Just to skip the row number
			traces[iTF][ti].at = (GetNextNumber(&line,&i) - BaseIAT) / 1000 / AT_SqueezeFactor; // us -> ms
			traces[iTF][ti].dn = GetNextNumber(&line,&i);
			traces[iTF][ti].fd = GetNextNumber(&line,&i);
			traces[iTF][ti].md = GetNextNumber(&line,&i);
			traces[iTF][ti].lba = (ull)(GetNextNumber(&line,&i) / BN_SqueezeFactor);
			traces[iTF][ti].nbs = (ull)(GetNextNumber(&line,&i));
			tIO++;
			if( (ti < MaxTraceNum) && (traces[iTF][ti].at < MaxTime) && (traces[iTF][ti].nbs) )// sometimes space == 0 !!!
			{
				// min Block number
				if(traces[iTF][ti].lba < minLBA)
					minLBA = traces[iTF][ti].lba;
				// max Block number
				if(traces[iTF][ti].lba + traces[iTF][ti].nbs > maxLBA)
					maxLBA = traces[iTF][ti].lba + traces[iTF][ti].nbs;
				// min size
				if(traces[iTF][ti].nbs < minNBS)
					minNBS = traces[iTF][ti].nbs;
				// max size
				if(traces[iTF][ti].nbs > maxNBS)
					maxNBS = traces[iTF][ti].nbs;
				// total bytes
				tByte += traces[iTF][ti].nbs;
				// space distribution set
				for(ull bnc = traces[iTF][ti].lba; bnc < traces[iTF][ti].lba + traces[iTF][ti].nbs; bnc++)
					(tracestats[iTF].SpaceDist[bnc])++;
				// increment the counter
				ti++;
			}
			//else
			//	printf("is zero!\n");
		} while ( (tIO < MaxTraceSize) && getline (logfile,line) );
		tracestats[iTF].TraceSize = ti;
	}
	logfile.close();

	//cout << (ldbl) ((1000 * tSize) / (traces[tracestats[iTF].TraceSize - 1].at - traces[0].at)) << "\n";
	//cout << (ldbl) ((1000 * tSize) / (traces[tracestats[iTF].TraceSize - 1].dt - traces[0].dt)) << "\n";

	// Primary Requirements
	tracestats[iTF].spaceratio = spaceratio;
	tracestats[iTF].space = maxLBA; // - min + 1;
	tracestats[iTF].MaxNBS = maxNBS;
	//tracestats[iTF].bw = (ldbl) (tByte * (1000 / (traces[iTF][tracestats[iTF].TraceSize - 1].at - traces[iTF][0].at))); // Bytes/s
	tracestats[iTF].f_time = (ldbl) (traces[iTF][tracestats[iTF].TraceSize - 1].at);
	tracestats[iTF].iops = (ldbl) (tracestats[iTF].TraceSize * (1000 / (traces[iTF][tracestats[iTF].TraceSize - 1].at - traces[iTF][0].at))); // IO/s
	tracestats[iTF].bw = (ldbl) (tByte * (1000 / (traces[iTF][tracestats[iTF].TraceSize].at - traces[iTF][0].at))); // Bytes/s
	tracestats[iTF].iops_t = (ldbl) (tIO * (1000 / (traces[iTF][tracestats[iTF].TraceSize].at - traces[iTF][0].at))); // IO/s
	cout << "App[" << iTF << "]: MaxBlockNum = " << tracestats[iTF].space
			<< ",	LastTime = " << tracestats[iTF].f_time
			<< ",	IOPS = " << tracestats[iTF].iops
			<< ",	Total BW = " << tracestats[iTF].bw
			<< ",	Total IOPS = " << tracestats[iTF].iops_t
			<< ",	Max Number of Blocks = " << tracestats[iTF].MaxNBS
			<< "\n";
	//tracestats[index].rt;
}

ldbl GetMaxThroughput(string tracefile, int iTF, int warmupnum)
{
	ifstream logfile (&tracefile[0]); // command line program gets log filename
	int itr = 0;
	ull tByte = 0;
	ull tIO = 0;
	//TRACE traces[1000000];//[MaxTraceNum];
	if (logfile.is_open())
	{
		while(warmupnum) {
			getline (logfile,line); // warmup
			warmupnum--;
		}
		while ( ( itr < MaxTraceNum ) && getline (logfile,line) )
		{
			int i = 0;
			traces[iTF][itr].id = GetNextNumber(&line, &i);
			traces[iTF][itr].at = GetNextNumber(&line, &i);
			traces[iTF][itr].dn = GetNextNumber(&line, &i);
			traces[iTF][itr].lba = GetNextNumber(&line, &i);
			traces[iTF][itr].nbs = GetNextNumber(&line, &i);
			traces[iTF][itr].md = GetNextNumber(&line, &i);
			traces[iTF][itr].dt = GetNextNumber(&line, &i);
			//if( traces[iTF][itr].at < 1) // just for debug
			//	traces[iTF][itr].at = 0; // just for debug
			tByte += traces[iTF][itr].nbs;
			tIO ++;
			itr++;
		}
	}
	logfile.close();
	//ldbl bps = (ldbl) ((1000 * tByte) / (traces[iTF][itr-1].dt - traces[iTF][0].at));
	ldbl iops = (ldbl) ((1000 * tIO) / (traces[iTF][itr-1].dt - traces[iTF][0].at));
	return iops;
}

ldbl GetNextNumber(string* Line, int* iL)
{
	//string line = (*Line);
	ldbl num = 0;
	unsigned mode = 0;
	while(  ((*iL)<(int)(Line->length())) &&
			( 		   (((int)((*Line)[*iL])) < (int)'0' )
					|| (((int)((*Line)[*iL])) > (int)'9' )
					//|| (((int)((*Line)[*iL])) == 'W')
					//|| (((int)((*Line)[*iL])) == 'R')
					) )
	{
		if( ((int)((*Line)[*iL])) == 'W')
			mode = (unsigned)('W');
		if( ((int)((*Line)[*iL])) == 'R')
			mode = (unsigned)('R');
		(*iL)++;
	}
	if(mode)
		return mode;

	while(  ((*iL)<(int)(Line->length()))
			&& ( ((int)((*Line)[*iL])) >= (int)'0' ) && (((int)((*Line)[*iL])) <= (int)'9' ) )
	{
		num *= 10;
		num += (ldbl)((unsigned)((*Line)[*iL])) - (unsigned)('0');
		(*iL)++;
	}
	if( ((int)((*Line)[*iL])) == '.' ) {
		(*iL)++;
		ldbl dec = 1;
		while(  ((*iL)<(int)(Line->length()))
				&& ( ((int)((*Line)[*iL])) >= (int)'0' ) && (((int)((*Line)[*iL])) <= (int)'9' ) )
		{
			dec /= 10;
			num += (dec * ((ldbl)((unsigned)((*Line)[*iL])) - (unsigned)('0')));
			(*iL)++;
		}
	}
	return num;
}

void RAMTSv1(string folderpath, int nTF, int Algorithm) // number of trace files
{
	ull totalspace = 0;
	for(int itf = 0; itf < nTF; itf++) { totalspace += tracestats[itf].space; }
	for(int itf = 0; itf < nTF; itf++) {
		tracestats[itf].SSDspace = ((totalspace < ssd.space) ? tracestats[itf].space : ssd.space * tracestats[itf].space / totalspace); //( (ssd.space < tracestats[itf].space) ? ssd.space : tracestats[itf].space );
		tracestats[itf].HDDspace = tracestats[itf].space - tracestats[itf].SSDspace;	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 0;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF;			 itf++) {
		if(traces[itf][itr].lba < tracestats[itf].SSDspace) {
			(hit[itf])++; }
		else
			(miss[itf])++; }

	ldbl HDDdemandediops = 0;
	ldbl SSDdemandediops = 0;
	for(int itf = 0; itf < nTF; itf++) {
		HDDdemandediops += tracestats[itf].iops * miss[itf] / MaxTraceNum;
		SSDdemandediops += tracestats[itf].iops * hit [itf] / MaxTraceNum; }
	ldbl HDDrestofiopsperApp = (hdd.iops > HDDdemandediops)?(hdd.iops - HDDdemandediops) / nTF:0;
	ldbl SSDrestofiopsperApp = (ssd.iops > SSDdemandediops)?(ssd.iops - SSDdemandediops) / nTF:0;

	for(int itf = 0; itf < nTF; itf++) {
		tracestats[itf].HDDiops = (tracestats[itf].iops * miss[itf] / MaxTraceNum) + HDDrestofiopsperApp;
		tracestats[itf].SSDiops = (tracestats[itf].iops * hit [itf] / MaxTraceNum) + SSDrestofiopsperApp; }

	for(int itf = 0; itf < nTF; itf++) {
		printf("RAMTS1:	App[%d]: HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops); }

	ArrangeTraces(folderpath + "//ramts1_hybrid_hdd.trace", folderpath + "//ramts1_hybrid_ssd.trace", nTF, Algorithm);
}

void BAAv1(string folderpath, int nTF, int Algorithm) // number of trace files
{// Parameters of HDD and SSD

	// HDD/SSD size allocation
	ull totalspace = 0;
	for(int itf = 0; itf < nTF; itf++) { totalspace += tracestats[itf].space; }
	ldbl space_hitbal = (ldbl)hdd.space / (ldbl)(ssd.space + hdd.space);
	for(int itf = 0; itf < nTF; itf++) {
		//if( (tracestats[itf].space/totalspace) <= space_hitbal ) { // (hit[itf] / (hit[itf] + miss[itf]))
		if( 0.5 <= space_hitbal ) {
			//cout << "BAA: App[" << itf << "] is bottlenecked on SSD size!\n";
			tracestats[itf].SSDspace = ((tracestats[itf].space * nTF < ssd.space) ? tracestats[itf].space : ssd.space / nTF);
			tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
		} else {
			//cout << "BAA: App[" << itf << "] is bottlenecked on HDD size!\n";
			tracestats[itf].HDDspace = ((tracestats[itf].space * nTF < hdd.space) ? tracestats[itf].space : hdd.space / nTF);
			tracestats[itf].SSDspace = ((tracestats[itf].space > tracestats[itf].HDDspace) ? tracestats[itf].space - tracestats[itf].HDDspace : 0); // * miss[itf] / hit[itf];
		}
		//cout << "BAA: App[" << itf << "]: SSDspace = " << tracestats[itf].SSDspace << " and HDDspace = " << tracestats[itf].HDDspace << "\n";
	}

	//for(int itf = 0; itf < nTF; 		itf++) {
		//tracestats[itf].HDDspace = tracestats[itf].space * hdd.space / (hdd.space + ssd.space) / 2; // size ratio ???
		//tracestats[itf].HDDspace = tracestats[itf].space * sqrt(hdd.space) / (sqrt(hdd.space) + sqrt(ssd.space)); // size ratio ???
		//tracestats[itf].SSDspace = tracestats[itf].space - tracestats[itf].HDDspace;	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	// Get hit/miss ratios based on size allocation
	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 1;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF; 		itf++) {
		if(traces[itf][itr].lba < tracestats[itf].HDDspace) {
			(miss[itf])++; }
		else
			(hit[itf])++; }

	// IOPS allocation
	double iops_hitbal = ssd.iops / (ssd.iops + hdd.iops);
	for(int itf = 0; itf < nTF; itf++) {
		if( (hit[itf] / (hit[itf] + miss[itf])) <= iops_hitbal) {
			//cout << "BAA: App[" << itf << "] is bottlenecked on HDD IOPS!\n";
			tracestats[itf].HDDiops = ((tracestats[itf].iops * nTF < hdd.iops) ? tracestats[itf].iops : hdd.iops / nTF);
			tracestats[itf].SSDiops = tracestats[itf].HDDiops * hit[itf] / miss[itf]; }
		else {
			//cout << "BAA: App[" << itf << "] is bottlenecked on SSD IOPS!\n";
			tracestats[itf].SSDiops = ((tracestats[itf].iops * nTF < ssd.iops) ? tracestats[itf].iops : ssd.iops / nTF);
			tracestats[itf].HDDiops = tracestats[itf].SSDiops * miss[itf] / hit[itf]; }
	}

	for(int itf = 0; itf < nTF; itf++) {
		printf("BAA:	App[%d]:	HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops); }

	ArrangeTraces(folderpath + "//baa_hybrid_hdd.trace", folderpath + "//baa_hybrid_ssd.trace", nTF, Algorithm);
}

void DRFv1(string folderpath, int nTF, int Algorithm) // number of trace files
{// Parameters of HDD and SSD
	// Rearrange input I/O traces based on DRF
	// (1) SSDspace0 + SSDspace1 = ssd.space (2) SSDspace0 / ssd.space = HDDspace1 / hdd.space
	// Note: we have two valid allocations, but I chose the more fair one!
	if(tracestats[1].space > tracestats[0].space) {
		tracestats[1].SSDspace = ((ull) (ssd.space * (hdd.space - tracestats[1].space) / (hdd.space - ssd.space) ) ) / PageSize * PageSize;
		tracestats[1].HDDspace = tracestats[1].space - tracestats[1].SSDspace;
		tracestats[0].SSDspace = ssd.space - tracestats[1].SSDspace;
		tracestats[0].HDDspace = tracestats[0].space - tracestats[0].SSDspace; }
	else {
		tracestats[0].SSDspace = ((ull) (ssd.space * (hdd.space - tracestats[0].space) / (hdd.space - ssd.space) ) ) / PageSize * PageSize;
		tracestats[0].HDDspace = tracestats[0].space - tracestats[0].SSDspace;
		tracestats[1].SSDspace = ssd.space - tracestats[0].SSDspace;
		tracestats[1].HDDspace = tracestats[1].space - tracestats[1].SSDspace; }

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	if(tracestats[1].iops > tracestats[0].iops) {
		tracestats[1].HDDiops = hdd.iops * (ssd.iops - tracestats[1].iops) / (ssd.iops - hdd.iops);
		tracestats[1].SSDiops = tracestats[1].iops - tracestats[1].HDDiops;
		tracestats[0].HDDiops = hdd.iops - tracestats[1].HDDiops;
		tracestats[0].SSDiops = tracestats[0].iops - tracestats[0].HDDiops; }
	else {
		tracestats[0].HDDiops = hdd.iops * (ssd.iops - tracestats[0].iops) / (ssd.iops - hdd.iops);
		tracestats[0].SSDiops = tracestats[0].iops - tracestats[0].HDDiops;
		tracestats[1].HDDiops = hdd.iops - tracestats[0].HDDiops;
		tracestats[1].SSDiops = tracestats[1].iops - tracestats[1].HDDiops; }

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;
		printf("DRF:	App[%d]:	HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//drf_hybrid_hdd.trace", folderpath + "//drf_hybrid_ssd.trace", nTF, Algorithm);
}

void BAAv2(string folderpath, int nTF, int Algorithm) // number of trace files
{
	//OverUtilFactor = 1.0; // This factor over utilize the resource to get any more advantage of the system

	// Parameters of HDD and SSD
	// HDD/SSD size allocation
	//ull totalspace = 0;
	//for(int itf = 0; itf < nTF; itf++) { totalspace += tracestats[itf].space; }
	//ldbl space_hitbal[nTF] = {0.4, 0.6}; //(ldbl)hdd.space / (ldbl)(ssd.space + hdd.space);
	//for(int itf = 0; itf < nTF; itf++) {
		//if( (tracestats[itf].space/totalspace) <= space_hitbal ) { // (hit[itf] / (hit[itf] + miss[itf]))
	//	if( itf%2 ) { //0.5 <= space_hitbal[itf] ) {
	//		tracestats[itf].SSDspace = ((tracestats[itf].space * nTF < ssd.space) ? tracestats[itf].space : ssd.space / nTF);
	//		tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * hit[itf] / miss[itf];
	//	} else {
	//		tracestats[itf].HDDspace = ((tracestats[itf].space * nTF < hdd.space) ? tracestats[itf].space : hdd.space / nTF);
	//		tracestats[itf].SSDspace = ((tracestats[itf].space > tracestats[itf].HDDspace) ? tracestats[itf].space - tracestats[itf].HDDspace : 0); // * miss[itf] / hit[itf];
	//	}
		//cout << "BAA: App[" << itf << "]: SSDspace = " << tracestats[itf].SSDspace << " and HDDspace = " << tracestats[itf].HDDspace << "\n";
	//}

	ldbl SumSize = 0;
	for(int itf = 0; itf < nTF; itf++)
		SumSize += (ldbl)(tracestats[itf].space);
	ldbl space_bal = SumSize / ((ldbl)(ssd.space + hdd.space)) * ((ldbl)(hdd.space));
	for(int itf = 0; itf < nTF; itf++) {
		if( ((ldbl)(tracestats[itf].space)) > space_bal)
		{
			//HDD space fit
			tracestats[itf].HDDspace = ((tracestats[itf].space * nTF < hdd.space) ? tracestats[itf].space : hdd.space / nTF);
			tracestats[itf].SSDspace = ((tracestats[itf].space > tracestats[itf].HDDspace) ? tracestats[itf].space - tracestats[itf].HDDspace : 0); // * miss[itf] / hit[itf];
		} else {
			//SSD space fit
			//cout << "BAA: App[" << itf << "] is bottlenecked on SSD IOPS!\n";
			tracestats[itf].SSDspace = ((tracestats[itf].space * nTF < ssd.space) ? tracestats[itf].space : ssd.space / nTF) * tracestats[itf].spaceratio;
			tracestats[itf].HDDspace = ((tracestats[itf].space > tracestats[itf].SSDspace) ? tracestats[itf].space - tracestats[itf].SSDspace : 0); // * miss[itf] / hit[itf];
		//cout << "BAA: App[" << itf << "]: SSDiops = " << tracestats[itf].SSDiops << " and HDDiops = " << tracestats[itf].HDDiops << "\n";
		}
	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].SpaceFreqThreshold = 0;
		ull sc;
		ull ss = 0;
		ull pre_ss =ss;
		while( ss + tracestats[itf].SSDspace < MAX_SPACE )
		{
			pre_ss = ss;
			for(sc = 0; sc < MAX_SPACE; sc++)
				if( tracestats[itf].SpaceDist[sc] == tracestats[itf].SpaceFreqThreshold )
					ss++;
			tracestats[itf].SpaceFreqThreshold++;
		}
		if( tracestats[itf].SSDspace > (ss - pre_ss) * 2 )
			tracestats[itf].SpaceFreqThreshold--;
	}

	// Get hit/miss ratios based on size allocation
	int miss[nTF]; // miss on SSD or hit on HDD
	int hit[nTF]; // hit on SSD or miss on HDD
	for(int i = 0; i < nTF; i++) { miss[i] = hit[i] = 0;}
	for(ull itr = 0; itr < MaxTraceNum; itr++)
	for(int itf = 0; itf < nTF; 		itf++) {
		if(traces[itf][itr].lba < tracestats[itf].HDDspace) {
			(miss[itf])++; }
		else
			(hit[itf])++; }

	// IOPS allocation
	double iops_hitbal = ssd.iops / (ssd.iops + hdd.iops);
	for(int itf = 0; itf < nTF; itf++) {
		if( (hit[itf] / (hit[itf] + miss[itf])) <= iops_hitbal) {
			//HDD iops fit
			//cout << "BAA: App[" << itf << "] is bottlenecked on HDD IOPS!\n";
			tracestats[itf].HDDiops = ((tracestats[itf].iops * nTF < hdd.iops) ? tracestats[itf].iops : hdd.iops / nTF);
			tracestats[itf].SSDiops = tracestats[itf].HDDiops * hit[itf] / miss[itf]; }
		else {
			//SSD iops fit
			//cout << "BAA: App[" << itf << "] is bottlenecked on SSD IOPS!\n";
			tracestats[itf].SSDiops = ((tracestats[itf].iops * nTF < ssd.iops) ? tracestats[itf].iops : ssd.iops / nTF);
			tracestats[itf].HDDiops = tracestats[itf].SSDiops * miss[itf] / hit[itf]; }
		//cout << "BAA: App[" << itf << "]: SSDiops = " << tracestats[itf].SSDiops << " and HDDiops = " << tracestats[itf].HDDiops << "\n";
	}

	for(int itf = 0; itf < nTF; itf++)
	{
		tracestats[itf].HDDiops *= UtilFactor;
		tracestats[itf].SSDiops *= UtilFactor;
		printf("BAA:	App[%d]: HDDspace = %llu, SSDspace = %llu, HDDiops = %Lf, SSDiops = %Lf \n", itf, tracestats[itf].HDDspace, tracestats[itf].SSDspace,
				tracestats[itf].HDDiops, tracestats[itf].SSDiops);
	}

	ArrangeTraces(folderpath + "//baa_hybrid_hdd.trace", folderpath + "//baa_hybrid_ssd.trace", nTF, Algorithm);
}
