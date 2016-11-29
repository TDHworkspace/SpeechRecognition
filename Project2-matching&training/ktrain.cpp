#include "ktrain.h"

vector<vector<vector<double> > > samples;

vector<vector<double> > means;
vector<vector<double> > covars;
vector<vector<double> > kdist;

double trans [SEGNUM+1][SEGNUM];
int segdiv [SAMPLENUM][SEGNUM+1];

void readSampleDCT(const char* filename, int flag){
	samples.resize(SAMPLENUM);
	int i;
	for (i = 0; i < SAMPLENUM; i++){
		samples[i].resize(39);
	}
	//Read MFCC Spectrum from file.
	ifstream in_DCT(filename);
	string temp;
	double DCT_value;
	int count = 0;
	if (!in_DCT.is_open()) {
		cout << "Open file failed! \n";
	}
	else {
		while (!in_DCT.eof()) {
			if (count >= 39){
				count = 0;
			}
		    in_DCT >> temp;
			DCT_value = std::stod(temp);
			samples[flag][count].push_back(DCT_value);
			count ++;
		}
	}
}

void initSeg()
{
	readSampleDCT("s-9-1.txt", 0);
    readSampleDCT("s-9-2.txt", 1);
	int i, j;
	int frameLen;
	for (i = 0; i <SAMPLENUM; i++){
		frameLen = samples[i][0].size();
		for (j = 0; j <SEGNUM; j++){
			segdiv[i][j] = (frameLen / 5) * (j);
			//cout << segdiv[i][j] << " ";
	    }
		segdiv[i][j] = frameLen;
		//cout << segdiv[i][j] << " ";
	}
	//cout << endl;
}

void initTrans(){
    int i, j, k, l;
	double count = 0.0;
	for (i = 0; i < SEGNUM + 1; i++){
		for(j = 0; j < SEGNUM; j++){
			if(i == 0){
				if (j == 0){
				    trans[i][j] = 0;
				}
				else{
					trans[i][j] = INF;
				}
			}
			else{
				if ((i - j) == 1 || (i == j)){

					for (k = 0; k < SAMPLENUM; k++){
						for (l = segdiv[k][i-1]; l < segdiv[k][i]; l++)
						{count++;}
					}

					if (( i - j ) == 1){
						trans[i][j] = -log((count-SAMPLENUM)/count);
					}
					else{
						trans[i][j] = -log(SAMPLENUM / count);
					}

					count = 0;

				}
				else{
					trans[i][j] = INF;
				}
			}
		}
		
	}
}

void calSegmentMeansCovar(){

    //initSegTrans();
	means.resize(39);
	covars.resize(39);
	double tmpMean = 0.0;
	int count = 0;
	int samp, dimens, frame, seg;
	//Calculate 5 means of each segment.
	for(seg = 0; seg < SEGNUM; seg ++){
	    //Calculate 39 means of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment means of all samples.
	        for(samp = 0; samp < SAMPLENUM; samp++){
			//Get the segment range.
			    for (frame = segdiv[samp][seg]; frame < segdiv[samp][seg+1]; frame++){
				    tmpMean += samples[samp][dimens][frame];
				    count++;
			    }
		    }
		    //cout << tmpMean /count << " " << endl;
		    means[dimens].push_back(tmpMean/count);
		    tmpMean = 0;
		    count = 0;
	    }
	}

    double tmpCovar = 0.0;

    //Calculate 5 covars of each segment.
	for(seg = 0; seg < SEGNUM; seg ++){
	    //Calculate 39 covars of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment covars of all samples.
	        for(samp = 0; samp < SAMPLENUM; samp++){
			//Get the segment range.
			    for (frame = segdiv[samp][seg]; frame < segdiv[samp][seg+1]; frame++){
				    tmpCovar += (samples[samp][dimens][frame] - means[dimens][seg]) * (samples[samp][dimens][frame] - means[dimens][seg]);
				    count++;
			    }
		    }
		    //cout << tmpMean /count << " " << endl;
		    covars[dimens].push_back(tmpCovar/count);
		    tmpCovar = 0;
		    count = 0;
	    }
	}
}

void calSampDistance (int num){
	//Warning! Call this function only after means and covar is fully initialized!
	ofstream disfile("distance.txt");

	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = samples[num][1].size();
	int meansize = means[1].size();

	kdist.resize(inputsize);
	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < meansize; j++)
		{
			for ( k = 0; k < 39; k++)
			{
				tmpDistance += log( 2 * 3.14159 * covars[k][j] ) + (samples[num][k][i] - means[k][j]) * (samples[num][k][i] - means[k][j]) / covars[k][j];
			}
			tmpDistance = 0.5 * tmpDistance;
			kdist[i].push_back(tmpDistance);
			//cout << tmpDistance << " ";
			tmpDistance = 0;
		}
	}
	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < inputsize; i++)
		{
			for (j = 0; j < meansize; j++)
			{
				disfile << kdist[i][j] << ' ';
			}
			disfile << endl;
		}
	}
	//cout << kdist.size() << endl;
}


double kmin3(double superdiag, double horizon, double diag){
	double temp = (superdiag < horizon) ? superdiag : horizon;
	return (temp < diag) ? temp : diag;
}

double kmin2(double horizon, double diag){
	return (horizon < diag) ? horizon : diag;
}


void kmeansDistance(){
	//Warning! Call this function only after kdist vector is fully initialized!
	ofstream disfile("keditdistance.txt");
	int insize = kdist.size();
	int tmpsize = kdist[0].size();
	int i, j;
	for (i = 0; i < insize; i++){
		for(j = 0; j < tmpsize; j++){
			//If i!= 0 and j!= 0;
			if (!(i == 0 && j == 0)) {
				//If (i,j) is in the range of calculate
				if ( j > i * 2 ){
					kdist[i][j] = INF;
				}
				else{
					if (j == 0){
						kdist[i][j] = kdist[i-1][j] + trans[j+1][j] + kdist[i][j];
					}
					else if (j == 1){

						kdist[i][j] = kmin2(kdist[i-1][j] + trans[j+1][j], kdist[i-1][j-1] + trans[j+1][j+1]) + kdist[i][j];
					}
					else{
						kdist[i][j] = kmin3(kdist[i-1][j] + trans[j+1][j], kdist[i-1][j-1] + trans[j+1][j+1], kdist[i-1][j-2] + trans[j+1][j+2]) + kdist[i][j];
					}
				}
			}
		}
	}
	/*for(i = 0; i < 39; i++){
		input[i].clear();
        sample[i].clear();
	}*/

	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < insize; i++)
		{
			for (j = 0; j < tmpsize; j++)
			{
				disfile << kdist[i][j] << ' ';
			}
			disfile << endl;
		}
	}
}

void kBackTrack(int num){
	//Warning! Only after kmeansdistance() invokes, the function can be called.
	//Will update segdiv[][].
	int insize = kdist.size();
	int tmpsize = kdist[0].size();
	int inp = insize - 1;
	int tmp = tmpsize - 1;
	while(inp >=0){
		segdiv[num][tmp] = inp;
		if (tmp == 0){
			inp--;
		}
		else if (tmp == 1){
			if (kmin2(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1]) == kdist[inp - 1][tmp]){
				inp--;
			}
			else{
				inp--;
				tmp--;
			}
		}
		else{
			if (kmin3(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1], kdist[inp - 1][tmp - 2]) == kdist[inp - 1][tmp]){
				inp--;
			}
			else if (kmin3(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1], kdist[inp - 1][tmp - 2]) == kdist[inp - 1][tmp - 1]){
				inp--;
				tmp--;
			}
			else{
				inp--;
				tmp=tmp-2;
			}
		}
	}
}

void clearforLoop(){
	means.clear();
    covars.clear();
    kdist.clear();
}


void testPrintSample(){
    int i, j;
	//cout << samples.size() << endl; //==2
	//cout << samples[0].size() << endl; //==39
	//cout << samples[1].size() << endl; //==39
	//cout << samples[0][0].size() << endl; //==帧长
	//这样输出来的每一个元素都是第i维的第j帧的元素
    /*for (i = 0; i <39; i++){
		for (j = 0; j < samples[0][0].size(); j++){
			cout << samples[1][i][j] << " ";
	    }
		cout << endl;
	}*/

	/*
	cout << means.size() << endl; //
	cout << means[0].size() << endl;
	for (i = 0; i <39; i++){
		for (j = 0; j < means[0].size(); j++){
			cout << covars[i][j] << " ";
	    }
		cout << endl;
	}*/

	/*for (i=0; i<SEGNUM+1; i++){
		for (j=0; j<SEGNUM; j++)
		{
			cout << trans[i][j] << " ";
		}
		cout << endl;
	}*/

  for(i = 0; i < SAMPLENUM; i++){
	for (j = 0; j <SEGNUM+1; j++){
		cout << segdiv[i][j] << " ";
	}
	cout << endl;
  }
}

void segKMeansProcess(){
	initSeg();
	int i, j;
	for (i = 0; i < 5; i++){
	    initTrans();
        calSegmentMeansCovar();
		for (j = 0; j < SAMPLENUM; j++){
	        calSampDistance(j);
            kmeansDistance();
	        kBackTrack(j);
			kdist.clear();
		}
	    testPrintSample();
		clearforLoop();
	}
}