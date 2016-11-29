#include "singleDTW.h"

vector<vector<double>> input(39);
vector<vector<double>> sample(39);
vector<vector<double>> dist;

int SPEECH = 1;
int SAMP = 2;

void readDCT(const char* filename, int flag){
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
			if (flag == SPEECH){
				input[count].push_back(DCT_value);
			}
			else if (flag == SAMP){
				sample[count].push_back(DCT_value);
			}
			else{
				break;
			}
			count ++;
		}
	}
	// input[i][j]: iMAX = 38 (0-38, 39 Dimensions); jMAX = Number of Frames.
	
}

void calDistance(const char* inputfilename, const char* samplefilename){

	readDCT(inputfilename, SPEECH);
	readDCT(samplefilename, SAMP);
	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = input[1].size();
	int samplesize = sample[1].size();
	dist.resize(inputsize);
	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < samplesize; j++)
		{
			for ( k = 0; k < 39; k++)
			{
				tmpDistance += (input[k][i] - sample[k][j]) * (input[k][i] - sample[k][j]); 
			}
			tmpDistance = sqrt(tmpDistance);
			dist[i].push_back(tmpDistance);
			tmpDistance = 0;
		}
	}
}

double min3(double superdiag, double horizon, double diag){
	double temp = (superdiag < horizon) ? superdiag : horizon;
	return (temp < diag) ? temp : diag;
}

double min2(double horizon, double diag){
	return (horizon < diag) ? horizon : diag;
}

double DTWeditDistance(){
	int insize = dist.size();
	int tmpsize = dist[0].size();
	int i, j;
	for (i = 0; i < insize; i++){
		for(j = 0; j < tmpsize; j++){
			//If i!= 0 and j!= 0;
			if (!(i == 0 && j == 0)) {
				//If (i,j) is in the range of calculate
				if ( j > i * 2 ){
					dist[i][j] = INF;
				}
				else{
					if (j == 0){
						dist[i][j] = dist[i-1][j] + dist[i][j];
					}
					else if (j == 1){

						dist[i][j] = min2(dist[i-1][j], dist[i-1][j-1]) + dist[i][j];
					}
					else{
						dist[i][j] = min3(dist[i-1][j], dist[i-1][j-1], dist[i-1][j-2]) + dist[i][j];
					}
				}
			}
		}
	}
	for(i = 0; i < 39; i++){
		input[i].clear();
        sample[i].clear();
	}
	return dist[insize-1][tmpsize-1];
}

void clearDist(){
	int i;
	for (i = 0; i < dist.size(); i++){
		dist[i].clear();
	}
}

void singleDTW(){
	int i, j;
    string sample[10][5] =
	{"s-0-1.txt","s-0-2.txt","s-0-3.txt","s-0-4.txt","s-0-5.txt",
	"s-1-1.txt","s-1-2.txt","s-1-3.txt","s-1-4.txt","s-1-5.txt",
	"s-2-1.txt","s-2-2.txt","s-2-3.txt","s-2-4.txt","s-2-5.txt",
    "s-3-1.txt","s-3-2.txt","s-3-3.txt","s-3-4.txt","s-3-5.txt",
	"s-4-1.txt","s-4-2.txt","s-4-3.txt","s-4-4.txt","s-4-5.txt",
	"s-5-1.txt","s-5-2.txt","s-5-3.txt","s-5-4.txt","s-5-5.txt",
	"s-6-1.txt","s-6-2.txt","s-6-3.txt","s-6-4.txt","s-6-5.txt",
	"s-7-1.txt","s-7-2.txt","s-7-3.txt","s-7-4.txt","s-7-5.txt",
	"s-8-1.txt","s-8-2.txt","s-8-3.txt","s-8-4.txt","s-8-5.txt",
	"s-9-1.txt","s-9-2.txt","s-9-3.txt","s-9-4.txt","s-9-5.txt"};
	double min[10][5] = {};
	double tmpMin = INF;
	int number = 0;
	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < 5; j++)
		{
			calDistance("input.txt", sample[i][j].c_str());
	        min[i][j] = DTWeditDistance();
	        cout << min[i][j] << " ";
	        clearDist();
			if (tmpMin > min[i][j]) 
			{ 
					tmpMin = min[i][j];
					number = i+1;
			}	
		}
		cout << endl;
	}
	
	cout << "You said: "<< number << endl;
}