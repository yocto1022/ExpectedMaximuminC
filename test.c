#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "matrix.h"

#define NumofDist 2
#define NumofSample 500
#define thretholdofConvergenve 0.0001

//分布パラメータ
#define StandardDist1 1.0	//分布1の標準偏差
#define StarndrdDist2 2.0	//分布2の標準偏差
#define MeanDist1 0.0		//分布1の期待値
#define MeanDist2 10.0		//分布2の期待値
#define Ratio 0.9		//分布1の割合


void optimize(double sample[NumofSample]);//パラメータ最適化
double calc_gaussprobability(double in,double mean, double variance);//ガウス分布でのサンプルの生起pdf
double gauss_distribution(double devitation);//いつもの


//最適化したいパラメータ
double theta[NumofDist][2];//母数パラメータ(分布数×パラメータ数(正規分布だと平均と分散))
double mixingRatio[NumofDist];//混合比率

int main(int argc, char const * argv[]){
	int sampleIndex;
	double sampleData[NumofSample];

	for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
		if(sampleIndex<(NumofSample*Ratio)){//前半部分に標準の
			sampleData[sampleIndex] = gauss_distribution(1.0);
		}else{
			sampleData[sampleIndex] = gauss_distribution(2.0) + 10.0;
		}
	}
	optimize(sampleData);

	printf("Dist1:%f/%f,Dist2:%f/%f\n",theta[0][0],theta[0][1],theta[1][0],theta[1][1]);
	/*
	for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
		printf("Sample%d:%f\n",sampleIndex,sampleData[sampleIndex]);
	}
	*/

	return 0;
}

void optimize(double sample[NumofSample]){
	int sampleIndex,distIndex;
	double burdenRate[NumofSample][NumofDist];//負担率
	double sumofburdenRate[NumofDist];//負担率の合計値
	//まず適当に初期値を設定(平均は適当，分散1,混合比率は均等)
	for(distIndex=0;distIndex<NumofDist;distIndex++){
		theta[distIndex][0] = gauss_distribution(1.0);
		theta[distIndex][1] = 1.0;
	}
	for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
		for(distIndex=0;distIndex<NumofDist;distIndex++){
			mixingRatio[distIndex] = 1.0/NumofDist;
		}
	}
	printf("Mixing:%f/%f\tDist1:%f/%f\tDist2:%f/%f\n",mixingRatio[0],mixingRatio[1],theta[0][0],theta[0][1],theta[1][0],theta[1][1]);

	//負担率周りの初期化
	for(distIndex=0;distIndex<NumofDist;distIndex++){
		sumofburdenRate[distIndex] = mixingRatio[distIndex] * NumofSample;
	}


	//ここからいい感じになるまでループ
	double temp[NumofDist],temp2;//負担率推定の際のtemp
	double temp3;//パラメータ推定の際のtemp
	double temp4;//収束判定用
	double logLikelihood,logLikelihoodinPreviousEpoch=INFINITY;//対数尤度(1個前の)
	int itr = 0;
	for(;;){
		//負担率の推定
		for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
			temp2 = 0.0;
			for(distIndex=0;distIndex<NumofDist;distIndex++){//まず分母を求めたい
				temp[distIndex] = mixingRatio[distIndex] * calc_gaussprobability(sample[sampleIndex],theta[distIndex][0],theta[distIndex][1]);
				temp2 += temp[distIndex];
			}
			for(distIndex=0;distIndex<NumofDist;distIndex++){//分母を求めれたので負担率を実際に計算
				burdenRate[sampleIndex][distIndex] = temp[distIndex]/temp2;
			}
		}

		//各パラメータの最尤推定
		//各分布の平均値
		for(distIndex=0;distIndex<NumofDist;distIndex++){
			temp3 = 0.0;
			for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
				temp3 += (burdenRate[sampleIndex][distIndex]*sample[sampleIndex]);
			}
			theta[distIndex][0] = temp3 / sumofburdenRate[distIndex];
		}
		//分散パラメータ
		for(distIndex=0;distIndex<NumofDist;distIndex++){
			temp3 = 0.0;
			for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
				temp3 += (burdenRate[sampleIndex][distIndex]*pow(sample[sampleIndex]-theta[distIndex][0],2.0));
			}
			theta[distIndex][1] = temp3 / sumofburdenRate[distIndex];
		}
		//最初にN_K(負担率の合計値)とmixingrate
		for(distIndex=0;distIndex<NumofDist;distIndex++){
			sumofburdenRate[distIndex] = 0.0;//初期化
			for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
				sumofburdenRate[distIndex] += burdenRate[sampleIndex][distIndex];
			}
			mixingRatio[distIndex] = sumofburdenRate[distIndex]/NumofSample;
		}


		//収束判定
		logLikelihood = 0.0;
		for(sampleIndex=0;sampleIndex<NumofSample;sampleIndex++){
			temp4 = 0.0;
			for(distIndex=0;distIndex<NumofDist;distIndex++){
				temp4 += (mixingRatio[distIndex]*calc_gaussprobability(sample[sampleIndex],theta[distIndex][0],theta[distIndex][1]));
			}
			logLikelihood += log(temp4);
		}
		//推定パラメータ
		printf("%d\tMixing:%f/%f\tDist1:%f/%f\tDist2:%f/%f\tlike:%f\n",itr,mixingRatio[0],mixingRatio[1],theta[0][0],theta[0][1],theta[1][0],theta[1][1],logLikelihood);
		if(fabs(logLikelihood-logLikelihoodinPreviousEpoch)<thretholdofConvergenve){//収束判定
			printf("%f\t%f\t%f\n",logLikelihood,logLikelihoodinPreviousEpoch,logLikelihood-logLikelihoodinPreviousEpoch);
			break;
		}else{
			logLikelihoodinPreviousEpoch = logLikelihood;
			itr++;
		}

	}
	
	return;

}

double calc_gaussprobability(double in,double mean, double variance){
	double out = exp(-pow(in-mean,2.0)/(2.0*variance))/sqrt(2.0*M_PI*variance);
	return out;
}

double gauss_distribution(double devitation){
	return sqrt(-2.0 * pow(devitation,2) * log(genrand_real3())) * cos(2.0 * M_PI * genrand_real3());
}
