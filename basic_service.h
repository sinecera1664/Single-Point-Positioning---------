/**
* 此文件包含所有数据类型定义，宏定义
* 包含方法：坐标转换，文件读取和解码
* 有效性：是否相应星历、星历超时、不是双频、高度角过低
*/

#pragma once

#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define MAX_SAT_NUM 36

#define MAXRAWLEN 10240
#define POLYCRC32 0xEDB88320u //CRC32 polynomial 
#define GPSSATNUM 32
#define BDSSATNUM 63
#define MAXRAWLEN_SOCKET 30000

#define LIGHT_SPEED 2.99792458e8
#define F -4.442807633e-10
#define PI 3.1415926535898

//单位GM
#define GM_GPS 3.986005e14
#define GM_BDS 3.986004418e14

//自转角速度
#define WE_GPS 7.2921151467e-5
#define WE_BDS 7.2921150e-5     

//椭球参数
#define WGS84_A  6378137.0
#define WGS84_f  1/298.257223565
#define WGS84_E2  WGS84_f*(2-WGS84_f)

#define CGCS2000_A 6378137.0
#define CGCS2000_f 1/298.257223563
#define CGCS2000_E2  CGCS2000_f*(2-CGCS2000_f)

#define L1 1575.42
#define L2 1227.6

#define B1 1561.098
#define B3 1268.52

using namespace mym;
using std::vector;
using std::cout;
using std::endl;

struct GPSTime
{
	int week;
	double sow;

	GPSTime() = default;
	GPSTime(int _week, double _sow)
	{
		week = _week;
		sow = _sow;
	}
	
};

struct BLH
{
	//单位rad
	double B;  
	double L;
	double H;

	BLH() = default;

	BLH(double _B, double _L, double _H)
	{
		B = _B;
		L = _L;
		H = _H;
	}
};

enum NavSys
{
	GPS = 1,
	BDS,
	GAL,
	GLO
};


struct Sat
{
	int PRN;          //卫星编号
	bool isDualBand;  //双频判断
	NavSys sys;       //卫星系统
	double psr[2];    //双频伪距观测值
	double adr[2];    //双频相位观测值
	float dopp[2];    //多普勒频移，
	float noiseRatio[2];//载噪比
	float psrSigma[2];//伪距精度
	float adrSigma[2];//载波相位精度
};

struct SatMedian
{
	bool isValid;      //有效性判断
	NavSys sys;        //卫星系统
	double satP[3];  
	double satV[3];
	double satClkD;    //钟差
	double satClkV;    //钟速
	double tgd;        //群延迟
	double pseu;       //IF组合伪距
};

struct Obs
{
	GPSTime tReceiverSurface;
	int satObserved;
	int validNum;

	Sat sat[36];
	SatMedian validSat[36];

	Obs()
	{
		memset(&tReceiverSurface, 0, sizeof(tReceiverSurface));
		satObserved = 0;
		validNum = 0;
		memset(&sat, 0, sizeof(sat));
		memset(&validSat, 0, sizeof(validSat));
	}
	
};

struct Ephem
{
	bool isExisted;
	NavSys sys;
	int PRN, iodc;
	int weekReference;   //toe周
	double tow;
	double toe;          //星历参考时间
	double A;            //长半轴
	double deltaN;       //卫星平均运动速率与计算值之差
	double M0;           //参考时间的平近角点
	double ecc;          //偏心率（无量纲）
	double w;            //近地点幅角
	double cuc, cus, crc, crs, cic, cis;//改正项的振幅
	double Idot;         //轨道倾角变化率
	double I0, W0, Wdot;
	double toc, tgd[2];
	double af0, af1, af2, URA;
	int IODE[2], zWeek;
	Ephem() = default;
};

struct Weather
{
	double t;  //绝对温度
	double p;  //气压
	double RH; //相对湿度
};

struct Orientation
{
	GPSTime tObserve;
	double pos[5]; //包含钟差
	double vel[4]; //包含钟速
	double sigmaPos;  
	double sigmaVel;
	double PDOP;

	Orientation()
	{
		memset((char*)this, 0, sizeof(tObserve)+12 * sizeof(double));
	}
	
};

double Distance(double* left, double* right);

BLH XYZToBLH(double* var, NavSys sys);

void XYZToENU(double* oriXYZ, double* tarXYZ, double* oriBLH, double* tarENU);

bool DecodeBuffer(unsigned char* buf, int& lengthTail, Obs& obs, Ephem* gps, Ephem* bds);
