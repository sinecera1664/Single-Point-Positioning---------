#include "basic_service.h"

// 辅助代码
// =========================================

/**
* 从buf中提取数据
* 只含有一个静态函数BufToData
*/
template <typename TypeOfBuf>
class BufOutput
{
private:
	/* data */
public:
	BufOutput(/* args */);
	~BufOutput();
	static TypeOfBuf BufToData(unsigned char* buf);
};

template <typename TypeOfBuf>
BufOutput<TypeOfBuf>::BufOutput(/* args */)
{
}

template <typename TypeOfBuf>
BufOutput<TypeOfBuf>::~BufOutput()
{
}

template <typename TypeOfBuf>
TypeOfBuf BufOutput<TypeOfBuf>::BufToData(unsigned char* buf)
{
	TypeOfBuf ret;
	memcpy(&ret, buf, sizeof(TypeOfBuf));
	return ret;
}

/**
* crc32位校验
*/
unsigned int CRC32(const unsigned char* buff, int len)
{
	int i, j;
	unsigned int crc = 0;

	for (i = 0; i < len; i++)
	{
		crc ^= buff[i];
		for (j = 0; j < 8; j++)
		{
			if (crc & 1)
			{
				crc = (crc >> 1) ^ POLYCRC32;
			}
			else
			{
				crc >>= 1;
			}
		}
	}
	return crc;
}


/**
* 提取观测数据
*/
void GetObs(unsigned char* buf, Obs& obs)
{
	obs.tReceiverSurface = { BufOutput<unsigned short>::BufToData(buf + 14),
							0.001 * BufOutput<unsigned int>::BufToData(buf + 16) };

	unsigned char* startPerObs = buf + 28;

	// 解码得到观测值数量
	int obsNum = BufOutput<unsigned int>::BufToData(startPerObs);

	NavSys sys;
	int i, j, n;
	int nSat, satNum = 0;
	int freq, PRN, sigType;
	unsigned int chTrStatus; 

	
	for (i = 0; i < obsNum; i++, startPerObs += 44)
	{
		//只接收GPS，BDS
		PRN = BufOutput<unsigned short>::BufToData(startPerObs + 4);
		chTrStatus = BufOutput<unsigned int>::BufToData(startPerObs + 44);
		n = (chTrStatus >> 16) & 0x07;
		if (n == 0)
		{
			sys = GPS;
		}
		else if (n == 4)
		{
			sys = BDS;
		}
		else
		{
			continue;
		}

		
		sigType = (chTrStatus >> 21) & 0x1F;
		if (sys == GPS)
		{
			if (sigType == 0)
			{
				freq = 0;
			}
			else if (sigType == 9)
			{
				freq = 1;
			}
			else
			{
				freq = 2;
			}
		}
		if (sys == BDS)
		{
			if (sigType == 0 || sigType == 4)
			{
				freq = 0;
			}
			else if (sigType == 2 || sigType == 6)
			{
				freq = 1;
			}
			else
			{
				freq = 2;
			}
		}

		if (freq > 1)
		{
			continue;

		}
			

		// 有效性判断 PRN和isDualBand
		nSat = satNum;
		for (j = 0; j < satNum; j++)
		{
			if (sys == obs.sat[j].sys && PRN == obs.sat[j].PRN)
			{
				nSat = j;
				satNum--;
				obs.sat[nSat].isDualBand = true; 
			}
		}

		obs.sat[nSat].PRN = PRN;
		obs.sat[nSat].sys = sys;

		obs.sat[nSat].psr[freq] = BufOutput<double>::BufToData(startPerObs + 8);
		obs.sat[nSat].psrSigma[freq] = BufOutput<float>::BufToData(startPerObs + 16);
		obs.sat[nSat].adr[freq] = BufOutput<double>::BufToData(startPerObs + 20);
		obs.sat[nSat].adrSigma[freq] = BufOutput<float>::BufToData(startPerObs + 28);
		obs.sat[nSat].dopp[freq] = BufOutput<float>::BufToData(startPerObs + 32);
		obs.sat[nSat].noiseRatio[freq] = BufOutput<float>::BufToData(startPerObs + 36);

		satNum++;
	}

	obs.satObserved = satNum;

	
}

void GetGPSEph(unsigned char* buf, Ephem* gps)
{
	unsigned char* ephemBegin = buf + 28;

	int PRN;
	memcpy(&PRN, ephemBegin, 4);

	if (PRN > 32)
	{
		std::cout << "prn in GetGPSEph() exceeded threshold" << std::endl;
	}
	else
	{
		int i = PRN - 1;
		gps[i].sys = GPS;
		gps[i].PRN = PRN;
		gps[i].isExisted = true;
		memcpy(&gps[i].tow, ephemBegin + 4, 8);
		memcpy(&gps[i].IODE[0], ephemBegin + 16, 4);
		memcpy(&gps[i].IODE[1], ephemBegin + 20, 4);
		memcpy(&gps[i].weekReference, ephemBegin + 24, 4);
		memcpy(&gps[i].zWeek, ephemBegin + 28, 4);
		memcpy(&gps[i].toe, ephemBegin + 32, 8);
		memcpy(&gps[i].A, ephemBegin + 40, 8);
		memcpy(&gps[i].deltaN, ephemBegin + 48, 8);
		memcpy(&gps[i].M0, ephemBegin + 56, 8);
		memcpy(&gps[i].ecc, ephemBegin + 64, 8);
		memcpy(&gps[i].w, ephemBegin + 72, 8);
		memcpy(&gps[i].cuc, ephemBegin + 80, 8);
		memcpy(&gps[i].cus, ephemBegin + 88, 8);
		memcpy(&gps[i].crc, ephemBegin + 96, 8);
		memcpy(&gps[i].crs, ephemBegin + 104, 8);
		memcpy(&gps[i].cic, ephemBegin + 112, 8);
		memcpy(&gps[i].cis, ephemBegin + 120, 8);
		memcpy(&gps[i].I0, ephemBegin + 128, 8);
		memcpy(&gps[i].Idot, ephemBegin + 136, 8);
		memcpy(&gps[i].W0, ephemBegin + 144, 8);
		memcpy(&gps[i].Wdot, ephemBegin + 152, 8);
		memcpy(&gps[i].iodc, ephemBegin + 160, 4);
		memcpy(&gps[i].toc, ephemBegin + 164, 8);
		memcpy(&gps[i].tgd[0], ephemBegin + 172, 8);
		memcpy(&gps[i].af0, ephemBegin + 180, 8);
		memcpy(&gps[i].af1, ephemBegin + 188, 8);
		memcpy(&gps[i].af2, ephemBegin + 196, 8);
		memcpy(&gps[i].URA, ephemBegin + 216, 8);
	}
}

void GetBDSEph(unsigned char* buf, Ephem* bds)
{
	unsigned char* ephemBegin = buf + 28;

	int PRN, toc, toe;
	memcpy(&PRN, ephemBegin, 4);

	double sqrtA;

	if (PRN > 63)
	{
		std::cout << "prn in GetBDSEph() exceeded threshold" << std::endl;
	}
	else
	{
		int i = PRN - 1;
		bds[i].sys = BDS;
		bds[i].PRN = PRN;
		bds[i].isExisted = true;
		memcpy(&bds[i].weekReference, ephemBegin + 4, 4);
		bds[i].weekReference += 1356;
		memcpy(&bds[i].URA, ephemBegin + 8, 8);
		memcpy(&bds[i].tgd[0], ephemBegin + 20, 8);
		memcpy(&bds[i].tgd[1], ephemBegin + 28, 8);
		memcpy(&toc, ephemBegin + 40, 4);
		bds[i].toc = toc + 14.0;
		memcpy(&bds[i].af0, ephemBegin + 44, 8);
		memcpy(&bds[i].af1, ephemBegin + 52, 8);
		memcpy(&bds[i].af2, ephemBegin + 60, 8);
		memcpy(&toe, ephemBegin + 72, 4);
		bds[i].toe = toe + 14.0;
		memcpy(&sqrtA, ephemBegin + 76, 8);
		bds[i].A = sqrtA * sqrtA;
		memcpy(&bds[i].ecc, ephemBegin + 84, 8);
		memcpy(&bds[i].w, ephemBegin + 92, 8);
		memcpy(&bds[i].deltaN, ephemBegin + 100, 8);
		memcpy(&bds[i].M0, ephemBegin + 108, 8);
		memcpy(&bds[i].W0, ephemBegin + 116, 8);
		memcpy(&bds[i].Wdot, ephemBegin + 124, 8);
		memcpy(&bds[i].I0, ephemBegin + 132, 8);
		memcpy(&bds[i].Idot, ephemBegin + 140, 8);
		memcpy(&bds[i].cuc, ephemBegin + 148, 8);
		memcpy(&bds[i].cus, ephemBegin + 156, 8);
		memcpy(&bds[i].crc, ephemBegin + 164, 8);
		memcpy(&bds[i].crs, ephemBegin + 172, 8);
		memcpy(&bds[i].cic, ephemBegin + 180, 8);
		memcpy(&bds[i].cis, ephemBegin + 188, 8);
	}
}
// 辅助代码结束
// =========================================

double Distance(double* left, double* right)
{
	return sqrt((left[0] - right[0]) * (left[0] - right[0]) + (left[1] - right[1]) * (left[1] - right[1]) + (left[2] - right[2]) * (left[2] - right[2]));
}

BLH XYZToBLH(double* var, NavSys sys)
{
	if (var[0] == 0 && var[1] == 0 && var[2] == 0)
	{
		return BLH{ 0, 0, -6.371e6 };
	}

	double A;
	double E2;

	if (sys == GPS)
	{
		A = WGS84_A;
		E2 = WGS84_E2;
	}
	else
	{
		A = CGCS2000_A;
		E2 = CGCS2000_E2;
	}

	double p = sqrt(var[0] * var[0] + var[1] * var[1]);

	double diff = 1.0;
	double tanB = 0;
	double cosB = 1 / sqrt(1 + tanB * tanB);
	double sinB = tanB * cosB;
	double N = A / sqrt(1 - E2 * sinB * sinB);

	//最多进行30次迭代
	for (int i = 0; i < 30; i++)
	{
		if (diff > 1e-10)
		{
			double NewTanB = (var[2] + E2 * N * sinB) / p;
			diff = fabs(NewTanB - tanB);
			tanB = NewTanB;
			cosB = 1 / sqrt(1 + tanB * tanB);
			sinB = tanB * cosB;
			N = A / sqrt(1 - E2 * sinB * sinB);
		}
		else
		{
			break;
		}
		//判断迭代是否完成
		if (i == 30 && diff > 1e-10)
		{
			std::cout<<"Iteration incomplete!\n"<<std::endl;
			exit(EXIT_FAILURE);
		}
	}

	double B = atan(tanB);
	double H = p / cosB - N;
	double L = atan2(var[1], var[0]);  

	return BLH(B, L, H);
}

void XYZToENU(double* oriXYZ, double* tarXYZ, double* oriBLH, double* tarENU)
{
	const double sinB = sin(oriBLH[0]);
	const double cosB = cos(oriBLH[0]);
	const double sinL = sin(oriBLH[1]);
	const double cosL = cos(oriBLH[1]);

	vector<vector<double> > s = { {(-1) * sinL, cosL, 0},
								  {(-1) * sinB * cosL, (-1) * sinB * sinL, cosB},
								  {cosB * cosL, cosB * sinL, sinB } };
	vector<vector<double> > ri = { {tarXYZ[0] - oriXYZ[0]}, {tarXYZ[1] - oriXYZ[1]} , {tarXYZ[2] - oriXYZ[2]} };
	Matrixd S(s);
	Matrixd Ri(ri);
	Matrixd enu = S * Ri;
	tarENU[0] = enu[0][0];
	tarENU[1] = enu[1][0];
	tarENU[2] = enu[2][0];
}

/**
* 读取buffer并解码
* @param buf 输入buffer
* @param lengthTail 上一次读取的尾巴的长度
* @param obs 输出观测数据
* @param gps 输出GPS星历
* @param bds 输出BDS星历
*/
bool DecodeBuffer(unsigned char* buf, int& lengthTail, Obs& obs, Ephem* gps, Ephem* bds)
{
	unsigned char* p = buf;

	bool isObserved = false;

	while (1)
	{

		//是否一条完整语句，否则下次读取

		while (lengthTail > 32)
		{
			if (buf[0] == 0xAA && buf[1] == 0x44 && buf[2] == 0x12)
			{
				break;
			}
			else
			{
				p++;
				lengthTail--;
			}
		}

		//p指向AA

		if (lengthTail <= 32)
		{
			memcpy(buf, p, lengthTail);
			return isObserved;
		}
		 
		unsigned short mesLength = BufOutput<unsigned short>::BufToData(p + 8);
		
		if (lengthTail < (32 + mesLength))
		{
			memcpy(buf, p, lengthTail); 
			return isObserved;
		}

		      

		// CRC检验
		if (CRC32(p, 28 + mesLength) != BufOutput<unsigned int>::BufToData(p + 28 + mesLength))
		{
			lengthTail -= (32 + mesLength);
			p += (32 + mesLength);
			continue;
		}

		           

		// 解码
		unsigned short mesID = BufOutput<unsigned short>::BufToData(p + 4);

		switch (mesID)
		{
		case 43:
			isObserved = true;
			GetObs(p, obs);
			break;
		case 7:
			GetGPSEph(p, gps);
			break;
		case 1696:
			GetBDSEph(p, bds);
			break;
		default:
			break;
		}
		
		lengthTail = lengthTail - (32 + mesLength);
		p = p + (32 + mesLength);
	}

	std::cout << "error in  DecodeBuffer!!! "<<endl;
	exit(EXIT_FAILURE);
}