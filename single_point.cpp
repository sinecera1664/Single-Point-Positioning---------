#include "single_point.h"
// 辅助代码
// =========================================
/**
* 星历有效性检测
* tSatSurface和toe
*/
bool IsSatVaild(const GPSTime& tSatSurface, const Ephem& ehm)
{
	if (!ehm.isExisted)
	{
		return false;
	}

	double tk = (double(tSatSurface.week) - ehm.weekReference) * 604800.0 + tSatSurface.sow - ehm.toe;

	if (fabs(tk) > (ehm.sys == GPS ? 7200 : 3600))
	{
		return false;
	}
	return true;
}

/**
* 卫星位置、速度、钟差、钟速
* @param tSatSurface 卫星表面时间
* @param eph 卫星星历
* @param res  卫星定位结果，包含位置、速度、钟差、钟速
* @return bool 定位是否成功
*/
bool SatPos(GPSTime& tSatSurface, Ephem& eph, SatMedian& res)
{

	res.isValid = false;
	if (!IsSatVaild(tSatSurface, eph))
	{
		return false;
	}

	res.sys = eph.sys;

	double GM;
	double WE;
	if (eph.sys == GPS)
	{
		GM = GM_GPS;
		WE = WE_GPS;
	}
	else
	{
		GM = GM_BDS;
		WE = WE_BDS;
	}

	//轨道长度
	double A = eph.A;

	//平均运动角速度
	double n0 = sqrt(GM / (A * A * A));

	//相对于星历参考历元的时间
	double tk = (double(tSatSurface.week) - eph.weekReference) * 604800.0 + tSatSurface.sow - eph.toe;
	double tk2 = (double(tSatSurface.week) - eph.weekReference) * 604800.0 + tSatSurface.sow - eph.toc;
	res.satClkD = eph.af0 + eph.af1 * tk2 + eph.af2 * tk2 * tk2;
	tk = tk - res.satClkD;

	//平均运动角速度改正
	double n = n0 + eph.deltaN;

	//平近角点
	double Mk = eph.M0 + n * tk;

	//迭代偏近角点
	double e = eph.ecc;
	double Ek = Mk;
	double Etmp;
	double diff = 1;
	//最多进行30次迭代
	for (int i = 0; i < 30; i++)
	{
		
		if (diff > 1e-10)//对卫星位置的影响 3e7 * 1e-10
		{
			Etmp = Ek;
			Ek = Mk + e * sin(Etmp);
		}

		else
		{
			break;
		}
		//判断迭代是否完成
		if (i == 30 && diff > 1e-10)
		{
			printf("Iteration incomplete!\n");
			exit(EXIT_FAILURE);
		}
		diff = fabs(Ek - Etmp);
	}
	//真近点角
	double Vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
	//升交角距
	double faiK = Vk + eph.w;
	//二阶调和改正
	double deltaU = eph.cus * sin(2 * faiK) + eph.cuc * cos(2 * faiK);
	double deltaR = eph.crs * sin(2 * faiK) + eph.crc * cos(2 * faiK);
	double deltaI = eph.cis * sin(2 * faiK) + eph.cic * cos(2 * faiK);
	//改正的升交角距、向径、轨道倾角
	double Uk = faiK + deltaU;
	double Rk = A * (1 - e * cos(Ek)) + deltaR;
	double Ik = eph.I0 + eph.Idot * tk + deltaI;
	//轨道平面上的位置
	double xk = Rk * cos(Uk);
	double yk = Rk * sin(Uk);
	double Wk; //OmegaK

	bool isGeo = ((res.sys == BDS) && (eph.PRN < 6 || 58 < eph.PRN));

	// 对Geo卫星特殊处理
	Wk = eph.W0 + (eph.Wdot - (1.0 - isGeo) * WE) * tk - WE * (res.sys == GPS ? eph.toe : eph.toe - 14.0);

	//地固坐标系下的位置
	res.satP[0] = xk * cos(Wk) - yk * cos(Ik) * sin(Wk);
	res.satP[1] = xk * sin(Wk) + yk * cos(Ik) * cos(Wk);
	res.satP[2] = yk * sin(Ik);

	//相对论改正
	res.satClkD += F * e * sqrt(A) * sin(Ek);

	res.tgd = eph.tgd[0];

	double Ekdot = n / (1 - e * cos(Ek));
	double Faikdot = sqrt((1 + e) / (1 - e)) * pow(cos(Vk / 2) / cos(Ek / 2), 2) * Ekdot;
	double Ukdot = 2 * (eph.cus * cos(2 * faiK) - eph.cuc * sin(2 * faiK)) * Faikdot + Faikdot;
	double Rkdot = A * e * sin(Ek) * Ekdot + 2 * (eph.crs * cos(2 * faiK) - eph.crc * sin(2 * faiK)) * Faikdot;
	double Ikdot = eph.Idot + 2 * (eph.cis * cos(2 * faiK) - eph.cic * sin(2 * faiK)) * Faikdot;

	//Geo卫星特殊处理
	double Wkdot = eph.Wdot - (1.0 - isGeo) * WE;

	Matrixd R({ {cos(Wk), (-1) * sin(Wk) * cos(Ik), (-1) * (xk * sin(Wk) + yk * cos(Wk) * cos(Ik)), yk * sin(Wk) * sin(Ik)},
								{sin(Wk), cos(Wk) * cos(Ik), (xk * cos(Wk) - yk * sin(Wk) * cos(Ik)), yk * cos(Wk) * sin(Ik)},
								{ 0, sin(Ik), 0, yk * cos(Ik)} });
	double xkdot = Rkdot * cos(Uk) - Rk * Ukdot * sin(Uk);
	double ykdot = Rkdot * sin(Uk) + Rk * Ukdot * cos(Uk);


	Matrixd R_v({ {xkdot}, {ykdot}, {Wkdot}, {Ikdot} });
	Matrixd satV = R * R_v;
	satV.OutToArray(3, res.satV);

	//Geo卫星特殊处理
	if (isGeo)
	{
		const double xRotate = -5 * PI / 180;
		const double sinXR = sin(xRotate);
		const double cosXR = cos(xRotate);
		const double zRotate = WE * tk;
		const double sinZR = sin(zRotate);
		const double cosZR = cos(zRotate);

		Matrixd RvBDS({ {res.satV[0]}, {res.satV[1]}, {res.satV[2]}, {WE} });
		Matrixd LvBDS({ {cosZR, sinZR * cosXR, sinZR * sinXR, -1 * sinZR * res.satP[0] + cosZR * cosXR * res.satP[1] + cosZR * sinXR * res.satP[2]},
						  { -1 * sinZR, cosZR * cosXR, cosZR * sinXR, -1 * cosZR * res.satP[0] - sinZR * cosXR * res.satP[1] - sinZR * sinXR * res.satP[2]},
						  {0, -1 * sinXR, cosXR, 0 } });
		Matrixd temp1 = LvBDS * RvBDS;
		temp1.OutToArray(3, res.satV);


		Matrixd Rx({ {1, 0, 0},
					 {0, cosXR, sinXR},
					 {0, -1 * sinXR, cosXR} });
		Matrixd Rz({ {cosZR, sinZR, 0},
					 {-1 * sinZR, cosZR,0},
					 {0, 0, 1 } });
		Matrixd pos({ {res.satP[0]},{res.satP[1]} ,{res.satP[2]} });
		Matrixd temp2 = Rz * Rx * pos;
		temp2.OutToArray(3, res.satP);

	}

	//钟速
	res.satClkV = eph.af1 + 2 * eph.af2 * tk2 + F * e * sqrt(A) * cos(Ek) * Ekdot;

	res.isValid = true;
	return true;
}

/**
* 对流层误差
* 超出40km，超出有效距离，舍去并返回0
* @param H 测站高度
* @param E 卫星高度角
* @return double 对流层改正数
*/
double Hopefield(double H, double E)
{
	if (fabs(H) > 40000)
	{
		return 0;
	}

	Weather stdM = { 288.16, 1013.25, 0.5 };
	Weather M = { stdM.t - 0.0065 * H, stdM.p * pow(1 - 0.0000226 * H, 5.225), stdM.RH * exp(-0.0006396 * H) };
	double hd = 40136 + 148.72 * (stdM.t - 273.16);
	double Kd = 155.2e-7 * M.p * (hd - H) / M.t;
	double e = M.RH * exp(-37.2465 + 0.213166 * M.t - 0.000256908 * M.t * M.t);
	double Kw = 155.2e-7 * 4810 * e * (11000 - H) / (M.t * M.t);
	return (Kd / sin(PI * sqrt(E * E + 6.25) / 180)) + (Kw / sin(PI * sqrt(E * E + 2.25) / 180));

	return 0;
}

/**
* 地球自转,坐标旋转
*/
void EarthRotate(double deltaT, double* posSatOri, NavSys& sys, double* posSatRes)
{
	double We = sys == GPS ? WE_GPS : WE_BDS;
	double angleRotate = We * deltaT;

	vector<vector<double> > le = { {cos(angleRotate), sin(angleRotate), 0},
								   {-sin(angleRotate), cos(angleRotate), 0},
								   {0, 0, 1 } };
	Matrixd Le(le);
	Matrixd pos({ {posSatOri[0]},{posSatOri[1]},{posSatOri[2] } });
	Matrixd temp = Le * pos;
	temp.OutToArray(3, posSatRes);
}


/**
* IF组合伪距
* @param psr 伪距
* @param sys 卫星系统
* @return double IF组合伪距
*
*/
double IF(double* psr, NavSys& sys)
{
    double freq1 = sys == GPS ? L1 : B1;
    double freq2 = sys == GPS ? L2 : B3;
    return (freq1 * freq1 * psr[0] - freq2 * freq2 * psr[1]) / (freq1 * freq1 - freq2 * freq2);
}
// 辅助代码结束
// =========================================



/**
* 调用SatPos，生成IF组合改正所有卫星定位结果，用于SPP的输入
* @param obs 输入观测数据，初步处理，并输出
* @param gps GPS星历
* @param bds BDS星历
*/
void InitRes(Obs& obs, Ephem* gps, Ephem* bds)
{
	double psrIF;
	GPSTime tSatSurface;
	obs.validNum = 0;

	for (int i = 0; i < obs.satObserved; i++)
	{
		//单频不参与计算
		if (!obs.sat[i].isDualBand)
		{
			obs.validSat[i].isValid = false;   
			continue;                 //跳到下一次循环
		}

		NavSys sys = obs.sat[i].sys;
		psrIF = IF(obs.sat[i].psr, sys);
		tSatSurface = { obs.tReceiverSurface.week, obs.tReceiverSurface.sow - (psrIF / LIGHT_SPEED) };

		if (SatPos(tSatSurface, sys == GPS ? gps[obs.sat[i].PRN - 1] : bds[obs.sat[i].PRN - 1], obs.validSat[i]))
		{
			obs.validNum++;
			obs.validSat[i].pseu = psrIF + (sys == GPS ? 0 : (LIGHT_SPEED * B1 * B1 * obs.validSat[i].tgd / (B3 * B3 - B1 * B1)));
		}
	}
}


/**
* SPP定位
* @param  obs 观测数据
* @param  resSPP 单点定位结果
* @return bool 定位是否成功
*/
bool GetSPP(Obs& obs, Orientation& resSPP)
{
	if (obs.validNum < 5)
	{
		return false;
	}

	NavSys sys;
	BLH orientationBLH;

	double W[MAX_SAT_NUM], B[5 * MAX_SAT_NUM], satPosAfterErotate[3];
	
	for (auto& var : W)
	{
		var = 0.0;
	}
	for (auto& var : B)
	{
		var = 0.0;
	}
	double deltaT, elev, trop, rho, diff = 1.0;

	//判断卫星数目有效性
	bool isSatNumValid = true;  	
	bool isSingleSys;
	Matrixd matQ, matW, matB, matX;

	//最多进行10次迭代
	for (size_t j = 0; j < 10; j++)
	{
		while (diff > 1e-8 )
		{			                            
			int m = 0;//与i不同步
			int gpsN = 0, bdsN = 0;

			for (int i = 0; i < obs.satObserved; i++)
			{
				if (!obs.validSat[i].isValid)
				{
					continue;
				}

				sys = obs.validSat[i].sys;
				//信号传播时间
				deltaT = (obs.validSat[i].pseu / LIGHT_SPEED) + obs.validSat[i].satClkD - (sys == GPS ? resSPP.pos[3] / LIGHT_SPEED : resSPP.pos[4] / LIGHT_SPEED);
				//坐标转换
				EarthRotate(deltaT, obs.validSat[i].satP, sys, satPosAfterErotate);
				orientationBLH = XYZToBLH(resSPP.pos, sys);

				//计算elev，高度角
				if (fabs(orientationBLH.H) < 40000)
				{
					double tarENU[3];
					double blh[3] = { orientationBLH.B,orientationBLH.L,orientationBLH.H };
					XYZToENU(resSPP.pos, obs.validSat[i].satP, blh, tarENU);
					elev = asin(tarENU[2] / sqrt(tarENU[0] * tarENU[0] + tarENU[1] * tarENU[1] + tarENU[2] * tarENU[2])) * 180 / PI;
				}
				else
				{
					elev = 90.0;
				}

				//有效性判定，高度角小于20°舍去
				if (elev < 20)
				{
					obs.validSat[i].isValid = false;
					--obs.validNum;
					if (obs.validNum < 5)
					{
						break;
					}
					continue;  
				}
				//对流层延迟
				trop = Hopefield(orientationBLH.H, elev);

				//根据用户初值得到的几何距离
				rho = Distance(resSPP.pos, satPosAfterErotate);

				//用于生成 W 矩阵            
				W[m] = obs.validSat[i].pseu - (rho + resSPP.pos[3 + sys - 1] - LIGHT_SPEED * obs.validSat[i].satClkD + trop);

				for (int j = 0; j < 3; j++)
				{
					B[m * 5 + j] = (resSPP.pos[j] - satPosAfterErotate[j]) / rho;
				}

				if (sys==GPS)
				{
					B[m * 5 + 3] = 1; 
					B[m * 5 + 4] = 0;

				}
				else
				{
					B[m * 5 + 3] = 0; 
					B[m * 5 + 4] = 1;
				}				
				sys == GPS ? gpsN++ : bdsN++;  
				m++;
			}
			if (obs.validNum < 5)
			{
				isSatNumValid = false;				
				std::cout << "satNum is not valid" << std::endl;
				break;
			}

			
			//开始最小二乘求解，矩阵行数由nValidNum决定
			matW = Matrixd(obs.validNum, 1);			
			for (size_t i = 0; i < obs.validNum; i++)
			{
				matW[i][0] = W[i];
			}

			//单系统判定,B矩阵单系统列数为4
			bool isGPSExisted = (gpsN != 0);
			bool isBDSExisted = (bdsN != 0);
			isSingleSys = !isGPSExisted || !isBDSExisted;
			
			if (isSingleSys)
			{
				matB = Matrixd(obs.validNum, 4);

				for (size_t i = 0; i < obs.validNum; i++)
				{
					matB[i][0] = B[i * 5 + 0];
					matB[i][1] = B[i * 5 + 1];
					matB[i][2] = B[i * 5 + 2];
					matB[i][3] = 1.0;
				}
			}
			else
			{
				matB = Matrixd(obs.validNum, 5);
				for (size_t i = 0; i < obs.validNum; i++)
				{					
					for (size_t j = 0; j < 5; j++)
					{
						matB[i][j] = B[i * 5 + j];
					}
				}

			}	
			matQ = (matB.Transposition() * matB).Inversion();  
			matX = matQ * matB.Transposition() * matW;
			
			//位置、钟差修正
			for (int k = 0; k < 3; k++)
			{				
				resSPP.pos[k] += matX[k][0];
			}

			if (isSingleSys)
			{
				resSPP.pos[3] += matX[3][0];
			}
			else
			{
				resSPP.pos[3] += matX[3][0];
				resSPP.pos[4] += matX[4][0];
			}

			//和上一次站点的距离
			diff = sqrt(matX[0][0] * matX[0][0] + matX[1][0] * matX[1][0] + matX[2][0] * matX[2][0]);

			//迭代未完成
			if (j == 10 && diff > 1e-8)
			{
				return false;
			}			
		}
	}
	//精度评定
	if (isSatNumValid)
	{
		resSPP.tObserve = obs.tReceiverSurface;
		Matrixd matV = matB * matX - matW;
		resSPP.sigmaPos = sqrt(matV[0][0]*matV[0][0]+ matV[1][0] * matV[1][0]+ matV[2][0] * matV[2][0]+ matV[3][0] * matV[3][0])
			/ sqrt(obs.validNum - (isSingleSys ? 4 : 5));
		resSPP.PDOP = sqrt(matQ[0][0] + matQ[1][1] + matQ[2][2]);
		return true;
	}
	return false;
}

/**
* 单点测速
* @param  obs 观测数据
* @param  resSPV 单点测速结果
*/
void GetSPV(Obs& obs, Orientation& resSPV)
{
	NavSys sys;
	double W[MAX_SAT_NUM], B[4 * MAX_SAT_NUM], satPosAfterErotate[3], satVelAfterErotate[3];
	double lambda, deltaT, rho;
	double rhoDerivative = 0.0;
	for (auto &var:W)
	{
		var = 0.0;
	}
	for (auto & var : B)
	{
		var = 0.0;
	}
	
	int m = 0;         //用于计数
	for (int i = 0; i < obs.satObserved; i++)
	{
		if (!obs.validSat[i].isValid)
		{
			continue; //下一颗卫星
		}

		sys = obs.validSat[i].sys;

		deltaT = (obs.validSat[i].pseu / LIGHT_SPEED) + obs.validSat[i].satClkD -
			     (sys == GPS ? resSPV.pos[3] / LIGHT_SPEED : resSPV.pos[4] / LIGHT_SPEED);
		
		//地球旋转
		EarthRotate(deltaT, obs.validSat[i].satP, sys, satPosAfterErotate);
		EarthRotate(deltaT, obs.validSat[i].satV, sys, satVelAfterErotate);

		rho = Distance(resSPV.pos, satPosAfterErotate);

		for (int j = 0; j < 3; ++j)
		{
			rhoDerivative += (satPosAfterErotate[j] - resSPV.pos[j]) * satVelAfterErotate[j];
		}
		rhoDerivative = rhoDerivative / rho;
		

		if (sys == GPS)
		{
			lambda = LIGHT_SPEED / (L1 * 1e6);
		}
		else
		{
			lambda = LIGHT_SPEED / (B1 * 1e6);
		}
		//dopp单位转换为m/s
		W[m] = -1 * obs.sat[i].dopp[0] * lambda - (rhoDerivative - LIGHT_SPEED * obs.validSat[i].satClkV);  

		for (int j = 0; j < 3; j++)
		{
			B[m * 4 + j] = (resSPV.pos[j] - satPosAfterErotate[j]) / rho;
		}
		B[m * 4 + 3] = 1;
		m++;
	}

	//最小二乘求解
	Matrixd matW(obs.validNum, 1);
	Matrixd matB(obs.validNum, 4);
	for (size_t i = 0; i < obs.validNum; i++)
	{
		matW[i][0] = W[i];
		for (size_t j = 0; j < 4; j++)
		{
			matB[i][j] = B[i * 4 + j];
		}
	}

	Matrixd matX = (matB.Transposition() * matB).Inversion() * matB.Transposition() * matW;
	matX.OutToArray(4, resSPV.vel);

	Matrixd matV = matB * matX - matW;
	double squareSum = 0;
	for (size_t i = 0; i < obs.validNum; i++)
	{
		squareSum += matV[i][0];
	}
	resSPV.sigmaVel = sqrt(squareSum) / sqrt(obs.validNum - 4);
}