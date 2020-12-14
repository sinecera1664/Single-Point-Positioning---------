#include "single_point.h"
// ��������
// =========================================
/**
* ������Ч�Լ��
* tSatSurface��toe
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
* ����λ�á��ٶȡ��Ӳ����
* @param tSatSurface ���Ǳ���ʱ��
* @param eph ��������
* @param res  ���Ƕ�λ���������λ�á��ٶȡ��Ӳ����
* @return bool ��λ�Ƿ�ɹ�
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

	//�������
	double A = eph.A;

	//ƽ���˶����ٶ�
	double n0 = sqrt(GM / (A * A * A));

	//����������ο���Ԫ��ʱ��
	double tk = (double(tSatSurface.week) - eph.weekReference) * 604800.0 + tSatSurface.sow - eph.toe;
	double tk2 = (double(tSatSurface.week) - eph.weekReference) * 604800.0 + tSatSurface.sow - eph.toc;
	res.satClkD = eph.af0 + eph.af1 * tk2 + eph.af2 * tk2 * tk2;
	tk = tk - res.satClkD;

	//ƽ���˶����ٶȸ���
	double n = n0 + eph.deltaN;

	//ƽ���ǵ�
	double Mk = eph.M0 + n * tk;

	//����ƫ���ǵ�
	double e = eph.ecc;
	double Ek = Mk;
	double Etmp;
	double diff = 1;
	//������30�ε���
	for (int i = 0; i < 30; i++)
	{
		
		if (diff > 1e-10)//������λ�õ�Ӱ�� 3e7 * 1e-10
		{
			Etmp = Ek;
			Ek = Mk + e * sin(Etmp);
		}

		else
		{
			break;
		}
		//�жϵ����Ƿ����
		if (i == 30 && diff > 1e-10)
		{
			printf("Iteration incomplete!\n");
			exit(EXIT_FAILURE);
		}
		diff = fabs(Ek - Etmp);
	}
	//������
	double Vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
	//�����Ǿ�
	double faiK = Vk + eph.w;
	//���׵��͸���
	double deltaU = eph.cus * sin(2 * faiK) + eph.cuc * cos(2 * faiK);
	double deltaR = eph.crs * sin(2 * faiK) + eph.crc * cos(2 * faiK);
	double deltaI = eph.cis * sin(2 * faiK) + eph.cic * cos(2 * faiK);
	//�����������Ǿࡢ�򾶡�������
	double Uk = faiK + deltaU;
	double Rk = A * (1 - e * cos(Ek)) + deltaR;
	double Ik = eph.I0 + eph.Idot * tk + deltaI;
	//���ƽ���ϵ�λ��
	double xk = Rk * cos(Uk);
	double yk = Rk * sin(Uk);
	double Wk; //OmegaK

	bool isGeo = ((res.sys == BDS) && (eph.PRN < 6 || 58 < eph.PRN));

	// ��Geo�������⴦��
	Wk = eph.W0 + (eph.Wdot - (1.0 - isGeo) * WE) * tk - WE * (res.sys == GPS ? eph.toe : eph.toe - 14.0);

	//�ع�����ϵ�µ�λ��
	res.satP[0] = xk * cos(Wk) - yk * cos(Ik) * sin(Wk);
	res.satP[1] = xk * sin(Wk) + yk * cos(Ik) * cos(Wk);
	res.satP[2] = yk * sin(Ik);

	//����۸���
	res.satClkD += F * e * sqrt(A) * sin(Ek);

	res.tgd = eph.tgd[0];

	double Ekdot = n / (1 - e * cos(Ek));
	double Faikdot = sqrt((1 + e) / (1 - e)) * pow(cos(Vk / 2) / cos(Ek / 2), 2) * Ekdot;
	double Ukdot = 2 * (eph.cus * cos(2 * faiK) - eph.cuc * sin(2 * faiK)) * Faikdot + Faikdot;
	double Rkdot = A * e * sin(Ek) * Ekdot + 2 * (eph.crs * cos(2 * faiK) - eph.crc * sin(2 * faiK)) * Faikdot;
	double Ikdot = eph.Idot + 2 * (eph.cis * cos(2 * faiK) - eph.cic * sin(2 * faiK)) * Faikdot;

	//Geo�������⴦��
	double Wkdot = eph.Wdot - (1.0 - isGeo) * WE;

	Matrixd R({ {cos(Wk), (-1) * sin(Wk) * cos(Ik), (-1) * (xk * sin(Wk) + yk * cos(Wk) * cos(Ik)), yk * sin(Wk) * sin(Ik)},
								{sin(Wk), cos(Wk) * cos(Ik), (xk * cos(Wk) - yk * sin(Wk) * cos(Ik)), yk * cos(Wk) * sin(Ik)},
								{ 0, sin(Ik), 0, yk * cos(Ik)} });
	double xkdot = Rkdot * cos(Uk) - Rk * Ukdot * sin(Uk);
	double ykdot = Rkdot * sin(Uk) + Rk * Ukdot * cos(Uk);


	Matrixd R_v({ {xkdot}, {ykdot}, {Wkdot}, {Ikdot} });
	Matrixd satV = R * R_v;
	satV.OutToArray(3, res.satV);

	//Geo�������⴦��
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

	//����
	res.satClkV = eph.af1 + 2 * eph.af2 * tk2 + F * e * sqrt(A) * cos(Ek) * Ekdot;

	res.isValid = true;
	return true;
}

/**
* ���������
* ����40km��������Ч���룬��ȥ������0
* @param H ��վ�߶�
* @param E ���Ǹ߶Ƚ�
* @return double �����������
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
* ������ת,������ת
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
* IF���α��
* @param psr α��
* @param sys ����ϵͳ
* @return double IF���α��
*
*/
double IF(double* psr, NavSys& sys)
{
    double freq1 = sys == GPS ? L1 : B1;
    double freq2 = sys == GPS ? L2 : B3;
    return (freq1 * freq1 * psr[0] - freq2 * freq2 * psr[1]) / (freq1 * freq1 - freq2 * freq2);
}
// �����������
// =========================================



/**
* ����SatPos������IF��ϸ����������Ƕ�λ���������SPP������
* @param obs ����۲����ݣ��������������
* @param gps GPS����
* @param bds BDS����
*/
void InitRes(Obs& obs, Ephem* gps, Ephem* bds)
{
	double psrIF;
	GPSTime tSatSurface;
	obs.validNum = 0;

	for (int i = 0; i < obs.satObserved; i++)
	{
		//��Ƶ���������
		if (!obs.sat[i].isDualBand)
		{
			obs.validSat[i].isValid = false;   
			continue;                 //������һ��ѭ��
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
* SPP��λ
* @param  obs �۲�����
* @param  resSPP ���㶨λ���
* @return bool ��λ�Ƿ�ɹ�
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

	//�ж�������Ŀ��Ч��
	bool isSatNumValid = true;  	
	bool isSingleSys;
	Matrixd matQ, matW, matB, matX;

	//������10�ε���
	for (size_t j = 0; j < 10; j++)
	{
		while (diff > 1e-8 )
		{			                            
			int m = 0;//��i��ͬ��
			int gpsN = 0, bdsN = 0;

			for (int i = 0; i < obs.satObserved; i++)
			{
				if (!obs.validSat[i].isValid)
				{
					continue;
				}

				sys = obs.validSat[i].sys;
				//�źŴ���ʱ��
				deltaT = (obs.validSat[i].pseu / LIGHT_SPEED) + obs.validSat[i].satClkD - (sys == GPS ? resSPP.pos[3] / LIGHT_SPEED : resSPP.pos[4] / LIGHT_SPEED);
				//����ת��
				EarthRotate(deltaT, obs.validSat[i].satP, sys, satPosAfterErotate);
				orientationBLH = XYZToBLH(resSPP.pos, sys);

				//����elev���߶Ƚ�
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

				//��Ч���ж����߶Ƚ�С��20����ȥ
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
				//�������ӳ�
				trop = Hopefield(orientationBLH.H, elev);

				//�����û���ֵ�õ��ļ��ξ���
				rho = Distance(resSPP.pos, satPosAfterErotate);

				//�������� W ����            
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

			
			//��ʼ��С������⣬����������nValidNum����
			matW = Matrixd(obs.validNum, 1);			
			for (size_t i = 0; i < obs.validNum; i++)
			{
				matW[i][0] = W[i];
			}

			//��ϵͳ�ж�,B����ϵͳ����Ϊ4
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
			
			//λ�á��Ӳ�����
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

			//����һ��վ��ľ���
			diff = sqrt(matX[0][0] * matX[0][0] + matX[1][0] * matX[1][0] + matX[2][0] * matX[2][0]);

			//����δ���
			if (j == 10 && diff > 1e-8)
			{
				return false;
			}			
		}
	}
	//��������
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
* �������
* @param  obs �۲�����
* @param  resSPV ������ٽ��
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
	
	int m = 0;         //���ڼ���
	for (int i = 0; i < obs.satObserved; i++)
	{
		if (!obs.validSat[i].isValid)
		{
			continue; //��һ������
		}

		sys = obs.validSat[i].sys;

		deltaT = (obs.validSat[i].pseu / LIGHT_SPEED) + obs.validSat[i].satClkD -
			     (sys == GPS ? resSPV.pos[3] / LIGHT_SPEED : resSPV.pos[4] / LIGHT_SPEED);
		
		//������ת
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
		//dopp��λת��Ϊm/s
		W[m] = -1 * obs.sat[i].dopp[0] * lambda - (rhoDerivative - LIGHT_SPEED * obs.validSat[i].satClkV);  

		for (int j = 0; j < 3; j++)
		{
			B[m * 4 + j] = (resSPV.pos[j] - satPosAfterErotate[j]) / rho;
		}
		B[m * 4 + 3] = 1;
		m++;
	}

	//��С�������
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