#pragma once
#include "basic_service.h"


void InitRes(Obs& obs, Ephem* gps, Ephem* bds);

bool GetSPP(Obs& obs, Orientation& resSPP);

void GetSPV(Obs& obs, Orientation& resSPV);

