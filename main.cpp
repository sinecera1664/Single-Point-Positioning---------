// Main.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。

#include "single_point.h"
#include "sockets.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iosfwd>
using namespace std;

int main()
{   
    SOCKET sc;
#ifndef SOCKETCHECK
#define SOCKETCHECK
    assert(OpenSocket(sc, "47.114.134.129", 6000));
#endif 

    ofstream output("output.txt");  //输出解算出的数据
    int lengthTail = 0;
    int lengthBufTCP = 0;
    unsigned char bufTCP[MAXRAWLEN_SOCKET], buf[2 * MAXRAWLEN_SOCKET];

    double oriXYZ[3] = { -2267800.0690, 5009341.6472, 3220989.9179 };
    double oriBLH[3] = { 30.528437963 * PI / 180, 114.356937696 * PI / 180, 38.075359 };
    double tarENU[3];
    
    //从0开始迭代,默认构造为0
    Obs obs;
    Ephem gps[GPSSATNUM], bds[BDSSATNUM];
    Orientation res;
    size_t count = 0;
    ios::sync_with_stdio(false);

    while (true)
    {
        Sleep(980);

        if ((lengthBufTCP = recv(sc, (char*)bufTCP, MAXRAWLEN_SOCKET, 0)) > 0)
        {
            //lengthTail=buf中指针位置
            if ((lengthTail + lengthBufTCP) > sizeof(buf))
            {
                lengthTail = 0;
            }

            memcpy(buf + lengthTail, bufTCP, lengthBufTCP);
            lengthTail += lengthBufTCP;
            
            if (DecodeBuffer(buf, lengthTail, obs, gps, bds) != true)
            {
                continue;
            }

            InitRes(obs, gps, bds);
            
           
            if (!GetSPP(obs, res))
            {
                continue;
            }

            //转换为enu坐标
            XYZToENU(oriXYZ, res.pos, oriBLH, tarENU);

            /*
            cout << setiosflags(ios_base::fixed);
            cout << int(res.tObserve.sow) << "   ";
            for (size_t i = 0; i < 3; i++)
            {
                cout << res.pos[i] << "  ";
            }      
            
            for (const auto &var : tarENU)
            {
                cout << var << "   ";
            }
          
            cout << endl;
            */

            cout << ++count << endl;
            output <<res.tObserve.sow << "   ";
            for (size_t i = 0; i < 3; i++)
            {
                output << res.pos[i] << "  ";
            }

            for (auto var : tarENU)
            {
                output << var << "   ";
            }

            output << endl;
            
        }
    }

    
    CloseSocket(sc);
    ios::sync_with_stdio(true);
    output.close();
    return 0;
}

