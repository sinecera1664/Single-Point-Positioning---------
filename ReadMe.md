

## 1、程序设计概述

​    本程序使用C++实现，并在MSVC上编译通过，可较长时间运行。实现了读取网口实时单点定位和测速的功能，并对中间过程进行封装，只对外暴露解码和单点定位测速方法，程序流程如下：

~~~flow
st1=>start: 开始
op1=>operation: 初始化socket对象
cond1=>condition: 连接成功(是或否?)
io1=>inputoutput: 读取buffer
op2=>operation: DecodeBuffer
op3=>operation: InitRes 
                生成卫星定位结果
cond2=>condition: SPP成功（是或否）
op4=>operation: GetSPV
io2=>inputoutput: 输出XYZ坐标和enu坐标
e=>end: 终止

st1->op1->cond1
cond1(yes)->io1->op2->op3(right)->cond2
cond2(no)->io1
cond2(yes)->op4(right)->io2(top)->io1
cond1(no)->e
~~~

​    

## 2.算法与程序设计

​     程序由5个部分组成：

### 2.1 main.cpp

​    程序开始执行和结束的地方，仅含有main函数，通过调用方法实现单点定位，具体内容如图1-1所述。

### 2.2 basic_service.h和basic_service.cpp

​    包含所有数据类型定义，宏定义，和基础的方法。

​     方法如下：分别用于求三维坐标的距离，坐标转换和读取buffer并解码

 ~~~c++
double Distance(double* left, double* right);

BLH XYZToBLH(double* var, NavSys sys);

void XYZToENU(double* oriXYZ, double* tarXYZ, double* oriBLH, double* tarENU);

bool DecodeBuffer(unsigned char* buf, int& lengthTail, Obs& obs, Ephem* gps, Ephem* bds);
 ~~~



### 2.3 matrix.h

​    仅含有用于矩阵计算的模板类（类名为Matrix）的定义,  实现了矩阵转置、求逆、运算符重载等功能，命名空间为MyMatrix，尽可能不污染变量域的命名，使用的命名仅有：

~~~c++
 template <typename T> class Matrix;
 using Matrixd = Matrix<double>;
 using Matrixf = Matrix<float>;
 using Matrixi = Matrix<int>;
~~~



### 2.4 socket.h

   用于初始化和关闭socket对象，实现网口的读取，含有以下两个方法：

~~~c++
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);

void CloseSocket(SOCKET& sock);
~~~



### 2.5 single_point.h和single_point.cpp

程序以及算法设计的主要部分， single_point.h仅对外暴露三个方法，分别用于卫星定位，SPP和SPV：

~~~c++
void InitRes(Obs& obs, Ephem* gps, Ephem* bds);   

bool GetSPP(Obs& obs, Orientation& resSPP);

void GetSPV(Obs& obs, Orientation& resSPV);
~~~

single_point.cpp含有的方法如下：

#### 2.5.1  星历有效性检测

~~~c++
/**
* 星历有效性检测
* 星历是否存在
* tSatSurface和toe差是否符合要求
* GPS限额7200，Bds限额3600
*/
bool IsSatVaild(const GPSTime& tSatSurface, const Ephem& ehm)
~~~

#### 2.5.2 Geo卫星定位

~~~c++
/**
* 获取Geo卫星位置、速度、钟差、钟速
* @param tSatSurface 卫星表面时间
* @param eph 卫星星历
* @param res  卫星定位结果，包含位置、速度、钟差、钟速
* @return bool 定位是否成功
*/
bool SatGeoPos(GPSTime& tSatSurface, Ephem& eph, SatMedian& res)
~~~

#### 2.5.3 卫星定位

~~~c++
/**
* 调用IsSatVaild和SatGeoPos方法
* 获取GPS和BDS卫星位置、速度、钟差、钟速
* @param tSatSurface 卫星表面时间
* @param eph 卫星星历
* @param res  卫星定位结果，包含位置、速度、钟差、钟速
* @return bool 定位是否成功
*/
bool SatPos(GPSTime& tSatSurface, Ephem& eph, SatMedian& res)
~~~

###### 1.位置计算：

计算轨道长半轴
$$
\quad A=(\sqrt{A})^{2}
$$
计算平均运动角速度
$$
\quad n_{0}=\sqrt{\frac{\mu}{A^{3}}}
$$
计算相对于星历参考历元的时间
$$
\quad t_{k}=t-t_{o e}
$$


对平均运动角速度进行改正
$$
n=n_{0}+\Delta n
$$
计算平近点角
$$
M_{k}=M_{0}+n t_{k}
$$
计算偏近点角（利用下面的开普勒方程，迭代 求解
$$
M_{k}=E_{k}-e \sin E_{k}
$$
计算真近点角
$$
v_{k}=\arctan \left(\frac{\sin v_{k}}{\cos v_{k}}\right)=\arctan \left(\frac{\sqrt{1-e^{2}} \sin E_{k} /\left(1-e \cos E_{k}\right)}{\left(\cos E_{k}-e\right) /\left(1-e \cos E_{k}\right)}\right)
$$
计算升交角距
$$
\Phi_{k}=v_{k}+\omega
$$
计算二阶调和改正数 

​      计算升交角距的改正数
$$
\delta u_{k}=C_{u s} \sin 2 \Phi_{k}+C_{u c} \cos 2 \Phi_{k}
$$
​      计算向径的改正数
$$
\delta r_{k}=C_{r s} \sin 2 \Phi_{k}+C_{r c} \cos 2 \Phi_{k}
$$
​      计算轨道倾角改正数
$$
\delta i_{k}=C_{i s} \sin 2 \Phi_{k}+C_{i c} \cos 2 \Phi_{k_{k}}
$$
计算经过改正的升交角距
$$
u_{k}=\Phi_{k}+\delta u_{k}
$$
计算经过改正的向径
$$
r_{k}=A\left(1-e \cdot \cos E_{k}\right)+\delta r_{k}
$$
计算经过改正的轨道倾角
$$
i_{k}=i_{0}+\delta i_{k}+(I D O T) \cdot t_{k}
$$
计算卫星在轨道平面上的位置
$$
\left\{\begin{array}{l}
x_{k}^{\prime}=r_{k} \cos u_{k} \\
y_{k}^{\prime}=r_{k} \sin u_{k}
\end{array}\right.
$$
计算改正后的升交点经度
$$
\Omega_{k}=\Omega_{0}+\left(\dot{\Omega}-\dot{\Omega}_{e}\right) \cdot t_{k}-\dot{\Omega}_{e} \cdot t_{o e}
$$
计算在地固坐标系下的位置
$$
\left\{\begin{array}{l}
x_{k}=x_{k}^{\prime} \cos \Omega_{k}-y_{k} ^{\prime} \cos i_k \sin \Omega_{k} \\
y_{k}=x_{k} ^{\prime}\sin \Omega_{k}+y_{k} ^{\prime} \cos i_k \cos \Omega_{k} \\
z_{k}=y_{k} ^{\prime} \sin i_k 
\end{array}\right.
$$
钟差计算
$$
\begin{array}{l}

t=t_{\mathrm{SV}}-\left(\Delta t_{\mathrm{SV}}\right)_{\mathrm{L} 1} \\
\left(\Delta t_{\mathrm{SV}}\right)_{\mathrm{L} 1}=C l k B i a s+C l k D r i f t \cdot\left(t-t_{o c}\right)+C l k D r i f t R a t e \cdot\left(t-t_{o c}\right)^{2}+\Delta t_{\mathrm{r}}-T G D \\
\Delta t_{\mathrm{r}}=F e \sqrt{A} \sin E_{k} \\
F=\frac{-2 \sqrt{\mu}}{c^{2}}=-4.442807633 \cdot 10^{-10}\left(\mathrm{sec} \cdot \mathrm{m}^{-1 / 2}\right)
\end{array}
$$
卫星速度计算
$$
\begin{array}{l}
\dot{\mathbf{R}}=\left[\begin{array}{cccc}
\cos \Omega_{k} & -\sin \Omega_{k} \cos i_{k} & -\left(x_{k}^{\prime} \sin \Omega_{k}+y_{k}^{\prime} \cos \Omega_{k} \cos i_{k}\right) & y_{k}^{\prime} \sin \Omega_{k} \sin i_{k} \\
\sin \Omega_{k} & \cos \Omega_{k} \cos i_{k} & \left(x_{k}^{\prime} \cos \Omega_{k}-y_{k}^{\prime} \sin \Omega_{k} \cos i_{k}\right) & y_{k}^{\prime} \cos \Omega_{k} \sin i_{k} \\
0 & \sin i_{k} & 0 & y_{k}^{\prime} \cos i_{k}
\end{array}\right]
\end{array}
$$

$$
\begin{array}{l}
x_{k}^{\prime}=r_{k} \cos u_{k} & \dot{x}_{k}^{\prime}=\dot{r}_{k} \cos u_{k}-r_{k} \dot{u}_{k} \sin u_{k} \\
y_{k}^{\prime}=r_{k} \sin u_{k} & \dot{y}_{k}^{\prime}=\dot{r}_{k} \sin u_{k}+r_{k} \dot{u}_{k} \cos u_{k}
\end{array}
$$

$$
\begin{array}{l}
{\left[\begin{array}{c}
\dot{x}_{k} \\
\dot{y}_{k} \\
\dot{z}_{k}
\end{array}\right]=\dot{\mathbf{R}}\left[\begin{array}{c}
\dot{x}_{k}^{\prime} \\
\dot{y}_{k}^{\prime} \\
\dot{\Omega}_{k} \\
\dot{I}_{k}
\end{array}\right]}
\end{array}
$$

###### 2. 钟速计算


$$
\begin{array}{l}
\left(\Delta t_{\mathrm{SV}}\right)_{\mathrm{Ll}}=C l k D r i f t+2 C l k D r i f t R a t e \cdot\left(t-t_{o c}\right)+\Delta t_{\mathrm{r}}^{\prime} \\
\Delta t_{\mathrm{r}}^{\prime}=F e \sqrt{A} \cos E_{k} E_{k}^{\prime} \\
F=\frac{-2 \sqrt{\mu}}{c^{2}}=-4.442807633 \cdot 10^{-10}\left(\mathrm{sec} \cdot \mathrm{m}^{-1 / 2}\right)
\end{array}
$$


#### 2.5.4 对流层改正

~~~c++
/**
* 对流层误差
* 超出40km，超出有效距离，舍去并返回0
* @param H 测站高度
* @param E 卫星高度角
* @return double 对流层改正数
*/
double Hopefield(double H, double E)
~~~

$$
\begin{array}{l}
\Delta_{\text {Trop}}=\Delta_{d}+\Delta_{w}=\frac{K_{d}}{\sin \sqrt{E^{2}+6.25}}+\frac{K_{w}}{\sin \sqrt{E^{2}+2.25}} 
\end{array}
$$

$$
\begin{array}{l}
K_{d}=155.2 \cdot 10^{-7} \cdot \frac{p}{T}\left(h_{d}-H\right) 
\end{array}
$$

$$
\begin{array}{l}
K_{w}=155.2 \cdot 10^{-7} \cdot \frac{4810}{T^{2}} e\left(h_{w}-H\right) 
\end{array}
$$

$$
\begin{array}{l}
h_{d}=40136+148.72\left(T_{0}-273.16\right) \quad[\mathrm{m}] 
\end{array}
$$

$$
\begin{array}{l}
h_{w}=11000 \mathrm{~m} 
\end{array}
$$

$$
\begin{array}{l}
e=R H \cdot \exp \left(-37.2465+0.213166 T-0.000256908 T^{2}\right) 
\end{array}
$$

$$
\begin{array}{l}
T=T_{0}-0.0065 \cdot\left(H-H_{0}\right) 
\end{array}
$$

$$
\begin{array}{l}
p=p_{0} \cdot\left(1-0.0000226 \cdot\left(H-H_{0}\right)\right)^{5.225} 
\end{array}
$$

$$
\begin{array}{l}
R H=R H_{0} \cdot \exp \left(-0.0006396 \cdot\left(H-H_{0}\right)\right)
\end{array}
$$

#### 2.5.5  IF组合伪距

~~~c++
/**
* IF组合伪距
* @param psr 伪距
* @param sys 卫星系统
* @return double IF组合伪距
*
*/
double IF(double* psr, NavSys& sys)
~~~


$$
\begin{array}{l}
P_{G L_{i}}^{j}=\rho_{G}^{j}+c \cdot \delta t_{G}-c \cdot \delta t_{\mathrm{G}}^{j}+I_{\mathrm{GL}_{i}}^{j}+T^{j}+\varepsilon_{G L_{i}}^{j} \\
P_{C B_{i}}^{k}=\rho_{C}^{k}+c \cdot \delta t_{C}-c \cdot \delta t_{\mathrm{C}}^{k}+I_{C B_{i}}^{k}+T^{k}+\varepsilon_{C B_{i}}^{k}
\end{array}
$$

$$
\begin{array}{l}
P_{G, I F}^{j}=\rho_{G}^{j}+c \cdot \delta t_{G}-c \cdot \delta t_{\mathrm{G}}^{j}+T^{j}+\varepsilon_{G, I F}^{j}
\end{array}
$$

$$
\begin{array}{l}
P_{C, I F}^{j}=\rho_{C}^{j}+c \cdot \delta t_{C}-c \cdot \delta t_{\mathrm{C}}^{j}+T^{j}+\varepsilon_{C, I F}^{j} 
\end{array}
$$

$$
\begin{array}{l}
P_{G, I F}^{j}=\frac{f_{L_{1}}^{2}}{f_{L_{1}}^{2}-f_{L_{2}}^{2}} P_{G, L_{1}}^{j}-\frac{f_{L_{2}}^{2}}{f_{L_{1}}^{2}-f_{L_{2}}^{2}} P_{G, L_{2}}^{j} 
\end{array}
$$

$$
\begin{array}{l}
P_{C, I F}^{j}=\frac{f_{B_{1}}^{2}}{f_{B_{1}}^{2}-f_{B_{3}}^{2}} P_{C, B_{1}}^{j}-\frac{f_{B_{3}}^{2}}{f_{B_{1}}^{2}-f_{B_{3}}^{2}} P_{C, B_{3}}^{j}
\end{array}
$$

#### 2.5.6 地球旋转改正

~~~c++
/**
* 地球自转,坐标旋转，大概0.1s
*/
void EarthRotate(double deltaT, double* posSatOri, NavSys& sys, double* posSatRes)
~~~

#### 2.5.6 卫星位置初始化

~~~c++
/**
* 调用SatPos，生成IF组合改正所有卫星定位结果，包含部分有效性检测，用于SPP的输入
* @param obs 输入观测数据，初步处理，并输出
* @param gps GPS星历
* @param bds BDS星历
*/
void InitRes(Obs& obs, Ephem* gps, Ephem* bds)
~~~

#### 2.5.7 单点定位

~~~c++
/**
* SPP定位
* @param  obs 观测数据
* @param  resSPP 单点定位结果
* @return bool 定位是否成功
*/
bool GetSPP(Obs& obs, Orientation& resSPP)
~~~

###### 1. 观测方程线性化

$$
\begin{array}{l}
P_{i, 1 F}^{j}(t)+v_{i, I F}^{j}(t)=\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)+\delta t_{i}^{0}(t)-c \cdot \delta t_{i}^{j}\left(t_{R, i T}-\Delta t\right)+T^{j}\left(t_{R}\right) \\
+\frac{X_{R}^{0}-X_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} x_{R}+\frac{Y_{R}^{0}-Y_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} y_{R}\\+\frac{Z_{R}^{0}-Z_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} z_{R}+\delta t_{G}\left(t_{R, G P S T}\right)+\delta t_{C}\left(t_{R, B D T}\right) 
\end{array}
$$


$$
\begin{array}{l}
v_{i, I F}^{j}(t)=\frac{X_{R}^{0}-X_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} x_{R}+\frac{Y_{R}^{0}-Y_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} y_{R}\\+\frac{Z_{R}^{0}-Z_{i}^{j}\left(t_{R, i T}-\Delta t\right)}{\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)} z_{R}+\delta t_{G}\left(t_{R, G P S T}\right)+\delta t_{C}\left(t_{R, B D T}\right) \\
-\left\{P_{i, I F}^{j}(t)-\left[\rho_{i}^{j 0}\left(t_{R}, t_{R, i T}-\Delta t\right)+\delta t_{i}^{0}(t)-c \cdot \delta t_{i}^{j}\left(t_{R, i T}-\Delta t\right)+T^{j}\left(t_{R}\right)\right]\right\}
\end{array}
$$




###### 2. 最小二乘求解


$$
\begin{array}{l}
{v}={B} {x}-{w} 
\end{array}
$$

$$
\begin{array}{l}
\mathbf{v}=\left[\begin{array}{lllll}
v_{G P S, I F}^{1} & \ldots & v_{G P S, I F}^{a} & v_{B D S, I F}^{1} & \ldots & v_{B D S, I F}^{b}
\end{array}\right]^{T}
\end{array}
$$

$$
\begin{array}{l}
\mathbf{B}=\left[\begin{array}{ccccc}
l_{G P S}^{1} & m_{G P S}^{1} & n_{G P S}^{1} & 1 & 0 \\
\ldots & \ldots & \ldots & 1 & 0 \\
l_{G P S}^{a} & m_{G P S}^{a} & n_{G P S}^{a} & 1 & 0 \\
l_{B D S}^{1} & m_{B D S}^{1} & n_{B D S}^{1} & 0 & 1 \\
\ldots & \ldots & \ldots & 0 & 1 \\
l_{B D S}^{a} & m_{B D S}^{a} & n_{B D S}^{a} & 0 & 1
\end{array}\right]
\end{array}
$$

$$
\begin{array}{l}
\overrightarrow{\mathbf{x}}_{R}=\left[\begin{array}{lllll}
x_{R} & y_{R} & z_{R} & \delta t_{G} & \delta t_{C}
\end{array}\right]^{T}
\end{array}
$$

$$
\begin{array}{l}
\quad \mathbf{w}=\left[\begin{array}{llllll}
w_{G}^{1} & \ldots & w_{G}^{a} & w_{C}^{1} & \ldots & w_{C}^{b}
\end{array}\right]^{T}
\end{array}
$$

$$
\begin{array}{l}
\hat{\mathbf{x}}_{R}=\left(\mathbf{B}^{T} \mathbf{P B}\right)^{-1} \mathbf{B} \mathbf{P} \mathbf{w} 
\end{array}
$$



###### 3. 精度评定


$$
\begin{array}{l}
\hat{\sigma}_{0}=\sqrt{\frac{\mathbf{v}^{T} \mathbf{P v}}{s-5}} 
\end{array}
$$

$$
\begin{array}{l}
\mathbf{D}_{\hat{x} \hat{x}}=\hat{\sigma}_{0} \cdot \mathbf{Q}_{\hat{x} \hat{x}}, \quad \mathbf{Q}_{\hat{x} \hat{x}}=\left(\mathbf{B}^{T} \mathbf{P} \mathbf{B}\right)^{-1}
\end{array}
$$

$$
P D O P=\sqrt{q_{\hat{x} \hat{x}}+q_{\hat{y} \hat{y}}+q_{\hat{z}}}
$$



####  2.5.8 单点测速

~~~c++
/**
* 单点测速
* @param  obs 观测数据
* @param  resSPV 单点测速结果
*/
void GetSPV(Obs& obs, Orientation& resSPV)
~~~



###### 1. 观测方程线性化


$$
\begin{array}{l}
D_{i}^{j}(t)+v_{i}^{j}(t)=\dot{\rho}_{i}^{j}(t, t-\Delta t)+\delta \dot{t}_{i}(t)-c \cdot \delta \dot{t}^{j}(t-\Delta t)+\dot{d}_{i r o p}^{j}(t)+\dot{d}_{i o n}^{j}(t) 
\end{array}
$$

$$
\begin{array}{l}
\dot{\rho}_{i}^{j}(t, t-\Delta t)=\frac{\left(X^{j}-X_{i}\right)\left(\dot{X}^{j}-\dot{X}_{i}\right)+\left(Y^{j}-Y_{i}\right)\left(\dot{Y}^{j}-\dot{Y}_{i}\right)+\left(Z^{j}-Z_{i}\right)\left(\dot{Z}^{j}-\dot{Z}_{i}\right)}{\rho_{i}^{j}(t, t-\Delta t)} 
\end{array}
$$

$$
\begin{array}{l}
\rho_{i}^{j}(t, t-\Delta t)=\sqrt{\left(X^{j}(t-\Delta t)-X_{i}\right)^{2}+\left(Y^{j}(t-\Delta t)-Y_{i}\right)^{2}+\left(Z^{j}(t-\Delta t)-Z_{i}\right)^{2}}
\end{array}
$$

$$
\begin{array}{l}
l_{i}^{j}(t)=\frac{X_{i}-X^{j}(t-\Delta t)}{\rho_{i}^{j}(t, t-\Delta t)}, \quad m_{i}^{j}(t)=\frac{Y-Y^{j}(t-\Delta t)_{i}}{\rho_{i}^{j}(t, t-\Delta t)}, \quad n_{i}^{j}(t)=\frac{Z_{i}-Z^{j}(t-\Delta t)}{\rho_{i}^{j}(t, t-\Delta t)} 
\end{array}
$$

$$
\begin{array}{l}
w_{i}^{j}(t)=D_{i}^{j}(t)-\left(\dot{\rho}_{i}^{j}(t, t-\Delta t)-c \cdot \delta \dot{t}^{j}(t-\Delta t)\right) 
\end{array}
$$

$$
\begin{array}{l}
\quad v_{i}^{j}(t)=l_{i}^{j}(t) x_{i}+m_{i}^{j}(t) y_{i}+n_{i}^{j}(t) z_{i}+\delta \dot{t}_{i}(t)-w_{i}^{j}(t)
\end{array}
$$



###### 2.最小二成求解和精度评定

$$
\begin{aligned}
&\mathbf{v}=\mathbf{B} \mathbf{x}-\mathbf{w}\\
\end{aligned}
$$

$$
\begin{aligned}
&\mathbf{v}=\left[\begin{array}{llll}
v_{i}^{1} & v_{i}^{2} & \ldots & v_{i}^{s}
\end{array}\right]^{T}, \quad \mathbf{B}=\left[\begin{array}{cccc}
l_{i}^{1} & m_{i}^{1} & n_{i}^{1} & 1 \\
l_{i}^{2} & m_{i}^{2} & n_{i}^{2} & 1 \\
. & . & . & . \\
l_{i}^{s} & m_{i}^{s} & n_{i}^{s} & 1
\end{array}\right]
\end{aligned}
$$

$$
\begin{aligned}
&\mathbf{x}=\left[\begin{array}{llll}
\dot{X}_{i} & \dot{Y}_{i} & \dot{Z}_{i} & \delta \dot{\rho}_{i}
\end{array}\right]^{T}
\end{aligned}
$$

$$
\begin{aligned}
\quad \mathbf{w}=\left[\begin{array}{llll}
w_{i}^{1} & w_{i}^{2} & \ldots & w_{i}^{s}
\end{array}\right]^{T}
\end{aligned}
$$

$$
\begin{aligned}
\delta \dot{\rho}_{i}=c \cdot \delta \dot{t}_{i}(t)\\
&\\
\end{aligned}
$$

$$
\begin{aligned}
&\hat{\mathbf{x}}=\left(\mathbf{B}^{T} \mathbf{P} \mathbf{B}\right)^{-1} \mathbf{B} \mathbf{P} \mathbf{w}\\
\end{aligned}
$$

$$
\begin{aligned}
&\hat{\sigma}_{0}=\sqrt{\frac{\mathbf{v}^{T} \mathbf{P v}}{s-4}}
\end{aligned}
$$

$$
\begin{aligned}
&\mathbf{D}_{\hat{x} \hat{x}}=\hat{\sigma}_{0} \cdot \mathbf{Q}_{\hat{x} \hat{x}}, \quad \mathbf{Q}_{\hat{x} \hat{x}}=\left(\mathbf{B}^{T} \mathbf{P} \mathbf{B}\right)^{-1}
\end{aligned}
$$



## 3.精度分析

### 3.1 绝对坐标XYZ误差



### 3.2 enu方向误差



### 3.3 误差分析

单点定位误差主要有三部分：

#### 3.3.1 观测误差

观察误差主要是接收机钟差及接收机天线相位误差。

接收机钟差是接收机的钟面时与G P S标准时间的差值。接收机钟差主要通过影响卫星位置和站星几何距离的计算来影响定位。由于接收机的钟差引起的卫星误差不同，只要估计好接收机的钟差就能消除卫星坐标计算的影响。
天线相位中心是指发射或者接收信号点，接收机天线相位误差是指天线相位中心与天线参考点之间的差值。天线相位中心的影响可以通过模型改正方法来消除。

地球自转改正。坐标参照系是随着地球自转而变化，如W G S-84属于地心地固坐标系，I T R F属于地固坐标系。卫星信号发射时刻和信号接收时刻所对应的地固系是不同的，所以在地固坐标系中计算卫星到接收机的几何距离时，就需要考虑地球自转的影响。

#### 3.3.2 卫星误差

卫星的误差主要包括卫星钟差和卫星轨道误差等。

卫星钟差是指卫星钟的频率漂移引起的卫星钟时间与标准G P S时间的差值。卫星的钟差会影响卫星坐标与站星几何距离的计算，一般要保证该值不大于一微秒。其改正方法一般是事先估计其大小，再用观测方程来消除其影响。卫星轨道误差是指卫星星历中给出的或者报据卫星星历计算出的卫星位罝与真实的卫星位置之间的差值。

卫星质量中心的坐标是精密星历给出的卫星坐标，但是卫星天线相位中心是指卫星发射信号的位置，这样就形成了误差，即卫星质心和卫星天线相位中心之间的偏差，这个变差就是卫星天线相位中心偏差。这个误差的改正方法类似接收机天线相位中心改正方法。

#### 3.3.3 传播误差

在传播过程中，对流层延迟和电离层延迟等会引起传播误差。

由于电磁波在电离层中传播的速度和路径会发生变化，因此利用信号传播时间和光速得到的距离观测值与信号源到接收机之间的真实几何距离就存在差异，这就引起了电离层延迟。目前用的比较的多的改正方法是采用国际电离层模型和K l o b u c h a r模型。本程序中使用了IF组合改正，故不需考虑。

此外，处于测站附近的反射物所反射的卫星信号进入接收机天线和直接来自卫星的信号产生干涉，从而引起误差，这个误差称为多路径效应误差。程序中未进行改正，但可通过选择合适的地址测量，让天线地点尽量远离反射体，另外也可以通过小波分析等方法来消除这个误差。



