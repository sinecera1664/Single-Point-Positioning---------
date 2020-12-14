#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>


namespace
{
    //define NDEBUG
    //#define _IS_SQUARE_
    //#define _NOT_EMPTY_
    bool PRINT_WITH_MORE_BLANK = true;//等待改进，改变输出的换行数   

    // 随机数生成器，使用 SetRand
    static std::mt19937* _genPtr;
    static std::normal_distribution<double> normDis(0, 1);
    static std::uniform_real_distribution<double> unifDoubleDis(0, 1);
    static std::uniform_real_distribution<float> unifFloatDis(0, 1);
    static std::uniform_int_distribution<int> unifIntDis(0, 65535);

    // 辅助代码
    // =========================================
    // 用来判断类别是否相同
    // =========================================
    template <typename T, typename U>
    struct SameType
    {
        static const bool isSame = false;
    };

    template <typename T>
    struct SameType<T, T>
    {
        static const bool isSame = true;
    };

    class Exception
    {
    public:
        explicit Exception(const std::string& _m) : message(_m)
        {
        }

        void printMessage() const
        {
            std::cout << message << std::endl;
        }

    private:
        std::string message;
    };

    // 辅助代码结束
    // =========================================

}



namespace MyMatrix
{
    template <typename T>
    class Matrix;

    template<typename S>
    std::ostream& operator << (std::ostream& out, const Matrix<S>& mat);


    template <typename T>
    class Matrix
    {
    public:
        //用于构造
        Matrix() = default;
        Matrix(Matrix<T>&& other);
        Matrix(const Matrix<T>& other);
        Matrix(size_t _x);
        Matrix(size_t _x, size_t _y);
        Matrix(std::vector<std::vector<T> > dvec);

        //工厂函数
        static Matrix<T> Eye(size_t _x, size_t _y);

        //生成随机矩阵
        static void inline SetRand();
        static void inline SetRand(std::mt19937* _p);

        static T inline Randn();
        static Matrix<T> inline Randn(size_t n);
        static Matrix<T> inline Randn(size_t r, size_t c);

        static T inline Rand();
        static Matrix<T> inline Rand(size_t n);
        static Matrix<T> inline Rand(size_t r, size_t c);

        void inline Clear();
        void inline Swap(Matrix<T>& rhs);
        void inline SetZeros();

        //用于获取数据
        std::vector<std::vector<T> >& GetData();
        const std::vector<std::vector<T> >& GetData() const;
        size_t inline Col() const;
        size_t inline Row() const;
        Matrix<double> ToDouble() const;
        T inline GetNum(size_t n) const;
        void OutToArray(size_t size, T* array)const;

        //用于判断
        bool inline IsEmpty() const;
        bool inline IsSquare() const;

        //向量运算
        static T inline Avg(const std::vector<T>& vec);
        static T VecDotProduct(const std::vector<T> lhs, std::vector<T> rhs);

        //矩阵运算
        T inline Max() const;
        T inline Min() const;
        T inline Avg() const;

        Matrix<T> Transposition() const;
        Matrix<double> Inversion() const;
        void Mul(Matrix<T>& ret, const Matrix<T>& other) const;
        void SMul(Matrix<T>& ret, const Matrix<T>& other) const;
        void StrassenMul(size_t rs, size_t re, size_t cs, size_t ce, const Matrix<T>& other, Matrix<T>& ret) const;

        //运算符重载
        std::vector<T>& operator [] (size_t index);
        const std::vector<T>& operator [] (size_t index) const;
        Matrix<T> operator = (const Matrix<T>& other);
        Matrix<T> operator = (Matrix<T>&& other);
        Matrix<T> inline operator + (const Matrix<T>& other) const;
        Matrix<T> inline operator - (const Matrix<T>& other) const;
        Matrix<T> inline operator += (const Matrix<T>& other);
        Matrix<T> inline operator -= (const Matrix<T>& other);
        Matrix<T> inline operator * (const Matrix<T>& other) const;
        Matrix<T> inline operator *= (const Matrix<T>& other);
        /**
        * 友元模板实现
        */
        template <typename S>
        friend std::ostream& operator << (std::ostream& out, const Matrix<S>& mat);

        //输出矩阵

        //void PrintMatrix() const;
        //static void PrintMatrix(const Matrix<T>& mat);




    private:
        std::vector<std::vector<T> > data;
    };

    /**
     * 构造函数：移动构造
     */
    template <typename T>
    Matrix<T>::Matrix(Matrix<T>&& other)
    {
        data.swap(other.data);
    }

    /**
     * 构造函数：拷贝构造
     */
    template <typename T>
    Matrix<T>::Matrix(const Matrix<T>& other)
    {
        data = other.GetData();
    }

    /**
     * 构造函数：建立一个n行的空矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(size_t _x)
    {
        std::vector<std::vector<T> > temp(_x);
        data = temp;
    }

    /**
     * 构造函数：根据行列数构造矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(size_t _x, size_t _y)
    {
        std::vector<std::vector<T> > temp(_x, std::vector<T>(_y));
        data = temp;
    }

    /**
     * 构造函数：调用vector拷贝方法，深拷贝
     */
    template <typename T>
    Matrix<T>::Matrix(std::vector<std::vector<T> > dvec)
    {
        data = dvec;
    }

    /**
     * 工厂函数：生成一个指定大小的单位阵
     * @param  _x 行数
     * @param  _y  列数
     * @return  一个单位矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Eye(size_t _x, size_t _y)
    {
        Matrix<T> mat(_x, _y);

        for (size_t i = 0; i < _x; i++)
        {
            for (size_t j = 0; j < _y; j++)
            {
                if (i == j)
                {
                    mat.data[i][j] = 1;
                }
            }
        }

        return mat;
    }

    /**
     * 初始化单例随机数生成器
     */
    template <typename T>
    void Matrix<T>::SetRand()
    {
        if (!_genPtr)
        {
            _genPtr = new std::mt19937(std::time(NULL));
        }
    }

    /**
     * 传入初始化单例随机数生成器
     */
    template <typename T>
    void Matrix<T>::SetRand(std::mt19937* _p)
    {
        if (_genPtr)
        {
            delete _genPtr;
        }

        _genPtr = _p;
    }

    /**
     * 如果矩阵是浮点数，生成 0 - 1 之间的一个均匀分布随机数
     * 如果矩阵是整形，生成 0 - 65525 之间的一个均匀分布随机数
     */
    template <typename T>
    T Matrix<T>::Rand()
    {
        SetRand();

        if (SameType<T, double>::isSame)
        {
            return unifDoubleDis(*_genPtr);
        }
        else if (SameType<T, float>::isSame)
        {
            return unifFloatDis(*_genPtr);
        }
        else
        {
            return unifIntDis(*_genPtr);
        }
    }

    /**
     * 生成一个 n * n 均匀分布随机矩阵
     * 输入参数为矩阵的边长
     */
    template <typename T>
    Matrix<T> Matrix<T>::Rand(size_t n)
    {
        SetRand();
        Matrix<T> mat(n, n);

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                mat[i][j] = Rand();
            }
        }

        return mat;
    }

    /**
     * 生成一个均匀分布随机矩阵
     * 输入参数为矩阵的行数和列数
     */
    template <typename T>
    Matrix<T> Matrix<T>::Rand(size_t r, size_t c)
    {
        SetRand();
        Matrix<T> mat(r, c);

        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                mat[i][j] = Rand();
            }
        }

        return mat;
    }

    /**
     * 生成 0 - 1 之间的一个正态随机数
     */
    template <typename T>
    T Matrix<T>::Randn()
    {
        SetRand();
        return normDis(*_genPtr);
    }

    /**
     * 生成一个 n * n 正态随机矩阵
     * 输入参数为矩阵的边长
     */
    template <typename T>
    Matrix<T> Matrix<T>::Randn(size_t n)
    {
        SetRand();
        Matrix<T> mat(n, n);

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                mat[i][j] = normDis(*_genPtr);
            }
        }

        return mat;
    }

    /**
     * 生成一个正态随机矩阵
     * 输入参数为矩阵的行数和列数
     */
    template <typename T>
    Matrix<T> Matrix<T>::Randn(size_t r, size_t c)
    {
        SetRand();
        Matrix<T> mat(r, c);
        auto fun = std::bind(normDis, *_genPtr);

        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                mat[i][j] = normDis(*_genPtr);
            }
        }

        return mat;
    }

    /**
     * 矩阵清空
     * 只删除数据 不清空内存
     */
    template <typename T>
    void Matrix<T>::Clear()
    {
        data.clear();
    }

    /**
     * 调用vector的swap方法，和右端矩阵元素整体调换
     */
    template <typename T>
    void Matrix<T>::Swap(Matrix<T>& rhs)
    {
        data.swap(rhs.GetData());
    }

    /**
     * 矩阵清零
     */
    template <typename T>
    void Matrix<T>::SetZeros()
    {
        size_t n = Row();
        size_t m = Col();

        for (size_t i = 0; i < n; i++)
        {
            std::memset(&data[i][0], 0, sizeof(T) * m);
        }
    }

    /**
     * 获取 data 成员，可用于整块更新
     */
    template <typename T>
    std::vector<std::vector<T> >& Matrix<T>::GetData()
    {
        return data;
    }

    /**
     * 以 const 方式获取成员，可用于安全读
     */
    template <typename T>
    const std::vector<std::vector<T> >& Matrix<T>::GetData() const
    {
        return data;
    }

    /**
     * 矩阵的行数
     */
    template <typename T>
    size_t Matrix<T>::Row() const
    {
        return data.size();
    }

    /**
     * 矩阵的列数
     */
    template <typename T>
    size_t Matrix<T>::Col() const
    {
        if (data.size())
        {
            return data[0].size();
        }
        else
        {
            return 0;
        }
    }

    /**
     * 转换为 double 矩阵
     */
    template <typename T>
    Matrix<double> Matrix<T>::ToDouble() const
    {
        Matrix<double> mat(Row(), Col());

        // 未实现：如果矩阵是复数，则只取实数部分，忽略虚数部分
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                mat[i][j] = static_cast<double>(data[i][j]);
            }
        }

        return mat;
    }

    /**
    * 返回第n个数
    */
    template <typename T>
    T Matrix<T>::GetNum(size_t n) const
    {
#ifndef NUMCHECK
#define NUMCHECK
        assert(data.size());
        assert((Col() * Row()) >= n);
#endif
        if (n % Col() == 0)
        {
            return data[n / Col() - 1][Col() - 1];

        }
        return data[n / Col()][n % Col() - 1];
    }

    /**
    * 提取矩阵前n个元素到数组
    */
    template<typename T>
    void Matrix<T>::OutToArray(size_t size, T* array) const
    {
#ifndef SIZECHECK
#define SIZECHECK
        assert(data.size());
        assert((Col() * Row()) >= size);
#endif
        size_t k = 0;

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                array[k] = data[i][j];
                k++;
                if (k == size)
                {
                    break;
                }
            }
            if (k == size)
            {
                break;
            }
        }
    }

    /**
     * 判断是否是空矩阵
     * @return 1 : 0 -> 空矩阵 : 不是空矩阵
     */
    template <typename T>
    bool Matrix<T>::IsEmpty() const
    {
        return !data.size();
    }

    /**
     * 判断矩阵是否是方阵
     * @return 1 : 0 -> 方阵 : 不是方阵
     */
    template <typename T>
    bool Matrix<T>::IsSquare() const
    {
        return IsEmpty() ? 0 : data.size() == data[0].size();
    }

    /**
     * 一个向量的均值
     */
    template <typename T>
    T Matrix<T>::Avg(const std::vector<T>& vec)
    {
        T sum = 0;

        for (T var : vec)
        {
            sum += var;
        }

        return sum / vec.size();
    }

    /**
     * 静态函数：获取两个向量的点乘结果
     * @param  lhs 向量1
     * @param  rhs 向量2
     * @return     double:点乘的结果
     */
    template <typename T>
    T Matrix<T>::VecDotProduct(const std::vector<T> lhs, std::vector<T> rhs)
    {
        T ans = 0;

        for (decltype(lhs.size()) i = 0; i != lhs.size(); i++)
        {
            ans += lhs[i] * rhs[i];
        }

        return ans;
    }

    /**
     * 返回矩阵中最大的元素
     */
    template <typename T>
    T Matrix<T>::Max() const
    {
        if (!data.size())
        {
            return static_cast<T>(0);
        }

        T maxv = data[0][0];

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                maxv = data[i][j] < maxv ? maxv : data[i][j];
            }
        }

        return maxv;
    }

    /**
     * 返回矩阵中最小的元素
     */
    template <typename T>
    T Matrix<T>::Min() const
    {
        if (!data.size())
        {
            return static_cast<T>(0);
        }

        T minv = data[0][0];

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                minv = data[i][j] < minv ? data[i][j] : minv;
            }
        }

        return minv;
    }

    /**
     * 矩阵的均值
     */
    template <typename T>
    T Matrix<T>::Avg() const
    {
        if (IsEmpty())
        {
            return static_cast<T>(0);
        }

        T sum = 0;

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                sum += data[i][j];
            }
        }

        return sum / (Row() * Col());
    }

    /**
     * 获得矩阵的转置
     * @return 新的矩阵，内容为原矩阵的转置
     */
    template <typename T>
    Matrix<T> Matrix<T>::Transposition() const
    {
        decltype(data.size()) sizeRow = data.size();

        if (sizeRow == 0)
        {
            std::cerr << "error** Matrix<T>::Transposition -> Empty Matrix!" << std::endl;
        }

        using size_t = decltype(data.size());
        size_t sizeCol = data[0].size();

        Matrix tran(sizeCol, sizeRow);

        for (size_t i = 0; i < sizeRow; i++)
        {
            for (size_t j = 0; j < sizeCol; j++)
            {
                tran.data[j][i] = data[i][j];
            }
        }

        return tran;
    }

    /**
     * 高斯约当法求逆
     */
    template <typename T>
    Matrix<double> Matrix<T>::Inversion() const
    {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!IsEmpty());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(IsSquare());
#endif
        size_t i, j, k, len = Row();
        double maxVal, temp;
        //将A矩阵存放在临时矩阵中
        Matrix<double> TMat;

        if (SameType<T, double>::isSame)
        {
            TMat = *this;
        }
        else
        {
            TMat = this->ToDouble();
        }

        //初始化ans矩阵为单位阵
        Matrix<double> ans = Matrix<double>::Eye(Row(), Col());

        for (i = 0; i < len; i++)
        {
            //寻找主元
            maxVal = TMat[i][i];
            k = i;

            for (j = i + 1; j < len; j++)
            {
                if (std::abs(TMat[j][i]) > std::abs(maxVal))
                {
                    maxVal = TMat[j][i];
                    k = j;
                }
            }

            //如果主元所在行不是第i行，进行行交换
            if (k != i)
            {
                TMat[i].swap(TMat[k]);
                ans[i].swap(ans[k]);
            }

            //assert(cond2().real() < 1e11);
            //判断主元是否为0, 若是, 则矩阵A不是满秩矩阵,不存在逆矩阵
            //未实现
            // if (cond2().real() > 1e10)
            //{
            //     throw (Exception("ERROR** Matrix::Inversion -> there is no inverse matrix!"));
            // }

            //消去A的第i列除去i行以外的各行元素
            temp = TMat[i][i];

            for (j = 0; j < len; j++)
            {
                TMat[i][j] = TMat[i][j] / temp; //主对角线上的元素变为1
                ans[i][j] = ans[i][j] / temp;   //伴随计算
            }

            // 遍历行
            for (j = 0; j < len; j++)
            {
                // 不是第i行
                if (j != i)
                {
                    temp = TMat[j][i];

                    // 第j行元素 - i行元素 * j列i行元素
                    for (k = 0; k < len; k++)
                    {
                        TMat[j][k] -= TMat[i][k] * temp;
                        ans[j][k] -= ans[i][k] * temp;
                    }
                }
            }
        }
        return ans;
    }

    /**
     * 矩阵相乘
     * 普通算法
     */
    template <typename T>
    void Matrix<T>::Mul(Matrix<T>& ret, const Matrix<T>& other) const
    {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(Col() == other.Row());
#endif

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t k = 0; k < Col(); k++)
            {
                for (size_t j = 0; j < other.Col(); j++)
                {
                    ret[i][j] += (data[i][k] * other[k][j]);
                }
            }
        }

        return;
    }

    /**
     * 斯特拉森乘法主函数，两个 n * n 矩阵
     */
    template <typename T>
    void Matrix<T>::SMul(Matrix<T>& ret, const Matrix<T>& other) const
    {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(Col() == other.Row());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(IsSquare());
#endif
#ifndef OTHERISSQUARE
#define OTHERISSQUARE
        assert(other.IsSquare());
#endif
        size_t n = Row();
        StrassenMul(0, n, 0, n, other, ret);
    }

    /**
     * 斯特拉森乘法
     * 时间复杂度 n ^ 2.80，性能更好，但参数复杂
     */
    template <typename T>
    void Matrix<T>::StrassenMul(size_t rs, size_t re, size_t cs, size_t ce, const Matrix<T>& other, Matrix<T>& ret) const
    {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(Col() == other.Row());
#endif

        if (re - rs == 2 && ce - cs == 2)
        {
            size_t rs1 = rs + 1;
            size_t cs1 = cs + 1;
            T P1 = data[rs][cs] * (other[rs][cs1] - other[rs1][cs1]);
            T P2 = (data[rs][cs] + data[rs][cs1]) * other[rs1][cs1];
            T P3 = (data[rs1][cs] + data[rs1][cs1]) * other[rs][cs];
            T P4 = data[rs1][cs1] * (other[rs1][cs] - other[rs][cs]);
            T P5 = (data[rs][cs] + data[rs1][cs1]) * (other[rs][cs] + other[rs1][cs1]);
            T P6 = (data[rs][cs1] - data[rs1][cs1]) * (other[rs1][cs] + other[rs1][cs1]);
            T P7 = (data[rs][cs] - data[rs1][cs]) * (other[rs][cs] + other[rs][cs1]);
            ret[rs][cs] = P5 + P4 - P2 + P6;
            ret[rs][cs1] = P1 + P2;
            ret[rs1][cs] = P3 + P4;
            ret[rs1][cs1] = P1 + P5 - P3 - P7;
        }
        else if (re - rs < 2 || rs - rs < 2)
        {
            for (size_t i = rs; i < re; i++)
            {
                for (size_t k = cs; k < ce; k++)
                {
                    for (size_t j = cs; j < ce; j++)
                    {
                        ret[i][j] += data[i][k] * other[k][j];
                    }
                }
            }
        }
        else
        {
            size_t rm = rs + ((re - rs) / 2);
            size_t cm = cs + ((ce - cs) / 2);
            StrassenMul(rs, rm, cs, cm, other, ret);
            StrassenMul(rm, re, cs, cm, other, ret);
            StrassenMul(rs, rm, cm, ce, other, ret);
            StrassenMul(rm, re, cm, ce, other, ret);
        }
    }

    /**
     * 取矩阵中的某一行
     */
    template <typename T>
    std::vector<T>& Matrix<T>::operator [] (size_t index)
    {
        assert(index >= 0 && index < data.size());
        return data[index];
    }

    /**
     * 常对象取矩阵中的某一行
     */
    template <typename T>
    const std::vector<T>& Matrix<T>::operator [] (size_t index) const
    {
        assert(index >= 0 && index < data.size());
        return data[index];
    }

    /**
     * 拷贝矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator = (const Matrix<T>& other)
    {
        data = other.GetData();
        return *this;
    }

    /**
     * 移动矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator = (Matrix<T>&& other)
    {
        data.swap(other.GetData());
        return *this;
    }

    /**
     * 矩阵相加
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator + (const Matrix<T>& other) const
    {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(Row() == other.Row() && Col() == other.Col());
#endif

        Matrix<T> ret(Row(), Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                ret[i][j] = data[i][j] + other[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵相加赋值
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator += (const Matrix<T>& other)
    {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(Row() == other.Row() && Col() == other.Col());
#endif

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] += other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵相减
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator - (const Matrix<T>& other) const
    {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(Row() == other.Row() && Col() == other.Col());
#endif

        Matrix<T> ret(Row(), Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                ret[i][j] = data[i][j] - other[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵相减赋值
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator -= (const Matrix<T>& other)
    {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(Row() == other.Row() && Col() == other.Col());
#endif

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] -= other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵相乘
     * 普通算法
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator * (const Matrix<T>& other) const
    {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(Col() == other.Row());
#endif
        Matrix<T> ret(Row(), other.Col());

        // 对称矩阵使用 斯特拉森乘法
        if ((Row() == Col()) && (other.Row() == other.Col()))
        {
            SMul(ret, other);
        }
        else
        {
            // 普通乘法
            Mul(ret, other);
        }

        return ret;
    }

    /**
     * 矩阵 *=
     * 普通算法
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator *= (const Matrix<T>& other)
    {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(Col() == other.Row());
#endif
        Matrix<T> ret(Row(), other.Col());

        // 大的对称矩阵使用 斯特拉森乘法
        if (Row() > 50 && Row() == Col())
        {
            SMul(ret, other);
        }
        else
        {
            // 普通乘法
            Mul(ret, other);
        }

        return ret;
    }

    /**
     * 友元模板实现输出
     */
    template <typename S>
    std::ostream& operator << (std::ostream& out, const Matrix<S>& mat)
    {
        //判断要输出的数是否是 double 类型
        bool dFlag = SameType<S, double>::isSame;


        if (dFlag)
        {
            using std::ios;
            using std::setprecision;
            out << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(4);
        }

        for (size_t i = 0; i != mat.data.size(); i++)
        {
            for (size_t j = 0; j != mat.data[i].size(); j++)
            {
                if (dFlag)
                {
                    out << std::setw(12) << mat.data[i][j] << ' ';
                }
                else
                {
                    out << mat.data[i][j] << ' ';
                }
            }

            if (i < mat.data.size() - 1)
            {
                out << '\n';
            }
        }

        if (PRINT_WITH_MORE_BLANK)
        {
            out << std::endl;
        }

        out << std::endl;
        return out;
    }



    /*
    template<typename T>
    void Matrix<T>::PrintMatrix() const
    {
        //判断要输出的数是否是 double 类型
        bool dFlag = SameType<T, double>::isSame;

        if (dFlag)
        {
            using std::ios;
            using std::setprecision;
            std::cout << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(4);
        }

        for (size_t i = 0; i != data.size(); i++)
        {
            for (size_t j = 0; j != data[i].size(); j++)
            {
                if (dFlag)
                {
                    std::cout << std::setw(12) << data[i][j] << ' ';
                }
                else
                {
                    std::cout << data[i][j] << ' ';
                }
            }

            if (i < data.size() - 1)
            {
                std::cout << '\n';
            }
        }

        std::cout << std::endl<<std::endl;
    }


    template<typename T>
    void Matrix<T>::PrintMatrix(const Matrix<T> & mat)
    {
        //判断要输出的数是否是 double 类型
        bool dFlag = SameType<T, double>::isSame;

        if (dFlag)
        {
            using std::ios;
            using std::setprecision;
            std::cout << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(4);
        }

        for (size_t i = 0; i != mat.data.size(); i++)
        {
            for (size_t j = 0; j != mat.data[i].size(); j++)
            {
                if (dFlag)
                {
                    std::cout << std::setw(12) << mat.data[i][j] << ' ';
                }
                else
                {
                    std::cout << mat.data[i][j] << ' ';
                }
            }

            if (i < mat.data.size() - 1)
            {
                std::cout << '\n';
            }
        }

        std::cout << std::endl<<std::endl;
    }
    */

    using Matrixd = Matrix<double>;
    using Matrixf = Matrix<float>;
    using Matrixi = Matrix<int>;

}

namespace mym = MyMatrix;

// ©2020 sinecera.All rights reserved