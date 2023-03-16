#ifndef FRETCALCULATOR_H
#define FRETCALCULATOR_H

#include <math.h>
#include <random>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <numeric>
//#include <bits/stdc++.h>

#include <QDebug>
#include <QString>
#include <QDir>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define ld long double

//Switches of Debug
#define DEBUG_ALL           1
#define DEBUG_RUNDETAIL     0
#define DEBUG_RUNRESULT     1
#define DEBUG_PICKPOINT     0
#define DEBUG_LOADIMAGE     0
#define DEBUG_BATCH         1
#define DEBUG_VIEW          0
//Define : Max Gray of Different Channels
#define MAXGRAY_C16         65535
#define MAXGRAY_C8          255
//Define : Calculation Modes
#define CALCMODE_AB         0
#define CALCMODE_CD         1
#define CALCMODE_GK_AG      2   //反代法求GK, Anti-Generation
#define CALCMODE_GK_MF      3   //多质粒拟合求GK, Multiple-plasmid fitting
#define CALCMODE_EFRET      4
#define CALCMODE_3FRET      5
//Index of Array : Additional Functions
#define ADDIFUNC_NUMBER         2
#define ADDIFUNC_GENERATECSV    0
#define ADDIFUNC_SAVEMASK       1
//Index of Array : Pick Point Parameters
#define PICKPARAM_NUMBER            2
#define PICKPARAM_POINTNUMBER       0   //取点的数量，若该项取-1，则会统计像素点然后取点
#define PICKPARAM_PICKMODE          1   //取点模式，取值为如下Pick Modes中的设定值
//Define : Pick Modes
#define PICKMODE_EP     1
#define PICKMODE_PC     2
//Index of Array : Improve Mask Additions
#define MASKPARAM_NUMBER     3
#define MASKPARAM_SCATTER    0
#define MASKPARAM_ERODE      1
#define MASKPARAM_CORNER     2

enum ratio_name{
    A, B, C, D,     //串扰系数
    G, K,           //校准因子
    E               //消光系数比
};
enum channel_name{AA, DA, DD};
enum result_name{avg1, std1, avg2, std2};

using namespace std;

class FRETCalculator
{
private:
    //串扰系数与校正因子
    double ratio[7];
    //由于这是计算的基本单元，因此只针对单个样本
    QString view_path, batch_path;
    //记录各个通道文件是否存在的标志位
    bool flag_exist[3];
    //各通道背景值
    double gray_bg[3];  //gray of background
    //各通道原图片
    cv::Mat mat_src[3]; //matrix of the source image
    //掩膜模板
    cv::Mat mask_sc[3]; //mask of each single channel
    cv::Mat mask_final, mask_org;
    //计算的结果矩阵
    cv::Mat mat_result[2];
    //均值与方差
    double result[4];
    //附加功能开关标志位
    bool additional_function[ADDIFUNC_NUMBER];
    //设置优化模板的参数
    int mask_param[MASKPARAM_NUMBER];
    //取点参数
    int pick_param[PICKPARAM_NUMBER];

    //线性回归数据存储
    vector<pair<ld, ld>> point;
    ld slope;
    ld w, b;
    ld approach;

    //线性回归函数库
    void loadData(QString csv_path, ld rmin, ld rmax);
    void showImage(int channel_name);
    ld gradientDescend(ld w, ld lr);
    vector<ld> gradientDescend(ld w, ld b, ld lr);
    ld computeError(ld w, ld b);
    ld computeError(ld w);

/*功能函数********/
    //通过直方图寻找背景值，适用信噪比大的图片
    double calcBackGround(cv::Mat mat);
    //通过背景区域模板计算背景值，适用于信噪比小的图片，但要借助信噪比大的通道
    double calcBackGround(cv::Mat mat, cv::Mat mask);
    //用扣除背景后的矩阵计算模板
    bool calcMask(cv::Mat input_mat, cv::Mat output_mask);
    //进行多通道模板合成
    bool andMat(cv::Mat, cv::Mat);

    //计算均值
    double calcAvg(cv::Mat mat, cv::Mat mask);
    //计算标准差
    double calcStd(cv::Mat mat, cv::Mat mask, double mean);

    //优化细胞模板：消除散点
    void improveMaskScatter(int kernel_size);
    //优化细胞模板：腐蚀图像
    void improveMaskErode(int kernel_size);
    //优化细胞模板：消除四边
    void improveMaskCorner(int radius);

    //读取某一个视野的三通道图片
    bool loadImage(QString viewpath);
    //输出与数据可视化单元
    bool saveMask(QString path);

public:
    FRETCalculator();

    //计算单视野三通道数据，是最小的计算单元
    bool calcData(int calc_mode);
    //读取某一批样本
    bool calcBatchData(QString path, int calc_mode);

    //设置相关系数
    bool setRatio(int ratio, double value);
    //设置附加功能
    bool setAdditionalFunction(int func_name, bool state);
    //设置取点时的参数
    bool setPickParam(int pickparam_name, int value);
    //设置优化模板的参数
    bool setMaskParam(int maskparam_name, int value);

    //获得结果
    double getResult(int);
    //空间随机均匀取点
    bool getERPairByEquallyPicking(int num);
    //像素合并平均取点
    bool getERPairByPixelCombining(int num);
    //像素合并：按R或者1/R进行区间划分，取各区间内均值
    bool binData(QString inputPath, QString outputPath, ld min, ld max, int num);
    //线性回归函数

    void linearRegression(QString, int epochs, ld lr, ld min, ld max);  //针对的是多质粒拟合时使用的算法
    void calcSlope(QString path, int epochs, ld lr, ld rmin, ld rmax);  //针对的是计算截距为0的直线的斜率
    void calcApproach(QString path, ld rmin, ld rmax);                  //计算一个x范围内的y的平均值
    ld getStoichiometry();
    ld getSlope();
    ld getApproach();
    //统计模板像素点
    int countMaskPixels();


};

#endif // FRETCALCULATOR_H
