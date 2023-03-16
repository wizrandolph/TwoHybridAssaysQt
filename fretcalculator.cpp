#include "fretcalculator.h"


FRETCalculator::FRETCalculator()
{

    //初始化参数
    for (int i = 0; i < MASKPARAM_NUMBER; ++ i)
    {
        mask_param[i] = -1;
    }
}

/**
 * @brief FRETCalculator::loadImage
 * @details
 * Load images from the files path
 * If all three channel images are detected and loaded return true
 * else return false
 * @param path
 * @return completed or not
 */
bool FRETCalculator::loadImage(QString path)
{
    this->view_path = path;

    flag_exist[0] = flag_exist[1] = flag_exist[2] = false;
    QDir dir(path);
    QFileInfoList fileinfolist = dir.entryInfoList(QStringList() << "*.tif" << "*.TIF");

    QString keywords[6] = {"AA", "aa", "DA", "da", "DD", "dd"};

    for ( int i = 0; i < fileinfolist.length(); ++ i )
    {
        QFileInfo fileinfo = fileinfolist.at(i);
        QString filename = fileinfo.fileName();

        //根据文件名字符串检索自动读取图片文件
        for (int j = 0; j < 6; ++ j)
        {
            //qDebug() << filename.indexOf(keywords[j]);
            if (filename.indexOf(keywords[j]) != -1)
            {
                mat_src[j / 2] = cv::imread(fileinfo.absoluteFilePath().toStdString(), -1);
                if (!mat_src[j / 2].empty())
                {
                    flag_exist[j / 2] = true;
                }
                break;
            }
        }
    }

    if (DEBUG_ALL && DEBUG_LOADIMAGE)
    {
        for (int i = 0; i < 3; ++ i)
        {
            if (flag_exist[i])
                qDebug() << "Channel[" + QString::number(i) + "] Image Loaded";
                qDebug() << mat_src[i].rows << "," << mat_src[i].cols;
        }
    }

    return flag_exist[0] && flag_exist[1] && flag_exist[2];
}

/**
 * @brief FRETCalculator::calcBackGround
 * Calculate the background gray value by histogram.
 * Find the gray value that most pixels' gray values are.
 * @param mat
 * @return the background gray value
 */
double FRETCalculator::calcBackGround(cv::Mat mat)
{
    using namespace cv;

    //计算直方图
    Mat hist;   //直方图
    float range[2] = {0, 500};  //统计灰度值范围
    const float* ranges[1] = {range};   //格式需要，指针的指针
    const int bins[1] = {500};  //宽度，即直方图的柱数，即横轴的分布
    calcHist(&mat, 1, 0, Mat(), hist, 1, bins, ranges);  //计算直方图

    //直方图求峰值
    double minValue, maxValue;
    Point minIdx, maxIdx;
    minMaxLoc(hist, &minValue, &maxValue, &minIdx, &maxIdx);

    double bg = maxIdx.y;   //直方图峰值
    return bg;
}

/**
 * @brief FRETCalculator::calcBackGround
 * calc the mean
 * @param mat
 * @param mask
 * @return the background gray value
 */
double FRETCalculator::calcBackGround(cv::Mat mat, cv::Mat mask)
{
    using namespace cv;

    int num = 0;
    double sum = 0;

    mask.convertTo(mask, CV_8U);

    for (int r = 0; r < mask.rows; r++)
    {
        for (int c = 0; c < mask.cols; c++)
        {
            if(mask.at<uchar>(r, c) == MAXGRAY_C8)
            {
                ++ num;
                sum += (mat.at<double>(r, c));
            }
        }
    }

    return sum / num;
}

/**
 * @brief FRETCalculator::calcMask
 * @param input_mat
 * @param output_mask
 * @return bool : done or not
 */
bool FRETCalculator::calcMask(cv::Mat input_mat, cv::Mat output_mask)
{
    using namespace cv;

    int type = input_mat.type();
    if (type == CV_64F)
    {
        for (int r = 0; r < input_mat.rows; r++)
            for (int c = 0; c < input_mat.cols; c++)
            {
                if (input_mat.at<double>(r, c) <= 0)
                {
                    output_mask.at<uchar>(r, c) = 0;
                }
                else
                {
                    output_mask.at<uchar>(r, c) = MAXGRAY_C8;
                }

            }
        return true;
    }

    else if (type == CV_16U)
    {
        for (int r = 0; r < input_mat.rows; r++)
            for (int c = 0; c < input_mat.cols; c++)
            {
                if (input_mat.at<ushort>(r, c) == 0)
                {
                    output_mask.at<uchar>(r, c) = 0;
                }
                else
                {
                    output_mask.at<uchar>(r, c) = MAXGRAY_C8;
                }
            }
        return true;
    }
    else
    {
        qDebug() << "Mat type not supported";
        return false;
    }
}

/**
 * @brief FRETCalculator::andMat
 * @param mat
 * @param mask
 * @return bool : Done or not
 */
bool FRETCalculator::andMat(cv::Mat mat, cv::Mat mask)
{
    //检查mat与mask大小是否匹配
    if (mat.rows != mask.rows || mat.cols != mask.cols)
    {
        qDebug() << "Mask size unmatch";
        return false;
    }

    if (mat.type() == mask.type())
    {
        //将mask与mat进行合成
        mat = mask.mul(mat);
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * @brief FRETCalculator::calcAvg
 * @param mat
 * @param mask
 * @return double : average value
 */
double FRETCalculator::calcAvg(cv::Mat mat, cv::Mat mask)
{
    int num = 0;
    double sum = 0;

    for (int r = 0; r < mask.rows; r ++)
    {
        for (int c = 0; c < mask.cols; c ++)
        {
            if(mask.at<uchar>(r, c) == MAXGRAY_C8)
            {
                num ++;
                sum = sum + mat.at<double>(r, c);
            }
        }
    }
    if (num == 0) return -1;

    double avg = sum / num;
    return avg;
}

/**
 * @brief FRETCalculator::calcStd
 * @param mat
 * @param mask
 * @param mean
 * @return double \ standard deviation
 */
double FRETCalculator::calcStd(cv::Mat mat, cv::Mat mask, double mean)
{
    int num = 0;
    double sum = 0;

    for (int r = 0; r < mask.rows; r++)
    {
        for (int c = 0; c < mask.cols; c++)
        {
            if( mask.at<uchar>(r, c) == MAXGRAY_C8 )
            {
                num ++;
                sum = sum + pow((mat.at<double>(r, c) - mean), 2);
            }
        }
    }

    if (num == 0) return -1;

    double std = sqrt(sum / num);
    return std;
}

/**
 * @brief FRETCalculator::calcData
 * Calculate FRET data from the images loaded to FRET values
 * Use calc_mode to change the mode of calculation:
 * CALCMODE_AB
 * @param calc_mode
 * @return bool : done or not
 */
bool FRETCalculator::calcData(int calc_mode)
{
    using namespace cv;

    //拷贝数值，同时进行类型转换
    Mat mat_copy[3];
    for (int i = 0; i < 3; ++ i)
    mat_src[i].convertTo(mat_copy[i], CV_64F);

/*计算背景值**************************************/
    if (calc_mode == CALCMODE_AB)
    {
        Mat mask_bg;   //背景位置为MAXGRAY，非背景位置是0
        gray_bg[AA] = calcBackGround(mat_src[AA]);
        threshold(mat_src[AA], mask_bg, gray_bg[AA] * 2, MAXGRAY_C8, THRESH_BINARY_INV);
        gray_bg[DA] = calcBackGround(mat_copy[DA], mask_bg);
        gray_bg[DD] = calcBackGround(mat_copy[DD], mask_bg);
    }
    else if (calc_mode == CALCMODE_CD)
    {
        Mat mask_bg;   //背景位置为MAXGRAY，非背景位置是0
        gray_bg[DD] = calcBackGround(mat_src[DD]);
        threshold(mat_src[DD], mask_bg, gray_bg[DD] * 2, MAXGRAY_C8, THRESH_BINARY_INV);
        gray_bg[AA] = calcBackGround(mat_copy[AA], mask_bg);
        gray_bg[DA] = calcBackGround(mat_copy[DA], mask_bg);
    }
    else
    {
        if (calc_mode == CALCMODE_EFRET || calc_mode == CALCMODE_GK_AG || calc_mode == CALCMODE_GK_MF || calc_mode == CALCMODE_3FRET)
        {
            gray_bg[AA] = calcBackGround(mat_src[AA]);
            gray_bg[DA] = calcBackGround(mat_src[DA]);
            gray_bg[DD] = calcBackGround(mat_src[DD]);
        }
        else
        {
            //出现未定义的计算模式，报错
            return false;
        }
    }
/*计算掩膜Mask******************************************/
    double multiple[3];
    if (calc_mode == CALCMODE_AB)
    {
        multiple[AA] = 3;
        multiple[DA] = 1;
        multiple[DD] = 1;
    }
    else if (calc_mode == CALCMODE_CD)
    {
        multiple[AA] = 1;
        multiple[DA] = 1;
        multiple[DD] = 3;
    }
    else if (calc_mode == CALCMODE_EFRET || calc_mode == CALCMODE_GK_AG || calc_mode == CALCMODE_GK_MF || calc_mode == CALCMODE_3FRET)
    {
        multiple[AA] = multiple[DA] = multiple[DD] = 4;
    }
    else
    {
        return false;
    }

    //计算各个通道单独的掩膜模板
    threshold(mat_src[AA], mask_sc[AA], gray_bg[AA] * multiple[AA], MAXGRAY_C8, THRESH_BINARY);
    threshold(mat_src[DA], mask_sc[DA], gray_bg[DA] * multiple[DA], MAXGRAY_C8, THRESH_BINARY);
    threshold(mat_src[DD], mask_sc[DD], gray_bg[DD] * multiple[DD], MAXGRAY_C8, THRESH_BINARY);

    for (int i = 0; i < 3; ++ i)
    {
        mask_sc[i].convertTo(mask_sc[i], CV_8U);
    }

    //合成掩膜
    bitwise_and(mask_sc[AA], mask_sc[DA], mask_org);
    bitwise_and(mask_org, mask_sc[DD], mask_org);

    //优化掩膜
    mask_org.convertTo(mask_final, CV_8U);
    //imshow("CV_8U", mask_final);
    mask_org.convertTo(mask_org, CV_8U);

    improveMaskScatter(mask_param[MASKPARAM_SCATTER]);
    improveMaskErode(mask_param[MASKPARAM_ERODE]);
    improveMaskCorner(mask_param[MASKPARAM_CORNER]);

/*扣除背景，得到校正的数据******************************/
    Mat mat_corr[3];
    for (int i = 0; i < 3; ++ i)
    {
        mat_corr[i] = mat_copy[i] - gray_bg[i];
        //qDebug() << gray_bg[i];
    }

/*应用FRET公式计算************************************/
    if (calc_mode == CALCMODE_EFRET)
    {
        Mat mat_fc = mat_corr[DA]
                     - ratio[A] * (mat_corr[AA] - ratio[C] * mat_corr[DD])
                     - ratio[D] * (mat_corr[DD] - ratio[B] * mat_corr[AA]);
        Mat mat_denominater = mat_fc + ratio[G] * mat_corr[DD];

        mat_result[0] = mat_fc / mat_denominater;           //Ed

        Mat mat_numerator = ratio[K] * mat_corr[AA];
        mat_denominater = mat_fc / ratio[G] + mat_corr[DD];
        mat_result[1] = mat_numerator / mat_denominater;    //Rc

        //去除负值
        Mat mask_pos(mat_result[0].rows, mat_result[0].cols, CV_8U);
        calcMask(mat_result[0], mask_pos);
        bitwise_and(mask_final, mask_pos, mask_final);


    }
    else if (calc_mode == CALCMODE_3FRET)
    {
        //3Cube-FRET计算模式
        Mat mat_fc = mat_corr[DA]
                     - ratio[A] * (mat_corr[AA] - ratio[C] * mat_corr[DD])
                     - ratio[D] * (mat_corr[DD] - ratio[B] * mat_corr[AA]);

        mat_result[0] = mat_fc * ratio[E] / mat_corr[AA] / ratio[A];   //Ea

        Mat mat_denominater = mat_fc + ratio[G] * mat_corr[DD];
        Mat mat_numerator = ratio[K] * mat_corr[AA];
        mat_denominater = mat_fc / ratio[G] + mat_corr[DD];
        mat_result[1] =  mat_denominater / mat_numerator;    //1/Rc

        //去除负值
        Mat mask_pos(mask_final.rows, mask_final.cols, CV_8U);
        calcMask(mat_result[0], mask_pos);
        bitwise_and(mask_final, mask_pos, mask_final);
        //qDebug() << countMaskPixels();

    }
    else if (calc_mode == CALCMODE_CD)
    {
        mat_result[0] = mat_corr[DA] / mat_corr[DD];        //d
        mat_result[1] = mat_corr[AA] / mat_corr[DD];        //c
    }
    else if (calc_mode == CALCMODE_AB)
    {
        mat_result[0] = mat_corr[DA] / mat_corr[AA];        //a
        mat_result[1] = mat_corr[DD] / mat_corr[AA];        //b
    }
    else if (calc_mode == CALCMODE_GK_MF)
    {
        Mat mat_fc = mat_corr[DA]
                     - ratio[A] * (mat_corr[AA] - ratio[C] * mat_corr[DD])
                     - ratio[D] * (mat_corr[DD] - ratio[B] * mat_corr[AA]);

        mat_result[0] = (mat_corr[DD] / mat_corr[AA]) / ratio[A];   //横轴x
        mat_result[1] = (mat_fc / mat_corr[AA]) / ratio[A]; //纵轴y

        //去除负值
        Mat mask_pos(mask_final.rows, mask_final.cols, CV_8U);
        calcMask(mat_result[1], mask_pos);
        bitwise_and(mask_final, mask_pos, mask_final);
    }
    else
    {
        return false;
    }

/*计算均值***************************************************/
    result[avg1] = calcAvg(mat_result[0], mask_final);
    result[std1] = calcStd(mat_result[0], mask_final, result[avg1]);
    result[avg2] = calcAvg(mat_result[1], mask_final);
    result[std2] = calcStd(mat_result[1], mask_final, result[avg2]);

/*双杂交分析取点******************************************************/
    if (additional_function[ADDIFUNC_GENERATECSV])
    {
        if (pick_param[PICKPARAM_PICKMODE] == PICKMODE_EP)
        {
            getERPairByEquallyPicking(pick_param[PICKPARAM_POINTNUMBER]);
        }
        else if (pick_param[PICKPARAM_PICKMODE == PICKMODE_PC])
        {
            getERPairByPixelCombining(pick_param[PICKPARAM_POINTNUMBER]);
        }
    }

    //返回正常的结果
    return true;
}

/**
 * @name saveMask
 * @author wiz randolph
 * @brief save mask matrix to local file
 * @param path
 * @return success or not
 */
bool FRETCalculator::saveMask(QString path)
{
    using namespace cv;
    if (imwrite(path.toStdString(), mask_final)) return true;
    return false;
}

/**
 * @brief FRETCalculator::setRatio
 * @param value
 * @param ratio_name
 * @return completed or not
 */
bool FRETCalculator::setRatio(int ratio_name, double value)
{
    ratio[ratio_name] = value;
    return true;
}

/**
 * @brief FRETCalculator::showImage
 * @param channel_name
 */
void FRETCalculator::showImage(int channel_name)
{
    using namespace cv;

    Mat mat_8U = mat_src[channel_name] / 256;
    mat_src[channel_name].convertTo(mat_8U, CV_8U);
    imshow("img", mat_8U);
}

/**
 * @brief FRETCalculator::improveMaskScatter
 * @param kernel_size
 */
void FRETCalculator::improveMaskScatter(int kernel_size)
{
    if (kernel_size < 0) return;

    int pxnum = (kernel_size * 2 + 1) * (kernel_size * 2 + 1);
    cv::Mat mask_copy;
    mask_final.copyTo(mask_copy);

    for (int r = kernel_size; r < mask_final.rows - kernel_size - 1; ++ r)
        for (int c = kernel_size; c < mask_final.cols - kernel_size - 1; ++ c)
        {
            int count = 0;
            for (int i = 0 - kernel_size; i <= kernel_size; ++ i)
                for (int j = 0 - kernel_size; j <= kernel_size; ++ j)
                {
                    if (mask_final.at<uchar>(r + i, c + j) == MAXGRAY_C8)
                        ++ count;
                }
            mask_copy.at<uchar>(r, c) = count < pxnum * 0.3 ? 0 : MAXGRAY_C8;
        }
    mask_final = mask_copy;

}

/**
 * @brief FRETCalculator::improveMaskErode
 * @param kernel_size
 */
void FRETCalculator::improveMaskErode(int kernel_size)
{
    if (kernel_size < 0) return;
    using namespace cv;
    int width = kernel_size * 2 + 1;
    Mat element = getStructuringElement(MORPH_RECT, Size(width, width));
    Mat mask_src, mask_dst;

    mask_src = mask_final;
    mask_src.convertTo(mask_src, CV_8U);
    erode(mask_src, mask_dst, element);
    mask_final = mask_dst;
}
/**
 * @brief FRETCalculator::improveMaskCorner
 * @param radius : The circle radius
 */
void FRETCalculator::improveMaskCorner(int radius)
{
    using namespace cv;

    //如果半径为负值，相当于该功能没有开启
    if (radius < 0) return;
    if (radius * 2 > mask_final.cols)
    {
        qDebug() << "window size out of range";
        return;
    }
    int center = mask_final.rows / 2;
    for (int r = 0; r < mask_final.rows; ++ r)
    {
        for (int c = 0; c < mask_final.cols; ++ c)
        {
            if (abs(r - center) * abs(c - center) > radius * radius)
            {
                mask_final.at<uchar>(r, c) = 0;
            }
        }
    }
}
/**
 * @brief FRETCalculator::getResult
 * @param result_name
 * @return
 */
double FRETCalculator::getResult(int result_name)
{
    return result[result_name];
}

/**
 * @brief FRETCalculator::calcBatchData
 * @details
 * This function will detect all possible view paths under the batch path inputed.
 * Will generate conclusion csv file.
 * @param path : The path of the data batch
 * @param calc_mode : See the calculation modes in #define
 * @return bool : completed or not
 */
bool FRETCalculator::calcBatchData(QString path, int calc_mode)
{
    //检查文件目录是否存在
    QDir sampledir(path);
    if (!sampledir.exists()) return false;

    if (DEBUG_ALL && DEBUG_BATCH)
    {
        qDebug() << path;
    }

    this->batch_path = path;
    QFile pointfile(batch_path + "/all.csv");
    if (pointfile.exists()) pointfile.remove();

    //打开csv文件
    QFile file(batch_path + "/CalcResult.csv");
    if (file.exists()) file.remove();
    file.open(QIODevice::Append);
    QTextStream out(&file);

    //对目录下的每个子文件夹进行处理
    QFileInfoList folderlist = sampledir.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (int i = 0; i < folderlist.size(); ++ i)
    {
        if (DEBUG_ALL && DEBUG_VIEW) qDebug() << folderlist.at(i).absoluteFilePath();
        if (loadImage(folderlist.at(i).absoluteFilePath()))
        {
            if (calcData(calc_mode))
            {
                if (DEBUG_ALL && DEBUG_VIEW) qDebug() << "Data calculation completed";
                if (additional_function[ADDIFUNC_SAVEMASK]) saveMask(folderlist.at(i).absoluteFilePath() + "/mask.tif");
                out << result[avg1] << ","
                    << result[std1] << ","
                    << result[avg2] << ","
                    << result[std2] << "\n";
            }
        }
        else
        {
            if (DEBUG_ALL && DEBUG_VIEW) qDebug() << "Data imcomplete";
        }
    }

    file.close();

    return true;
}
/**
 * @brief FRETCalculator::getERPairByEquallyPicking
 * @details
 * This function should be called after you done the calculation.
 * Expect the final mask, result mat are all ready.
 * @param num : point numbers expected
 * @return bool : enough points have been picked or not
 */
bool FRETCalculator::getERPairByEquallyPicking(int num)
{
    if (num == -1) num = countMaskPixels();

    if (DEBUG_ALL && DEBUG_PICKPOINT && DEBUG_RUNDETAIL)
    {
        qDebug() << "Starting picking Ed-Rc points...";
    }
    int count = 0;
    int idx[16][16];
    memset(idx, 0, sizeof(idx));

    QFile file(batch_path + "/all.csv");
    file.open(QIODevice::Append);
    QTextStream out(&file);

    cv::Mat mask;
    mask_final.convertTo(mask, CV_8U);

    cv::imwrite((view_path + "/maskpick.tif").toStdString(), mask);

    bool flag = true;
    while (count < num && flag)
    {
        flag = false;
        for (int x = 0; x < 16; x ++)
        {
            for (int y = 0; y < 16; y ++)
            {
                while(idx[x][y] < 128 * 128)
                {
                    int rr, rc;
                    rr = x * 128 + idx[x][y] / 128;
                    rc = y * 128 + idx[x][y] % 128;
                    idx[x][y] ++;
                    if (mask.at<uchar>(rr, rc) == MAXGRAY_C8)
                    {
                        flag = true;
                        out << mat_result[1].at<double>(rr, rc)
                            << "," << mat_result[0].at<double>(rr, rc) << "\n";
                        count ++;
                        break;
                    }

                }
            }
        }
    }

    file.close();

    if (DEBUG_ALL && DEBUG_PICKPOINT && DEBUG_RUNRESULT)
    {
        qDebug() << count << "points are picked.";
    }

    if (count < num)
        return false;
    else
        return true;
}

/**
 * @brief FRETCalculator::getERPairByPixelCombining
 * @param num : pixel numbers expected
 * @return bool : enough point are picked or not
 */
bool FRETCalculator::getERPairByPixelCombining(int num)
{
    if (num == -1) num = countMaskPixels();
    if(DEBUG_ALL && DEBUG_PICKPOINT && DEBUG_RUNDETAIL)
    {
        qDebug() << "Starting picking Ed-Rc points...";
    }

    int count = 0;
    QFile file(batch_path + "/all.csv");
    file.open(QIODevice::Append);
    QTextStream out(&file);

    cv::Mat mask;
    mask_final.convertTo(mask, CV_8U);

    for (int x = 0; x < 1024; x ++)
    {
        for (int y = 0; y < 1024; y ++)
        {
            bool flag = true;
            double rc_sum = 0;
            double ed_sum = 0;
            for (int i = 0; i < 4; ++ i)
            {
                int rr, rc;
                rr = x * 2 + i / 2;
                rc = y * 2 + i % 2;
                if (mask.at<uchar>(rr, rc) == MAXGRAY_C8)
                {
                    rc_sum += mat_result[1].at<double>(rr, rc);
                    ed_sum += mat_result[0].at<double>(rr, rc);
                }
                else
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                out << rc_sum / 4 << "," << ed_sum / 4 << "\n";
                count ++;
            }

        }
    }
    file.close();

    if (DEBUG_ALL && DEBUG_PICKPOINT && DEBUG_RUNRESULT)
    {
        qDebug() << count << "points are picked";
    }
    if (count >= num) return true;
    else return false;
}
/**
 * @brief FRETCalculator::setAdditionalFunction
 * @param func_name
 * @param state
 * @return bool : done or not
 */
bool FRETCalculator::setAdditionalFunction(int func_name, bool state)
{
    QString strFunc;

    if (func_name == ADDIFUNC_GENERATECSV)
    {
       strFunc = "[GENERATE ER POINT CSV]";
    }
    else if (func_name == ADDIFUNC_SAVEMASK)
    {
        strFunc = "[SAVE FINAL MASK]";
    }

    if (DEBUG_ALL && DEBUG_RUNRESULT)
    qDebug() << "Set Additional Function" + strFunc;
    if (func_name < ADDIFUNC_NUMBER)
    {
        if (DEBUG_ALL && DEBUG_RUNDETAIL)
        {
            qDebug() << "Completed";
        }
        additional_function[func_name] = state;
    }
    else
    {
        return false;
    }
    return true;
}
/**
 * @brief FRETCalculator::setPickParam
 * @param pickparam_name
 * @param value
 * @return bool : done or not
 */
bool FRETCalculator::setPickParam(int pickparam_name, int value)
{
    if (pickparam_name >= PICKPARAM_NUMBER)
        return false;
    pick_param[pickparam_name] = value;
    return true;
}
bool FRETCalculator::setMaskParam(int maskparam_name, int value)
{
    if (maskparam_name >= PICKPARAM_NUMBER)
        return false;
    pick_param[maskparam_name] = value;
    return true;
}
/**
 * @brief FRETCalculator::loadData
 * @param csvpath
 * @param rmin
 * @param rmax
 */
void FRETCalculator::loadData(QString csvpath, ld rmin, ld rmax)
{
    //清空vector，但是clear()清除后并没有在内存上清空
    point.clear();


    QFile file(csvpath);
    if (!file.open(QIODevice::ReadOnly))
    {
        qDebug() << "Csv file missed";
        return;
    }

    QTextStream * read = new QTextStream(&file);
    QStringList Data = read->readAll().split("\n", Qt::SkipEmptyParts);

    for (int i = 0; i < Data.count(); ++ i)
    {
        QStringList strLine = Data.at(i).split(",");
        ld x = (strLine.at(0)).toDouble();
        if (x < rmin || x > rmax) continue;
        ld y = (strLine.at(1)).toDouble();
        point.push_back(make_pair(x, y));
    }

}

/**
 * @brief FRETCalculator::binData
 * @param csv_path
 * @param bin_path
 * @param min
 * @param max
 * @param num
 * @return bool : completed or not
 */
bool FRETCalculator::binData(QString inputPath, QString outputPath, ld min, ld max, int num)
{
    QFile csvfile(inputPath);
    if (!csvfile.exists()) return false;
    csvfile.open(QIODevice::ReadOnly);
    ld step = (max - min) / num;
    QTextStream * read = new QTextStream(&csvfile);
    QStringList data = read->readAll().split("\n", Qt::SkipEmptyParts);
    vector<vector<ld> > vecE(num), vecR(num);
    for (int i = 0; i < data.count(); ++ i)
    {
        QStringList strLine = data.at(i).split(",");
        ld x = (strLine.at(0)).toDouble();
        if (x < min || x > max) continue;
        int index = (x - min) / step;
        ld y = (strLine.at(1)).toDouble();
        vecE[index].push_back(y);
        vecR[index].push_back(x);
    }

    QFile binfile(outputPath);
    if (binfile.exists()) binfile.remove();
    binfile.open(QIODevice::Append);
    QTextStream out(&binfile);
    for (int i = 0; i < num; ++ i)
    {
        if(!vecE[i].size()) continue;
        ld sum_y = accumulate(begin(vecE[i]), end(vecE[i]), 0.0);
        ld y = sum_y / (ld)vecE[i].size();
        ld sum_x = accumulate(begin(vecR[i]), end(vecR[i]), 0.0);
        ld x = sum_x / (ld)vecR[i].size();
        out << (double)x << "," << (double)y << "\n";
    }

    binfile.close();
    return true;
}

/**
 * @brief FRETCalculator::gradientDescend
 * @param w
 * @param lr
 * @return w updated
 */
ld FRETCalculator::gradientDescend(ld w, ld lr)
{
    ld dw = 0;
    int n = point.size();
    for (int i = 0; i < n; ++ i)
    {
        ld x = point[i].first;
        ld y = point[i].second;
        dw += -(2.0 / n) * x * (y - w * x);
    }
    return w - dw * lr;
}
vector<ld> FRETCalculator::gradientDescend(ld w, ld b, ld lr)
{
    ld dw = 0, db = 0;
    vector<ld> tmp;
    int n = point.size();
    for (int i = 0; i < n; ++ i)
    {
        ld x = point[i].first;
        ld y = point[i].second;
        db += -(2.0 / n) * (y - (w * x + b));
        dw += -(2.0 / n) * x * (y - (w * x + b));
    }
    tmp.push_back(w - (dw * lr));
    tmp.push_back(b - (db * lr));
    return tmp;
}
/**
 * @brief FRETCalculator::computeError
 * @param w
 * @return error
 */
ld FRETCalculator::computeError(ld w)
{
    ld loss = 0;
    int n = point.size();
    for (int i = 0; i < n; ++ i)
    {
        ld x = point[i].first;
        ld y = point[i].second;
        loss += (y - w * x) * (y - w * x);
    }
    return loss / (ld)n;
}
ld FRETCalculator::computeError(ld w, ld b)
{
    ld loss = 0;
    int n = point.size();
    for (int i = 0; i < n; ++ i)
    {
        ld x = point[i].first;
        ld y = point[i].second;
        loss += (y - (w * x + b)) * (y - (w * x + b));
    }
    return loss / (ld)n;
}
void FRETCalculator::linearRegression(QString path, int epochs, ld lr, ld min, ld max)
{
    loadData(path, min, max);
    w = -5, b = 20;
    for (int epoch = 0; epoch < epochs; ++ epoch)
    {
        vector<ld> t = gradientDescend(w, b, lr);
        w = t[0];
        b = t[1];
    }

    qDebug() << epochs << "iterations" << "w = " << (double)w << "b = " << (double)b;
}
/**
 * @brief FRETCalculator::linearRegression
 * @param path
 * @param epochs
 * @param lr
 * @param rmin
 * @param rmax
 */
void FRETCalculator::calcSlope(QString path, int epochs, ld lr, ld rmin, ld rmax)
{
    loadData(path, rmin, rmax);
    slope = 0.3;    //斜率初始值
    if (DEBUG_ALL && DEBUG_RUNDETAIL) qDebug() << "Starting gradient descent at slope = " << (double)slope << ", error = " << (double)computeError(slope);
    if (DEBUG_ALL && DEBUG_RUNDETAIL) qDebug() << "Running...";
    for (int epoch = 0; epoch < epochs; ++ epoch)
    {
        slope = gradientDescend(slope, lr);
    }
    if (DEBUG_ALL && DEBUG_RUNDETAIL) qDebug() << "Iteration " << epochs << ", error is " << (double)computeError(slope);
    if (DEBUG_ALL && DEBUG_RUNRESULT) qDebug() << "Slope = " << (double)slope;
}
/**
 * @brief FRETCalculator::calcApproach
 * @param path
 * @param rmin
 * @param rmax
 */
void FRETCalculator::calcApproach(QString path, ld rmin, ld rmax)
{
    loadData(path, rmin, rmax);
    if (DEBUG_ALL && DEBUG_RUNDETAIL) qDebug() << "Starting calculating the asymptote between " << (double)rmin << "and" << (double)rmax;
    if (DEBUG_ALL && DEBUG_RUNDETAIL) qDebug() << "Running...";
    ld sum = 0;
    for (uint i = 0; i < point.size(); ++ i)
    {
        sum += point[i].second;
    }
    approach = sum / point.size();
    if (DEBUG_ALL && DEBUG_RUNRESULT) qDebug() << "Approach = " << (double)approach;
}
/**
 * @brief FRETCalculator::stoichiometry
 * @return stoichiometry
 */
ld FRETCalculator::getStoichiometry()
{
    return slope / approach;
}
ld FRETCalculator::getApproach()
{
    return approach;
}
ld FRETCalculator::getSlope()
{
    return slope;
}
/**
 * @brief FRETCalculator::countMaskPixels
 * @return pixel number
 */
int FRETCalculator::countMaskPixels()
{
    int sum = 0;
    cv::Mat mask;
    mask_final.convertTo(mask, CV_8U);
    for (int r = 0; r < mask.rows; ++ r)
    {
        for (int c = 0; c < mask.cols; ++ c)
        {
            if (mask.at<uchar>(r, c) == MAXGRAY_C8)
                sum ++;
        }
    }
    return sum;
}
