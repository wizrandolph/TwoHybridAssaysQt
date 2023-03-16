#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    initTableHeader();

    calcBatchED();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initTableHeader()
{
    QStandardItemModel *table_model = new QStandardItemModel();
    table_model->setHorizontalHeaderItem(0, new QStandardItem(QObject::tr("视野名")));
    table_model->setHorizontalHeaderItem(1, new QStandardItem(QObject::tr("就绪")));
    table_model->setHorizontalHeaderItem(2, new QStandardItem(QObject::tr("处理")));
    ui->tableView->setModel(table_model);

    ui->tableView->setColumnWidth(0,300);
    ui->tableView->setColumnWidth(1,100);
    ui->tableView->setColumnWidth(2,100);
}

void MainWindow::calcBatchEA()
{
    //设好参数
    calculator->setRatio(A, 0.186046857);
    calculator->setRatio(B, 0.002056633);
    calculator->setRatio(C, 0.001193742);
    calculator->setRatio(D, 0.781544552);
    calculator->setRatio(G, 4.6658);
    calculator->setRatio(K, 0.645755859);
    calculator->setRatio(E, 0.061747);

    //设置模板功能
    //calculator->setAdditionalFunction(SAVE_MASK, true);

    //设置取点参数
    calculator->setAdditionalFunction(ADDIFUNC_GENERATECSV, true);
    calculator->setPickParam(PICKPARAM_POINTNUMBER, -1);    //将取点数量设置为-1将会按照最大数量取点
    calculator->setPickParam(PICKPARAM_PICKMODE, PICKMODE_PC);  //设置取点模式

    //双杂交批处理
    //设置附加功能
    calculator->setAdditionalFunction(ADDIFUNC_SAVEMASK, true);

    //设置取点参数
    calculator->setPickParam(PICKPARAM_POINTNUMBER, -1);    //将取点数量设置为-1将会按照最大数量取点
    calculator->setPickParam(PICKPARAM_PICKMODE, PICKMODE_PC);  //设置取点模式

    QString batch_paths[2];
    batch_paths[0] = "D:/FretData/test/THA/hole";
    batch_paths[1] = "D:/FretData/test/THA/mhole";

    QString outpath = "D:/FretData/test/THA/thaResult.csv";
    QFile outfile(outpath);
    if (outfile.exists()) outfile.remove();
    outfile.open(QIODevice::Truncate | QIODevice::WriteOnly | QIODevice::Append);
    QTextStream ts(&outfile);

    ts << "Eamax," << "Edmax," << "nD/nA\n";

    for (int j = 0; j < 2; ++ j)
    {
        for (int i = 1; i <= 4; ++ i)
        {

            QString batch_path = batch_paths[j] + QString::number(i);
            ts << batch_path << "\n";
            //计算
            calculator->calcBatchData(batch_path, CALCMODE_3FRET);
            //源数据直接计算
            calculator->calcSlope(batch_path + "/all.csv", 100000, 0.0001, 0, 0.5);
            calculator->calcApproach(batch_path + "/all.csv", 2, 10);
            ts << (double)calculator->getSlope() << ",";
            ts << (double)calculator->getApproach() << ",";
            ts << (double)calculator->getStoichiometry() << "\n";
            //全区间合并
            calculator->binData(batch_path + "/all.csv", batch_path + "/allbin.csv", 0, 10, 200);
            calculator->calcSlope(batch_path + "/allbin.csv", 100000, 0.0001, 0, 0.5);
            calculator->calcApproach(batch_path + "/allbin.csv", 2, 10);
            ts << (double)calculator->getSlope() << ",";
            ts << (double)calculator->getApproach() << ",";
            ts << (double)calculator->getStoichiometry() << "\n";
        }
    }

    outfile.close();
}

void MainWindow::calcBatchED()
{
    //设好参数
    calculator->setRatio(A, 0.186046857);
    calculator->setRatio(B, 0.002056633);
    calculator->setRatio(C, 0.001193742);
    calculator->setRatio(D, 0.781544552);
    calculator->setRatio(G, 4.6658);
    calculator->setRatio(K, 0.645755859);
    calculator->setRatio(E, 0.061747);

    //设置模板功能
    //calculator->setAdditionalFunction(SAVE_MASK, true);
//    calculator->setMaskParam(MASKPARAM_SCATTER, 3);
//    calculator->setMaskParam(MASKPARAM_ERODE, 3);
//    calculator->setMaskParam(MASKPARAM_CORNER, 1000);

    //设置取点参数
    calculator->setAdditionalFunction(ADDIFUNC_GENERATECSV, true);
    calculator->setPickParam(PICKPARAM_POINTNUMBER, -1);    //将取点数量设置为-1将会按照最大数量取点
    calculator->setPickParam(PICKPARAM_PICKMODE, PICKMODE_PC);  //设置取点模式

    //设置附加功能
    calculator->setAdditionalFunction(ADDIFUNC_SAVEMASK, true);


    QString batch_paths[2];
    batch_paths[0] = "D:/FretData/test/THA/hole";
    batch_paths[1] = "D:/FretData/test/THA/mhole";

    QString outpath = "D:/FretData/test/THA/thaResult.csv";
    QFile outfile(outpath);
    if (outfile.exists()) outfile.remove();
    outfile.open(QIODevice::Truncate | QIODevice::WriteOnly | QIODevice::Append);
    QTextStream ts(&outfile);

    ts << "Eamax," << "Edmax," << "nD/nA\n";

    for (int j = 0; j < 2; ++ j)
    {
        for (int i = 1; i <= 4; ++ i)
        {

            QString batch_path = batch_paths[j] + QString::number(i);
            ts << batch_path << "\n";
            //计算
            calculator->calcBatchData(batch_path, CALCMODE_EFRET);
            //源数据直接计算
            calculator->calcSlope(batch_path + "/all.csv", 100000, 0.0001, 0, 0.5);
            calculator->calcApproach(batch_path + "/all.csv", 3, 10);
            ts << (double)calculator->getSlope() << ",";
            ts << (double)calculator->getApproach() << ",";
            ts << (double)calculator->getStoichiometry() << "\n";
            //全区间合并
            calculator->binData(batch_path + "/all.csv", batch_path + "/allbin.csv", 0, 10, 200);
            calculator->calcSlope(batch_path + "/allbin.csv", 100000, 0.0001, 0, 0.5);
            calculator->calcApproach(batch_path + "/allbin.csv", 3, 10);
            ts << (double)calculator->getSlope() << ",";
            ts << (double)calculator->getApproach() << ",";
            ts << (double)calculator->getStoichiometry() << "\n";
        }
    }

    outfile.close();
}
