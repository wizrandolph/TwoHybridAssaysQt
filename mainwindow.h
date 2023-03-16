#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <fretcalculator.h>
#include <QStandardItemModel>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    FRETCalculator *calculator = new FRETCalculator();
    void calcBatchED();
    void calcBatchEA();
    void calcBatchGK();

/*GUI部分**************************/
    void initTableHeader();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
