#pragma once

#include <QMainWindow>
#include <QTabWidget>
#include <QPushButton>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QThread>
#include "Worker.h"
#include "simulation/NTL_sim.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private:
    // UI Elements
    QTabWidget* m_plotTabs;
    QPushButton* m_startButton;
    QTextEdit* m_console;
    QWidget* m_centralWidget;
    QVBoxLayout* m_mainLayout;

    // Threading
    QThread* m_workerThread;

    // Helper function to build the NTL_opt_setup from the (future) UI
    NTL::NTL_opt_setup getSetupFromUI();

private slots:
    // Slot to handle the "Start" button click
    void onStartOptimization();

    // Slot to receive the NTL object when the worker is finished
    void onOptimizationFinished(const NTL::NTL& ntl);

    // Slot to receive console messages from the worker
    void onConsoleMessage(ConsoleMessage message);
};