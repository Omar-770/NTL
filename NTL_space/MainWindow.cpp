#include "MainWindow.h"
#include <QChartView> // To get the chart view from the window

// Make sure NTL type is known to Qt's signal/slot system
// This MUST be done before connecting signals
Q_DECLARE_METATYPE(NTL::NTL);

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    // Register the types
    qRegisterMetaType<NTL::NTL>("NTL::NTL");
    qRegisterMetaType<ConsoleMessage>("ConsoleMessage");

    // 1. Set up the UI
    setWindowTitle("NTL Optimizer GUI");
    m_centralWidget = new QWidget(this);
    m_mainLayout = new QVBoxLayout(m_centralWidget);

    m_startButton = new QPushButton("Start Optimization", this);
    m_plotTabs = new QTabWidget(this);
    m_console = new QTextEdit(this);
    m_console->setReadOnly(true);
    m_console->setFontFamily("Courier");

    m_mainLayout->addWidget(m_startButton);
    m_mainLayout->addWidget(m_plotTabs, 1); // Give tabs more stretch
    m_mainLayout->addWidget(m_console, 0); // Give console less stretch

    setCentralWidget(m_centralWidget);
    resize(1200, 800);

    // 2. Connect the "Start" button
    connect(m_startButton, &QPushButton::clicked, this, &MainWindow::onStartOptimization);
}

MainWindow::~MainWindow()
{
    // Clean up the thread
    if (m_workerThread) {
        m_workerThread->quit();
        m_workerThread->wait();
    }
}

NTL::NTL_opt_setup MainWindow::getSetupFromUI()
{
    // --- THIS IS WHERE YOU WOULD READ FROM UI WIDGETS ---
    // For now, we'll just hard-code the setup from your main.cpp

    NTL::NTL_opt_setup setup;
    setup.Z0 = 50; setup.Zs = 50; setup.Zl = 100; setup.er = 4.6; setup.d = 80e-3;
    setup.freqs = { 0.5e9, 1.3e9, 2e9 };
    setup.N = 10; setup.K = 50; setup.lb = std::vector<double>(setup.N, -1);
    setup.ub = std::vector<double>(setup.N, 1);
    setup.toll_bounds = std::vector<double>(2, 1e-6);
    setup.toll_z = std::vector<double>(setup.K, 1e-6);
    setup.GBL_MAX = 350e3; setup.LCL_MAX = 50e3; setup.accepted_error = 6e-3;
    setup.m_Z_min = 0.1 * setup.Z0; setup.m_Z_max = 2.6 * setup.Z0;
    setup.max_attempts = 1;

    return setup;
}

void MainWindow::onStartOptimization()
{
    m_startButton->setEnabled(false);
    m_console->clear();
    m_plotTabs->clear();
    m_console->append("Starting...");

    // 1. Get the setup
    NTL::NTL_opt_setup setup = getSetupFromUI();

    // 2. Create the worker and thread
    m_workerThread = new QThread(this);
    Worker* worker = new Worker(setup);
    worker->moveToThread(m_workerThread);

    // 3. Connect signals and slots
    //    - Start the worker when the thread starts
    connect(m_workerThread, &QThread::started, worker, &Worker::runOptimization);

    //    - Receive the final NTL object
    connect(worker, &Worker::optimizationFinished, this, &MainWindow::onOptimizationFinished);

    //    - Receive console messages
    connect(worker, &Worker::consoleMessage, this, &MainWindow::onConsoleMessage);

    //    - Clean up when the worker is done
    connect(worker, &Worker::optimizationFinished, m_workerThread, &QThread::quit);
    connect(worker, &Worker::optimizationFinished, worker, &Worker::deleteLater);
    connect(m_workerThread, &QThread::finished, m_workerThread, &QThread::deleteLater);
    connect(m_workerThread, &QThread::finished, [this]() {
        m_startButton->setEnabled(true);
        m_console->append("...Finished.");
        });

    // 4. Start the thread
    m_workerThread->start();
}

void MainWindow::onOptimizationFinished(const NTL::NTL& ntl)
{
    m_console->append("Optimization finished. Generating plots...");

    // --- THIS IS YOUR SIMULATION LOGIC ---
    // It runs on the main thread, but it's very fast.

    NTL::NTL_sim sim(ntl);
    sim.set_f_sweep(0.1e9, 3e9); // Wider sweep for plotting

    // Run all the simulations
    sim.z_profile();
    sim.w_h_profile();

    // Get Zs/Zl from the UI setup
    NTL::NTL_opt_setup setup = getSetupFromUI();
    sim.s_matrix(setup.Zs, setup.Zl);

    // Get the windows from the sim
    std::vector<QMainWindow*> windows = sim.get_windows();

    // Now, add them to our tab widget
// In MainWindow.cpp, inside onOptimizationFinished
    for (QMainWindow* window : windows)
    {
        if (!window) continue;

        QChartView* chartView = qobject_cast<QChartView*>(window->centralWidget());
        if (!chartView) {
            delete window;
            continue;
        }

        // "Steal" the chart view
        chartView->setParent(m_plotTabs);

        // --- FIX ---
        // Get the title from the chart, not the window
        QString title = chartView->chart()->title();
        m_plotTabs->addTab(chartView, title);
        // --- END FIX ---

        // We don't need the temporary QMainWindow anymore
        delete window;
    }
}

void MainWindow::onConsoleMessage(ConsoleMessage message)
{
    if (message.isError) {
        m_console->setTextColor(Qt::red);
    }
    else {
        m_console->setTextColor(Qt::black);
    }
    m_console->append(message.msg);
}