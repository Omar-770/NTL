#pragma once

#include <QObject>
#include <QString>
#include "optimisation/NTL_opt.h"
#include "models/ntl.h"

// Define a simple struct to pass console messages
// This is necessary because std::cout can't be used in a thread
struct ConsoleMessage {
    QString msg;
    bool isError;
};
Q_DECLARE_METATYPE(ConsoleMessage); // Register for signals/slots

class Worker : public QObject
{
    Q_OBJECT

public:
    // The worker takes the setup struct when created
    explicit Worker(NTL::NTL_opt_setup setup, QObject* parent = nullptr)
        : QObject(parent), m_setup(setup) {
    }

signals:
    // Signal to send the final NTL object back to the main thread
    void optimizationFinished(const NTL::NTL& ntl);

    // Signal to send console output back
    void consoleMessage(ConsoleMessage message);

public slots:
    // This is the main function that will be run on the new thread
    void runOptimization();

private:
    NTL::NTL_opt_setup m_setup;
};