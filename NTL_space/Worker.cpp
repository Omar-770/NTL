#include "Worker.h"
#include <QDebug> // For thread-safe console output (optional)

// We need to capture std::cout. This is an advanced trick.
// For simplicity, we'll just emit our own messages.
// In a real app, you'd redirect std::cout to emit a signal.

void Worker::runOptimization()
{
    // Register the types we're going to pass
    qRegisterMetaType<NTL::NTL>("NTL::NTL");
    qRegisterMetaType<ConsoleMessage>("ConsoleMessage");

    NTL::NTL ntl;
    try
    {
        // This is a placeholder for redirecting cout.
        // For now, we'll just emit our own messages.
        // NTL::console::active will print to the *debug* console, not the GUI.
        // To fix this, you would modify NTL_opt to take a callback.

        emit consoleMessage({ "Starting optimization...", false });

        NTL::NTL_opt opt(m_setup);

        // --- THIS IS THE LONG-RUNNING TASK ---
        ntl = opt.optimise(ntl, NTL::console::inactive); // Run silently
        // ------------------------------------

        emit consoleMessage({ "Optimization complete.", false });

        // Send the final result back to the main window
        emit optimizationFinished(ntl);
    }
    catch (const std::exception& e)
    {
        emit consoleMessage({ QString("Error: %1").arg(e.what()), true });
    }
    catch (...)
    {
        emit consoleMessage({ "An unknown critical error occurred.", true });
    }
}