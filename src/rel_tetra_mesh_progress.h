#ifndef REL_TETRA_MESH_PROGRESS_H
#define REL_TETRA_MESH_PROGRESS_H

#include <chrono>
#include <cstdint>

#include <rbl_logger.h>
#include <rbl_progress.h>

// Internal phase progress. Adaptive queues have no fixed denominator; keep
// those phases pulsing and throttle text updates to avoid flooding the GUI.
class RTetraMeshProgress
{
    using Clock = std::chrono::steady_clock;
    Clock::time_point lastUpdate;
    bool active = false;
    bool logBlockOpen = false;

    public:

        explicit RTetraMeshProgress(const QString &message, bool pulse = true)
        {
            initialize(message, pulse);
        }

        ~RTetraMeshProgress()
        {
            finalize();
        }

        RTetraMeshProgress(const RTetraMeshProgress &) = delete;
        RTetraMeshProgress &operator=(const RTetraMeshProgress &) = delete;

        void initialize(const QString &message, bool pulse = true)
        {
            finalize();
            RLogger::info("%s\n", message.toUtf8().constData());
            active = true;
            RProgressInitialize(message, pulse);
            lastUpdate = Clock::now();
        }

        // Open lazily so a phase with only a status line has no empty braces.
        void beginLog()
        {
            if (active && !logBlockOpen)
            {
                RLogger::indent();
                logBlockOpen = true;
            }
        }

        void finalize()
        {
            if (active)
            {
                RProgressFinalize();
                if (logBlockOpen)
                {
                    RLogger::unindent();
                    logBlockOpen = false;
                }
                active = false;
            }
        }

        void print(uint64_t step, uint64_t total)
        {
            const auto now = Clock::now();
            if (total > 0 && (step == total || now - lastUpdate >= std::chrono::milliseconds(250)))
            {
                RProgressPrint(step, total);
                lastUpdate = now;
            }
        }

        template <typename Message>
        void pulse(Message message)
        {
            if (Clock::now() - lastUpdate >= std::chrono::milliseconds(250))
            {
                // A work-count update belongs to the current phase, rather than
                // starting another timed log section.
                const QString text = message();
                beginLog();
                RLogger::info("%s\n", text.toUtf8().constData());
                RProgressInitialize(text, true);
                lastUpdate = Clock::now();
            }
        }
};

#endif
