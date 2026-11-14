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

public:
    explicit RTetraMeshProgress(const QString &message, bool pulse = true)
    {
        initialize(message,pulse);
    }

    ~RTetraMeshProgress() { RProgressFinalize(); }
    RTetraMeshProgress(const RTetraMeshProgress &) = delete;
    RTetraMeshProgress &operator=(const RTetraMeshProgress &) = delete;

    void initialize(const QString &message, bool pulse = true)
    {
        RLogger::info("%s\n",message.toUtf8().constData());
        RProgressInitialize(message,pulse);
        lastUpdate = Clock::now();
    }

    void print(uint64_t step, uint64_t total)
    {
        const auto now = Clock::now();
        if (total > 0 && (step == total || now-lastUpdate >= std::chrono::milliseconds(250)))
        {
            RProgressPrint(step,total);
            lastUpdate = now;
        }
    }

    template <typename Message>
    void pulse(Message message)
    {
        if (Clock::now()-lastUpdate >= std::chrono::milliseconds(250))
        {
            initialize(message());
        }
    }
};

#endif
