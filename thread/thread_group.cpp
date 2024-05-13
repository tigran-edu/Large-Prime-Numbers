#include "thread_group.hpp"
#include <thread>

namespace lpn
{

ThreadGroup::ThreadGroup() : threads_() {}

void ThreadGroup::ComputeAllTasks()
{
    for (auto & thread : threads_)
    {
        thread.join();
    }
}

int ThreadGroup::GetThreadAmount() { return std::max<int>(2, std::thread::hardware_concurrency()) - 1; }

}  // namespace lpn
