#pragma once

#include <thread>
#include <vector>

namespace lpn
{
class ThreadGroup
{
   public:
    explicit ThreadGroup();

    void ComputeAllTasks();

    template <typename Task>
    void AddTask(Task task)
    {
        threads_.emplace_back(std::move(task));
    }

    static int GetThreadAmount();

   private:
    std::vector<std::thread> threads_;
};
};  // namespace lpn
