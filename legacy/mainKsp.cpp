#include "rocketGuidance.hpp"

int main()
{
    RocketGuidance guidance;
    guidance.run();
    while (true)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
}
