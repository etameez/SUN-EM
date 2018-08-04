
#include "timer.h"

Timer::Timer() {};

void Timer::startTimer()
{
	 this->start_time = std::chrono::high_resolution_clock::now(); 
}

void Timer::endTimer()
{
	this->end_time = std::chrono::high_resolution_clock::now();

	if(this->timer_flag)
	{
		this->elapsed_time = 
        std::chrono::duration_cast<std::chrono::duration<double>>(this->end_time - this->start_time);
		this->timer_flag = false;
	}
	else
	{
		this->elapsed_time += 
        std::chrono::duration_cast<std::chrono::duration<double>>(this->end_time - this->start_time);
	}
}

void Timer::printTime()
{
    std::cout << this->elapsed_time.count() << "seconds" << std::endl;
}

double Timer::saveTime()
{
	return this->elapsed_time.count();
}