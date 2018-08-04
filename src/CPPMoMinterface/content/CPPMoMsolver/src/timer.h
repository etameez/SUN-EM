#ifndef TIMER
#define TIMER

#include <chrono>
#include <iostream>

class Timer
{
	public:
		Timer();

		void startTimer();
		void endTimer();
		void printTime();
		double saveTime();

	protected:
		std::chrono::high_resolution_clock::time_point start_time;
		std::chrono::high_resolution_clock::time_point end_time;
		std::chrono::duration<double> elapsed_time;
		bool timer_flag = true;
};

#endif