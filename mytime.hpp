#ifndef _MYTIME_H_
	#define _MYTIME_H_
	
	#include <iostream>
	#include <sys/time.h>

	using namespace std;
 
	class Timerx
	{
		private:
			timeval t1, t2;

		public:
			void start()
			{
				gettimeofday(&t1, NULL);
			}

			void stop(const char* msg)
			{
				gettimeofday(&t2, NULL);

				double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    				elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
   
				cout << msg << "; Time = [" << elapsedTime << "] ms" << endl;
			}

	};
#endif
