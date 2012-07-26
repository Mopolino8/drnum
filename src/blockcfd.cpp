#include "blockcfd.h"

unsigned long int global_flops     = 0;
unsigned long int global_flops_x86 = 0;

time_t global_start_time = time(NULL);

void startTiming()
{
  global_flops     = 0;
  global_flops_x86 = 0;
  cout << "\nTiming started.\n" << endl;
}

void stopTiming()
{
  int secs = time(NULL) - global_start_time;
  cout << "\nTiming finished:\n";
  if (global_flops > 0) {
    cout << "  floating point operations                : " << global_flops/1000000     << "M\n";
    cout << "  floating point operations (X86 weighted) : " << global_flops_x86/1000000 << "M\n";
    cout << "  seconds                                  : " << secs << "\n";
  } else {
    cout << "  seconds : " << secs << "\n";
  }
  cout << endl;
}


