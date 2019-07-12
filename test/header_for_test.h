#include <time.h>
#if OS_WIN
#include <windows.h>
#endif

#if OS_UNIX
long long stopwatch(int flag) {
  enum clock_unit { nano = 0, micro, milli, sec } unit;

  const long long NANOS = 1000000000LL;
  static struct timespec startTS, endTS;
  static long long diff = 0;

  /*
          여기서 단위 조정
          nano, micro, milli, sec
  */
  unit = micro;

  // start
  if (flag == 0) {
    diff = 0;
    if (-1 == clock_gettime(CLOCK_MONOTONIC, &startTS))
      printf("Failed to call clock_gettime\n");
  }
  // end
  else if (flag == 1) {
    if (-1 == clock_gettime(CLOCK_MONOTONIC, &endTS))
      printf("Failed to call clock_gettime\n");
    diff = NANOS * (endTS.tv_sec - startTS.tv_sec) +
           (endTS.tv_nsec - startTS.tv_nsec);

    switch (unit) {
      case nano:
        printf("elapsed time : % lld nano sec\n", diff);
        break;
      case micro:
        printf("elapsed time : % lld micro sec\n", diff / 1000);
        break;
      case sec:
        printf("elapsed time : % lld sec\n", diff / 1000000000);
        break;
      default:
        printf("elapsed time : % lld milli sec\n", diff / 100000);
        break;
    }
  } else {
    printf("wrong flag | 0 : start, 1 : end\n");
  }
}
#else

LARGE_INTEGER getFILETIMEoffset() {
  SYSTEMTIME s;
  FILETIME f;
  LARGE_INTEGER t;
  s.wYear = 1970;
  s.wMonth = 1;
  s.wDay = 1;
  s.wHour = 0;
  s.wMinute = 0;
  s.wSecond = 0;
  s.wMilliseconds = 0;
  SystemTimeToFileTime(&s, &f);
  t.QuadPart = f.dwHighDateTime;
  t.QuadPart <<= 32;
  t.QuadPart |= f.dwLowDateTime;
  return (t);
}

int clock_gettime(struct timeval *tv) {
  LARGE_INTEGER t;
  FILETIME f;
  double microseconds;
  static LARGE_INTEGER offset;
  static double frequencyToMicroseconds;
  static int initialized = 0;
  static BOOL usePerformanceCounter = 0;
  if (!initialized) {
    LARGE_INTEGER performanceFrequency;
    initialized = 1;
    usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
    if (usePerformanceCounter) {
      QueryPerformanceCounter(&offset);
      frequencyToMicroseconds =
          (double)performanceFrequency.QuadPart / 1000000.;
    } else {
      offset = getFILETIMEoffset();
      frequencyToMicroseconds = 10.;
    }
  }
  if (usePerformanceCounter)
    QueryPerformanceCounter(&t);
  else {
    GetSystemTimeAsFileTime(&f);
    t.QuadPart = f.dwHighDateTime;
    t.QuadPart <<= 32;
    t.QuadPart |= f.dwLowDateTime;
  }
  t.QuadPart -= offset.QuadPart;
  microseconds = (double)t.QuadPart / frequencyToMicroseconds;
  t.QuadPart = microseconds;
  tv->tv_sec = t.QuadPart / 1000000;
  tv->tv_usec = t.QuadPart % 1000000;
  return (0);
}

void stopwatch(int flag) {
  static struct timeval startTV, endTV;
  static long long diff;
  const long long MICRO = 1000000LL;

  if (flag == 0)
    clock_gettime(&startTV);
  else {
    clock_gettime(&endTV);
    diff = MICRO * (endTV.tv_sec - startTV.tv_sec) +
           (endTV.tv_usec - startTV.tv_usec);
    printf("elapsed time : %lld micro seconds\n", diff);
  }
}
#endif

long long get_micro_sec() {
  struct utime;
#if OS_UNIX
  clock_gettime(CLOCK_MONOTONIC, &utime);
#elif OS_WIN
  clock_gettime(&utime);
#endif
  return utime.tv_usec;
}
