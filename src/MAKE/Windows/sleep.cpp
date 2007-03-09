#include "sleep.h"
#include "windows.h"

void usleep (int x)
{
	int y = x;
	y = x/1000;

	Sleep(y);
}
