#include "Utility.h"
#include <stdarg.h>
#include <stdio.h>

void debug_info( int level, const char * fmt, ... )
{
	va_list arg_ptr;
	va_start( arg_ptr, fmt );

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL DEBUG_DEFAULT
#endif
	if(level > DEBUG_LEVEL)
	{
		vprintf( fmt, arg_ptr );
	}
	va_end( arg_ptr );
}
