#pragma once

#include <stdio.h>
#include <stdarg.h>

#define DEBUG_NEVER			12
#define DEBUG_ERROR			10
#define DEBUG_WARNING		9
#define DEBUG_IMPORTANT		8
#define DEBUG_DETAIL		6
#define DEBUG_MORE_DETAIL	4
#define DEBUG_ALWAYS		1
#define DEBUG_DEFAULT		1

void debug_info( int level, const char * fmt, ... );
