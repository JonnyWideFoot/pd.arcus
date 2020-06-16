////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2005 by Andrei Alexandrescu
// Copyright (c) 2006 Peter Kümmel
// Permission to use, copy, modify, distribute, and sell this software for any
//     purpose is hereby granted without fee, provided that the above copyright
//     notice appear in all copies and that both that copyright notice and this
//     permission notice appear in supporting documentation.
// The author makes no representations about the suitability of this software 
//     for any purpose. It is provided "as is" without express or implied 
//     warranty.
// Modified by Jon Rea for incorperation into MMLib as the type-safe formatting 
// engine on 30/11/06
// 
////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include <cassert>
#include "stringbuilder.h"
#include "safeformat.h"

LokiCore::PrintfState<std::FILE*, char> Printf(const char* format) {
    return LokiCore::PrintfState<std::FILE*, char>(stdout, format);
}

LokiCore::PrintfState<std::FILE*, char> Printf(const std::string format) {
    return LokiCore::PrintfState<std::FILE*, char>(stdout, format.c_str());
}

LokiCore::PrintfState<std::FILE*, char> FPrintf(FILE* f, const char* format) {
    return LokiCore::PrintfState<std::FILE*, char>(f, format);
}

LokiCore::PrintfState<std::FILE*, char> FPrintf(FILE* f, const std::string& format) {
    return LokiCore::PrintfState<std::FILE*, char>(f, format.c_str());
}

LokiCore::PrintfState<std::string&, char> SPrintf(std::string& s, const char* format) {
    return LokiCore::PrintfState<std::string&, char>(s, format);
}

LokiCore::PrintfState<std::string&, char> SPrintf(std::string& s, const std::string& format) {
    return LokiCore::PrintfState<std::string&, char>(s, format.c_str());
}

namespace LokiCore
{	
	// Crude writing method: writes straight to the file, unbuffered
	// Must be combined with a buffer to work properly (and efficiently)

	void write(std::FILE* f, const char* from, const char* to) {
		assert(from <= to);
		fwrite(from, 1, to - from, f);
	}

	// Write to a string
	void write(std::string& s, const char* from, const char* to) {
		assert(from <= to);
		s.append(from, to);
	}

	// Write to a StringBuilder
	void write(StringBuilder& sb, const char* from, const char* to) {
		assert(from <= to);
		sb.append(from,to - from);
	}
}
