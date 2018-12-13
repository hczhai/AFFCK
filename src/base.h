#pragma once
#include <string>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#if defined(__WINDOWS__) || defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#define _access access
#define _mkdir(x) mkdir(x, S_IRWXU)
#define _getcwd getcwd
#endif

using namespace std;
namespace affck {

	class Global {
	public:
		static bool file_override;
		static bool path_override;
		static string def_unit;
		static const string def_uniti;
		static int random_seed;
	};


	const string DoubleToString(const double d);
	string ToString(const double d);
	const string IntToString(const int i);
	int StringToInt(const string& x);
	string Join(const set<string>& a, const string& x);
	string Join(const vector<int>& a, const string& x);
	string Join(const vector<string>& a, const string& x);
	const string ToLower(const string& s);
	const string Capital(const string& s);

	class ParseException { // exception
	public:
		const string x;
		ParseException(const string& x) : x(x) {}
	};

	bool MakeDirectory(const string& path);
	bool FileOrPathExists(const string& filename);
	const string GetFileDirectory(const string& filename);
	const string GetUseableName(const string& name, bool ovr, bool extend = true);
	const string GetUseableFileName(const string& name, bool ovr, bool extend = true);
	void WriteFileContent(const string& filename, const string& cont, bool extend = true);
}