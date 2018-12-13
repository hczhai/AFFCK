
#include "base.h"

namespace affck {

	bool Global::file_override = false;
	bool Global::path_override = false;
	int Global::random_seed = 1;
	const string Global::def_uniti = "angstrom";
	string Global::def_unit = "angstrom";

	const string DoubleToString(const double d) {
		stringstream ss;
		ss.precision(10);
		ss << d;
		return ss.str();
	}

	string ToString(const double d) {
		stringstream ss;
		ss.precision(4);
		ss << d;
		return ss.str();
	}

	const string IntToString(const int i) {
		stringstream ss;
		ss << i;
		return ss.str();
	}

	int StringToInt(const string& x) { return atoi(x.c_str()); }

	string Join(const set<string>& a, const string& x) {
		string r = "";
		for (set<string>::const_iterator i = a.begin(); i != a.end(); i++)
			r += *i + x;
		if (r.size() != 0) r.resize(r.size() - x.length());
		return r;
	}

	string Join(const vector<int>& a, const string& x) {
		string r = "";
		for (vector<int>::const_iterator i = a.begin(); i != a.end(); i++)
			r += IntToString(*i) + x;
		if (r.size() != 0) r.resize(r.size() - x.length());
		return r;
	}

	string Join(const vector<string>& a, const string& x) {
		string r = "";
		for (vector<string>::const_iterator i = a.begin(); i != a.end(); i++)
			r += *i + x;
		if (r.size() != 0) r.resize(r.size() - x.length());
		return r;
	}

	const string ToLower(const string& s) {
		string g = "";
		for (size_t i = 0; i < s.length(); i++)
			g += (s[i] >= 'a' && s[i] <= 'z') ? s[i] : (char)tolower(s[i]);
		return g;
	}

	const string Capital(const string& s) {
		if (s.length() == 0 || (s[0] >= 'A' && s[0] <= 'Z')) return s;
		return string("") + (char)toupper(s[0]) + s.substr(1);
	}


	bool MakeDirectory(const string& path) {
		return _mkdir(path.c_str()) == 0;
	}

	bool FileOrPathExists(const string& filename) {
		return _access(filename.c_str(), 6) == 0;
	}

	const string GetFileDirectory(const string& filename) {
		int cur = filename.length() - 1;
		for (; cur >= 0; cur--) {
			if (filename[cur] == '\\' || filename[cur] == '/')
				break;
		}
		if (cur >= 0) {
			if (cur == 0) return string(filename, 0, 1);
			else return string(filename, 0, cur);
		} else {
			char *buffer = _getcwd(NULL, 0);
			if (buffer == NULL)
				throw ParseException("Directory Error: Cannot get current work directory: \n" + filename);
			const string x = string(buffer);
			free(buffer);
			return x;
		}
	}

	const string GetUseableName(const string& name, bool ovr, bool extend) {
		string nn = name;
		if (nn == "" || nn == "/" || ovr || !extend) {
			if (nn != "" && nn != "/") {
				if (!FileOrPathExists(GetFileDirectory(nn))) {
					bool b = MakeDirectory(GetFileDirectory(nn));
					if (!b) throw ParseException("Directory Error: Cannot make directory: \n" + GetFileDirectory(nn));
				}
			}
			return nn;
		}
		if (nn.length() > 1 && (nn[nn.length() - 1] == '\\' || nn[nn.length() - 1] == '/'))
			nn = string(nn, 0, nn.length() - 1);
		int cur = nn.length() - 1;
		for (; cur >= 0; cur--) {
			if (nn[cur] == '\\' || nn[cur] == '/')
				break;
			else if (nn[cur] == '.')
				break;
		}
		string pre = nn, end = "";
		if (cur >= 0 && nn[cur] == '.') {
			pre = string(nn, 0, cur);
			end = string(nn, cur, string::npos);
		}
		cur = pre.length() - 1;
		for (; cur >= 0; cur--) {
			if (pre[cur] == '_' || pre[cur] < '0' || pre[cur] > '9')
				break;
		}
		int num = 0;
		if (cur >= 0 && pre[cur] == '_' && cur != pre.length() - 1) {
			num = atoi(string(pre, cur + 1, pre.length() - cur - 1).c_str());
			pre = string(pre, 0, cur);
		}
		bool k = nn == pre + "_0" + end;
		while (FileOrPathExists(pre + (num == 0 && !k ? "" : "_" + IntToString(num)) + end))
			num++;
		string ok = pre + (num == 0 && !k ? "" : "_" + IntToString(num)) + end;
		if (!FileOrPathExists(GetFileDirectory(ok))) {
			bool b = MakeDirectory(GetFileDirectory(ok));
			if (!b) throw ParseException("Directory Error: Cannot make directory: \n" + GetFileDirectory(ok));
		}
		return ok;
	}

	const string GetUseableFileName(const string& name, bool ovr, bool extend) {
		string u = GetUseableName(name, ovr, extend);
		string dir = GetFileDirectory(u);
		if (u.length() > dir.length() && u.substr(0, dir.length()) == dir)
			u = u.substr(dir.length());
		if (u.length() != 0 && (u[0] == '/' || u[0] == '\\'))
			u = u.substr(1);
		return u;
	}

	void WriteFileContent(const string& filename, const string& cont, bool extend) {
		const string cs = GetUseableName(filename, Global::file_override, extend);
		ofstream ofs(cs.c_str());
		ofs.write(cont.c_str(), cont.length());
		ofs.close();
	}

}