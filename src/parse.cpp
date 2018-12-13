
#include "parse.h"

namespace affck {

	vector<string> Execute; // exe list
	vector<ParamTerm> ParamSet; // type and initial values of parameters

	void InitParams() {

		ParamSet.push_back(ParamTerm(&Global::file_override, "bool", "file_override", "fileor", "false"));
		ParamSet.push_back(ParamTerm(&Global::path_override, "bool", "path_override", "pathor", "false"));
		ParamSet.push_back(ParamTerm(&Global::random_seed, "int", "random_seed", "rands", "1"));
		ParamSet.push_back(ParamTerm(&Global::def_unit, "string", "default_unit", "defu", "angstrom"));

		ParamSet.push_back(ParamTerm(&GLJ::lj_power, "int", "lj_power", "ljp", "9"));
		ParamSet.push_back(ParamTerm(&GLJ::lj_distance, "int", "lj_distance", "ljd", "1"));

		ParamSet.push_back(ParamTerm(&CKRun::ck_name, "string", "ck_name", "ckname", ""));
		ParamSet.push_back(ParamTerm(&CKRun::ck_symmetry, "string", "ck_symmetry", "cksym", "C1"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_number, "int", "ck_number", "cknum", "100"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_dimension, "string", "ck_dimension", "ckdim", "3"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_surface_file, "string", "ck_surface_file", "cksf", ""));
		ParamSet.push_back(ParamTerm(&CKRun::ck_output_directory, "string", "ck_output_directory", "ckod", "ck_structs"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_output_xyz_files, "string", "ck_output_xyz_files", "ckox", "pos_#.xyz"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_output_start_number, "int", "ck_output_start_number", "ckos", "0"));
		ParamSet.push_back(ParamTerm(&CKRun::ck_fails_limit, "int", "ck_fails_limit", "ckfl", "10000"));

		ParamSet.push_back(ParamTerm(&CK::ShiftLength, "length", "ck_shift_length", "ckshift", "0.2 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::RotateAngle, "double", "ck_rotate_angle", "ckra", "8"));
		ParamSet.push_back(ParamTerm(&CK::MinFactor, "double", "ck_min_distance_factor", "ckmf", "0.67"));
		ParamSet.push_back(ParamTerm(&CK::DMax, "length", "ck_dmax", "ckdmax", "0.7 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::DRel, "double", "ck_drel", "ckdrel", "0.03"));
		ParamSet.push_back(ParamTerm(&CK::BSizeFactor, "double", "ck_boxsize_factor", "ckbf", "4"));
		ParamSet.push_back(ParamTerm(&CK::WriteTrajectory, "bool", "ck_write_trajectory", "ckwt", "false"));
		ParamSet.push_back(ParamTerm(&CK::EPSF, "length", "ck_fragment_pos_eps", "ckfpe", "1e-6 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::EPS, "length", "ck_top_layer_eps", "cktle", "1e-5 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::ZGap, "length", "ck_zgap", "ckzgap", "1.0 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::ZDown, "length", "ck_zdown", "ckzdown", "0.04 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::XShift, "length", "ck_xshift", "ckxsh", "0.0 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::YShift, "length", "ck_yshift", "ckysh", "0.0 angstrom"));
		ParamSet.push_back(ParamTerm(&CK::XYCOM, "bool", "ck_xy_center_of_mass", "ckxycom", "true"));
		ParamSet.push_back(ParamTerm(&CK::TrajOutDir, "string", "ck_trajectory_dir", "cktjdir", "trajs"));
		ParamSet.push_back(ParamTerm(&CK::OutSurface, "bool", "ck_out_surface", "ckouts", "true"));
		ParamSet.push_back(ParamTerm(&CK::SurfFix, "ilist", "ck_surface_fix_layer", "cksfix", ""));
		ParamSet.push_back(ParamTerm(&CK::SurfUFix, "ilist", "ck_surface_unfix_layer", "cksufix", ""));
		ParamSet.push_back(ParamTerm(&CK::SurfUOut, "ilist", "ck_surface_ignore_layer", "cksig", ""));

		ParamSet.push_back(ParamTerm(&FFRun::ff_input_file, "string", "ff_input_file", "ffipf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::ff_output_file, "string", "ff_output_file", "ffopf", ""));

		ParamSet.push_back(ParamTerm(&Cluster::BLMinimum, "length", "minimum_bond_length", "minbl", "0.40 angstrom"));
		ParamSet.push_back(ParamTerm(&Cluster::BLTolerance, "length", "bond_length_tolerance", "torbl", "0.45 angstrom"));
		ParamSet.push_back(ParamTerm(&Cluster::OpCode, "ilist", "functions", "fun", ""));

		ParamSet.push_back(ParamTerm(&CG::cg_imax, "int", "cg_imax", "cgimax", "200"));
		ParamSet.push_back(ParamTerm(&CG::cg_jmax, "int", "cg_jmax", "cgjmax", "10"));
		ParamSet.push_back(ParamTerm(&CG::cg_epsilon, "double", "cg_epsilon", "cgeps", "1e-7"));
		ParamSet.push_back(ParamTerm(&CG::cg_epsilon_n, "double", "cg_epsilon_n", "cgepsn", "1e-6"));

		ParamSet.push_back(ParamTerm(&FFRun::standard_list_file, "string", "standard_list_file", "stdlf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::standard_xyz_files, "string", "standard_xyz_files", "stdxf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_file_id_column, "int", "standard_list_file_id_column", "stdli", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_file_energy_column, "int", "standard_list_file_energy_column", "stdle", "1"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_start_number, "int", "standard_list_start_number", "stdstn", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_accept_number, "int", "standard_list_accept_number", "stdacn", "99999"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_accept_ratio, "double", "standard_list_accept_ratio", "stdacr", "1.00"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_list_sort, "int", "standard_list_sort", "stdsort", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::standard_fitted_output, "string", "standard_fitted_output", "stdout", ""));
		ParamSet.push_back(ParamTerm(&FFRun::standard_fitted_output_sort, "int", "standard_fitted_output_sort", "stdosort", "0"));

		ParamSet.push_back(ParamTerm(&FFRun::test_list_file, "string", "test_list_file", "testlf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::test_xyz_files, "string", "test_xyz_files", "testxf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_file_id_column, "int", "test_list_file_id_column", "testli", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_file_energy_column, "int", "test_list_file_energy_column", "testle", "1"));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_start_number, "int", "test_list_start_number", "teststn", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_accept_number, "int", "test_list_accept_number", "testacn", "99999"));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_accept_ratio, "double", "test_list_accept_ratio", "testacr", "1.00"));
		ParamSet.push_back(ParamTerm(&FFRun::test_list_sort, "int", "test_list_sort", "testsort", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::test_fitted_output, "string", "test_fitted_output", "testout", ""));
		ParamSet.push_back(ParamTerm(&FFRun::test_fitted_output_sort, "int", "test_fitted_output_sort", "testosort", "0"));

		ParamSet.push_back(ParamTerm(&FFRun::opt_list_file, "string", "opt_list_file", "optlf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::opt_xyz_files, "string", "opt_xyz_files", "optxf", ""));
		ParamSet.push_back(ParamTerm(&FFRun::opt_list_file_id_column, "int", "opt_list_file_id_column", "optli", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::opt_list_start_number, "int", "opt_list_start_number", "optstn", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::opt_list_accept_number, "int", "opt_list_accept_number", "optacn", "99999"));
		ParamSet.push_back(ParamTerm(&FFRun::opt_list_accept_ratio, "double", "opt_list_accept_ratio", "optacr", "1.00"));
		ParamSet.push_back(ParamTerm(&FFRun::opt_list_sort, "int", "opt_list_sort", "optsort", "0"));
		ParamSet.push_back(ParamTerm(&FFRun::opt_output_directory, "string", "opt_output_directory", "optod", ""));
		ParamSet.push_back(ParamTerm(&FFRun::opt_output_xyz_files, "string", "opt_output_xyz_files", "optox", ""));
		ParamSet.push_back(ParamTerm(&FFRun::opt_output_list_file, "string", "opt_output_list_file", "optol", ""));
		ParamSet.push_back(ParamTerm(&FFRun::opt_output_list_sort, "int", "opt_output_list_sort", "optosort", "0"));

	}

	const string ReadStream(istream *input) {
		string h;
		string r = "";
		for (; !input->eof();) {
			getline(*input, h);
			if (h == "EOF") break;
			size_t ind = h.find("!");
			if (ind != string::npos)
				h = string(h, 0, ind);
			while ((ind = h.find("\r")) != string::npos)
			    h.replace(ind, 1, "");
			r += h + " ";
		}
		return r;
	}

	string& Trim(string& s) {
		if (s.empty()) return s;
		s.erase(0, s.find_first_not_of(" \t"));
		s.erase(s.find_last_not_of(" \t") + 1);
		return s;
	}


	vector<string> Split(const string& s, const string& delim) {
		vector<string> ret = vector<string>();
		size_t last = 0;
		size_t index = s.find_first_of(delim, last);
		while (index != string::npos) {
			ret.push_back(s.substr(last, index - last));
			last = index + 1;
			index = s.find_first_of(delim, last);
		}
		if (index - last > 0)
			ret.push_back(s.substr(last, index - last));
		return ret;
	}

	const vector<string> TrimRemoveEmpty(const vector<string>& vv) {
		vector<string> v = vv;
		for_each(v.begin(), v.end(), Trim);
		vector<string>::const_iterator pend = remove(v.begin(), v.end(), "");
		v.resize(pend - v.begin());
		return v;
	}

	const vector<pair<string, MapSV> > RawParse(const string& x) {
		vector<pair<string, MapSV> > r;
		const vector<string>& progs = TrimRemoveEmpty(Split(x, "{}"));
		for (size_t i = 0; i < progs.size(); i++) {
			vector<string> ptbody = TrimRemoveEmpty(Split(progs[i], ":"));
			string title = "", body = "";
			if (ptbody.size() == 2)
				title = ptbody[0], body = ptbody[1];
			else if (ptbody.size() == 1)
				title = "parameter", body = ptbody[0];
			else {
				title = ptbody[0];
				ptbody.erase(ptbody.begin());
				body = Join(ptbody, ":");
			}
			const vector<string>& pterms = TrimRemoveEmpty(Split(body, ";"));
			vector<pair<string, vector<string> > > mterms;
			for (size_t j = 0; j < pterms.size(); j++) {
				const vector<string>& pterm = TrimRemoveEmpty(Split(pterms[j], "="));
				vector<string> pvals;
				if (pterm.size() == 2)
					pvals = TrimRemoveEmpty(Split(pterm[1], "|"));
				else if (pterm.size() == 1);
				else {
					vector<string> ptermg = pterm;
					ptermg.erase(ptermg.begin());
					pvals = TrimRemoveEmpty(Split(Join(ptermg, "="), "|"));
				}
				mterms.push_back(pair<string, vector<string> >(ToLower(pterm[0]), pvals));
			}
			r.push_back(pair<string, MapSV>(ToLower(title), MapSV(mterms)));
		}
		return r;
	}

	string PrintInt(int x, size_t len) {
		string r = IntToString(x);
		while (r.length() < len) r = "0" + r;
		return r;
	}

	vector<int> ParseIntList(const string& val) {
		const vector<string>& vg = TrimRemoveEmpty(Split(val, " \t"));
		vector<int> vi;
		for (size_t i = 0; i < vg.size(); i++) {
			const vector<string>& vgg = TrimRemoveEmpty(Split(vg[i], "~"));
			if (vgg.size() == 1) vi.push_back(atoi(vg[i].c_str()));
			else if (vgg.size() == 2) {
				int a = atoi(vgg[0].c_str());
				int b = atoi(vgg[0].c_str());
				for (int j = a; j <= b; j++)
					vi.push_back(j);
			} else throw ParseException("Parse Error: Value error: \n" + vg[i]);
		}
		return vi;
	}

	bool ParseBool(const string& val) {
		if (val == "F" || val == "f" || val == "false" || val == "0")
			return false;
		else
			return true;
	}

	double ParseUnit(const string& val, bool has_unit) {
		const vector<string>& vg = TrimRemoveEmpty(Split(val, " \t"));
		if (vg.size() != 1 && vg.size() != 2)
			throw ParseException("Parse Error: Value error: \n" + val);
		double d = atof(vg[0].c_str());
		const string& h = vg.size() == 2 ? ToLower(vg[1]) : (has_unit ? Global::def_unit : "a");
		if (h == "a" || h == "ang" || h == "angstrom");
		else if (h == "bohr" || h == "au" || h == "a.u.")
			d = d * 0.5291772108;
		else if (h == "pm")
			d = d * 0.01;
		else if (h == "m")
			d = d * 1e10;
		return d;
	}

	class CombineSet {
		vector<int> vi;
		vector<int> gi;
		bool is_start;
		bool is_valid;
	public:
		CombineSet() :vi(), gi() {}
		void Add(int i) { vi.push_back(i); }
		void Reset() {
			gi.clear();
			for (size_t i = 0; i < vi.size(); i++)
				gi.push_back(0);
			is_start = true;
			is_valid = vi.size() != 0;
		}

		CombineSet& operator++() {
			is_start = false;
			for (int i = vi.size() - 1; i >= 0; i--) {
				if (gi[i] < vi[i] - 1) {
					gi[i]++; 
					for (size_t j = i + 1; j < vi.size(); j++)
						gi[j] = 0;
					return *this;
				}
			}
			is_valid = false;
			return *this;
		}

		bool IsValid() const { return is_valid; }

		int operator[](int i) const { return gi[i]; }
		int Length() const { return vi.size(); }
		bool IsStart() const { return is_start; }
	};

	const vector<vector<pair<string, MapSV> > > GroupExecuteSolve(const vector<pair<string, MapSV> >& vpm) {
		vector<vector<pair<string, MapSV> > > r;
		vector<int> v;
		for (size_t i = 0; i < vpm.size(); i++) {
			if (vpm[i].first == "execute" || vpm[i].first == "exe") {
				if (v.size() == 0 || v[v.size() - 1] != i - 1)
					v.push_back(i);
				else
					v[v.size() - 1] = i;
			}
		}
		int k = 0;
		for (size_t i = 0; i < v.size(); i++) {
			vector<pair<string, MapSV> > g;
			for (; k <= v[i]; k++) g.push_back(vpm[k]);
			r.push_back(g);
		}
		return r;
	}

	const vector<vector<pair<string, vector<pair<string, string> > > > > MultiSolve(const vector<pair<string, MapSV> >& vpm) {
		CombineSet cs;
		for (size_t i = 0; i < vpm.size(); i++) {
			const vector<pair<string, vector<string> > >& msv = vpm[i].second.mapsv;
			for (size_t j = 0; j < msv.size(); j++) {
				const vector<string>& vs = msv[j].second;
				if (vs.size() > 1) cs.Add(vs.size());
			}
		}
		vector<vector<pair<string, vector<pair<string, string> > > > > r;
		if (cs.Length() == 0) cs.Add(1);
		cs.Reset();
		while (cs.IsValid()) {
			vector<pair<string, vector<pair<string, string> > > > rr;
			int k = 0;
			for (size_t i = 0; i < vpm.size(); i++) {
				const vector<pair<string, vector<string> > >& msv = vpm[i].second.mapsv;
				vector<pair<string, string> > vps;
				for (size_t j = 0; j < msv.size(); j++) {
					const vector<string>& vs = msv[j].second;
					if (vs.size() == 0 && (cs.IsStart() || vpm[i].first == "execute" || vpm[i].first == "exe"))
						vps.push_back(pair<string, string>(msv[j].first, ""));
					else if (vs.size() == 1 && (cs.IsStart() || vpm[i].first == "execute" || vpm[i].first == "exe"))
						vps.push_back(pair<string, string>(msv[j].first, vs[0]));
					else if (vs.size() > 1)
						vps.push_back(pair<string, string>(msv[j].first, vs[cs[k++]]));
				}
				rr.push_back(pair<string, vector<pair<string, string> > >(vpm[i].first, vps));
			}
			r.push_back(rr);
			++cs;
		}
		return r;
	}

	const string MultiPrint(const vector<pair<string, vector<pair<string, string> > > >& vps) {
		string r = "";
		for (size_t i = 0; i < vps.size(); i++) {
			string t = vps[i].first + ":\n";
			const vector<pair<string, string> >& ps = vps[i].second;
			for (size_t i = 0; i < ps.size(); i++)
				t += "   " + ps[i].first + (ps[i].second != "" ? " = " + ps[i].second : "") + ";\n";
			if (ps.size() != 0) r += "{ " + t + "}\n\n";
		}
		return r;
	}

	const string WriteEnergyList(const vector<pair<int, pair<double, double> > >& vpip) {
		string r = "";
		for (size_t i = 0; i < vpip.size(); i++)
			r += IntToString(vpip[i].first) + " " + DoubleToString(vpip[i].second.first)
			+ " " + DoubleToString(vpip[i].second.second) + "\n";
		return r;
	}

	const string WriteEnergyList(const vector<pair<int, double> >& vpi) {
		string r = "";
		for (size_t i = 0; i < vpi.size(); i++)
			r += IntToString(vpi[i].first) + " " + DoubleToString(vpi[i].second) + "\n";
		return r;
	}

	const string ReadFileContent(const string& filename) {
		const vector<string>& vg = TrimRemoveEmpty(Split(filename, " \t"));
		if (vg.size() == 1) {
			ifstream file(filename.c_str());
			if (!file)
				throw ParseException("File Parse Error: Cannot open input file: \n" + filename);
			string h;
			string r = "";
			for (; !file.eof();) {
				getline(file, h);
				r += h + "\n";
			}
			file.close();
			return r;
		} else if (vg.size() > 1 && vg[0] == "seq") {
			if (vg.size() == 2) {
				int a = atoi(vg[1].c_str());
				return IntToString(a) + "\n";
			} else if (vg.size() == 3) {
				int a = atoi(vg[1].c_str());
				int b = atoi(vg[2].c_str());
				if (b < a) throw ParseException("File Parse Error: Bad seq parameters b < a: \n" + filename);
				string r = "";
				for (int i = a; i <= b; i++) r += IntToString(i) + "\n";
				return r;
			} else if (vg.size() == 4) {
				int a = atoi(vg[1].c_str());
				int b = atoi(vg[2].c_str());
				int c = atoi(vg[3].c_str());
				if (c == 0) throw ParseException("File Parse Error: Bad seq parameters c == 0: \n" + filename);
				if (c > 0 && b < a) throw ParseException("File Parse Error: Bad seq parameters c > 0, b < a: \n" + filename);
				if (c < 0 && b > a) throw ParseException("File Parse Error: Bad seq parameters c < 0, b > a: \n" + filename);
				string r = "";
				for (int i = a; (c > 0 && i <= b) || (c < 0 && i >= b); i+=c) r += IntToString(i) + "\n";
				return r;
			} else throw ParseException("File Parse Error: Too many seq parameters: \n" + filename);
		} else {
			string r = "";
			for (size_t i = 0; i < vg.size(); i++) r += vg[i] + "\n";
			return r;
		}
	}

	const string ExpandFileName(const string& a, int i) {
		string b = a;
		for (;;) {
			size_t k = b.find('#'), l = k;
			for (; l != string::npos && b[l] == '#'; l++);
			if (k == string::npos) break;
			b.replace(k, l - k, PrintInt(i, l - k));
		}
		return b;
	}

	const vector<string> ExpandFileName(const string& a, const vector<int>& vi) {
		vector<string> r;
		for (size_t i = 0; i < vi.size(); i++)
			r.push_back(ExpandFileName(a,vi[i]));
		return r;
	}

	bool CmpTPairF(const pair<int, pair<double, double> >& a, const pair<int, pair<double, double> >& b) {
		return a.second.first < b.second.first;
	}

	bool CmpTPairS(const pair<int, pair<double, double> >& a, const pair<int, pair<double, double> >& b) {
		return a.second.second < b.second.second;
	}

	bool CmpPairF(const pair<int, double>& a, const pair<int, double>& b) {
		return a.first < b.first;
	}

	bool CmpPairS(const pair<int, double>& a, const pair<int, double>& b) {
		return a.second < b.second;
	}

	void PairSort(pair<vector<int>, vector<double> >& list, bool cp(const pair<int, double>&, const pair<int, double>&)) {
		vector<pair<int, double> > vpi;
		for (size_t i = 0; i < list.first.size(); i++)
		if (list.second.size() == list.first.size())
			vpi.push_back(pair<int, double>(list.first[i], list.second[i]));
		else
			vpi.push_back(pair<int, double>(list.first[i], 0.0));
		sort(vpi.begin(), vpi.end(), cp);
		for (size_t i = 0; i < list.first.size(); i++) {
			list.first[i] = vpi[i].first;
			if (list.second.size() == list.first.size())
				list.second[i] = vpi[i].second;
		}
	}

	Cluster* ReadCluster(const string& cont, const double energy) {
		const vector<string>& cx = TrimRemoveEmpty(Split(cont, "\n\r"));
		vector<Atom> va;
		for (size_t j = 2; j < cx.size(); j++) {
			const vector<string>& al = TrimRemoveEmpty(Split(cx[j], " \t"));
			if (al.size() != 4 && al.size() != 7)
				throw ParseException("XYZ File Read Error: number of column error: \n" + cx[j]);
			Atom at(al[0], ParseUnit(al[1], false), 
				ParseUnit(al[2], false), ParseUnit(al[3], false));
			if (al.size() == 7) {
				at.G = true; 
				at.GX = ParseUnit(al[4], false), at.GY = ParseUnit(al[5], false), at.GZ = ParseUnit(al[6], false);
			}
			va.push_back(at);
		}
		Cluster *clu = new Cluster(va, energy);
		return clu;
	}

	const pair<vector<int>, vector<double> > ReadEnergyList(const string& cont, int idc, int enc) {
		vector<int> ri;
		vector<double> rd;
		const vector<string>& vg = TrimRemoveEmpty(Split(cont, "\n\r"));
		size_t mc = max(idc, enc);
		for (size_t i = 0; i < vg.size(); i++) {
			const vector<string>& vgi = TrimRemoveEmpty(Split(vg[i], " \t"));
			if (vgi.size() <= mc)
				throw ParseException("File Read Error: No enough column ( " + IntToString(idc) +
				(enc != -1 ? ", " + IntToString(enc) : "") + " ) in list file: \n" + vg[i]);
			ri.push_back(atoi(vgi[idc].c_str()));
			if (enc != -1) rd.push_back(atof(vgi[enc].c_str()));
		}
		return pair<vector<int>, vector<double> >(ri, rd);
	}

	void SetParametersSingle(const pair<string, vector<pair<string, string> > >& pv) {
		if (pv.first == "execute" || pv.first == "exe") {
			for (size_t i = 0; i < pv.second.size(); i++)
				Execute.push_back(pv.second[i].first);
		} else if (pv.first == "covalent" || pv.first == "cov") {
			for (size_t i = 0; i < pv.second.size(); i++)
				Cluster::BLCovalent[pv.second[i].first] = ParseUnit(pv.second[i].second);
		} else if (pv.first == "weight" || pv.first == "wei") {
			for (size_t i = 0; i < pv.second.size(); i++)
				CK::AtomWeights[pv.second[i].first] = ParseUnit(pv.second[i].second, false);
		} else if (pv.first == "parameter") {
			for (size_t i = 0; i < pv.second.size(); i++) {
				string title = pv.second[i].first;
				bool handled = false;
				for (size_t j = 0; j < ParamSet.size(); j++)
					if (title == ParamSet[j].namea || title == ParamSet[j].nameb) {
						if (ParamSet[j].type == "int")
							*(int*)ParamSet[j].ptr = StringToInt(pv.second[i].second);
						else if (ParamSet[j].type == "double")
							*(double*)ParamSet[j].ptr = ParseUnit(pv.second[i].second, false);
						else if (ParamSet[j].type == "length")
							*(double*)ParamSet[j].ptr = ParseUnit(pv.second[i].second);
						else if (ParamSet[j].type == "string")
							*(string*)ParamSet[j].ptr = pv.second[i].second;
						else if (ParamSet[j].type == "bool")
							*(bool*)ParamSet[j].ptr = ParseBool(pv.second[i].second);
						else if (ParamSet[j].type == "ilist")
							*(vector<int>*)ParamSet[j].ptr = ParseIntList(pv.second[i].second);
						else
							throw ParseException("Parse Error: Bad parameter type: \n" + ParamSet[j].type);
						handled = true;
						break;
					}
				if (!handled) throw ParseException("Parse Error: Bad parameter name: \n" + title);
			}
		} else
			throw ParseException("Parse Error: Bad command name: \n" + pv.first);
	}

	ParamFitting::ParamFitting(const string& input, vector<double>& f) : LC(), Indices() {
		const vector<pair<string, MapSV> >& vps = RawParse(input);
		const vector<vector<pair<string, vector<pair<string, string> > > > >& vvp = MultiSolve(vps);
		if (vvp.size() != 1)
			throw ParseException("FF I/O Error: Multi-value in FF Fitting file.");
		for (size_t i = 0; i < vvp[0].size(); i++) {
			if (vvp[0][i].first == "final") {
				for (size_t j = 0; j < vvp[0][i].second.size(); j++)
					f.push_back(ParseUnit(vvp[0][i].second[j].first, false));
			} else if (vvp[0][i].first == "pn" && vvp[0][i].second.size() == 1) {
				PN = StringToInt(vvp[0][i].second[0].first);
			} else if (vvp[0][i].first == "indices") {
				for (size_t j = 0; j < vvp[0][i].second.size(); j++) {
					int fi = StringToInt(vvp[0][i].second[j].first);
					vector<int> se = ParseIntList(vvp[0][i].second[j].second);
					Indices[fi] = se;
				}
			} else
				throw ParseException("FF I/O Error: Bad data: \n" + vvp[0][i].first);
		}
	}

	void ResetParametersSingle(pair<string, vector<pair<string, string> > >& pv) {
		if (pv.first == "execute" || pv.first == "exe") {
			for (size_t i = 0; i < pv.second.size(); i++)
				pv.second[i].first = Execute[i];
		} else if (pv.first == "covalent" || pv.first == "cov") {
			for (size_t i = 0; i < pv.second.size(); i++)
				pv.second[i].second = DoubleToString(Cluster::BLCovalent[pv.second[i].first]) + " " + Global::def_uniti;
		} else if (pv.first == "weight" || pv.first == "wei") {
			for (size_t i = 0; i < pv.second.size(); i++)
				pv.second[i].second = DoubleToString(CK::AtomWeights[pv.second[i].first]);
		} else if (pv.first == "parameter") {
			for (size_t i = 0; i < pv.second.size(); i++) {
				string title = pv.second[i].first;
				bool handled = false;
				for (size_t j = 0; j < ParamSet.size(); j++)
					if (title == ParamSet[j].namea || title == ParamSet[j].nameb) {
						if (ParamSet[j].type == "int")
							pv.second[i].second = IntToString(*(int*)ParamSet[j].ptr);
						else if (ParamSet[j].type == "double")
							pv.second[i].second = DoubleToString(*(double*)ParamSet[j].ptr);
						else if (ParamSet[j].type == "length")
							pv.second[i].second = DoubleToString(*(double*)ParamSet[j].ptr) + " " + Global::def_uniti;
						else if (ParamSet[j].type == "string")
							pv.second[i].second = *(string*)ParamSet[j].ptr;
						else if (ParamSet[j].type == "bool")
							pv.second[i].second = *(bool*)ParamSet[j].ptr ? "true" : "false";
						else if (ParamSet[j].type == "ilist")
							pv.second[i].second = Join(*(vector<int>*)ParamSet[j].ptr, " ");
						else
							throw ParseException("Parse Error: Bad parameter type: \n" + ParamSet[j].type);
						handled = true;
						break;
					}
				if (!handled) throw ParseException("Parse Error: Bad parameter name: \n" + title);
			}
		} else
			throw ParseException("Parse Error: Bad command name: \n" + pv.first);
	}

}
