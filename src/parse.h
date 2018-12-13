#pragma once
#include "cluster.h"
#include "ck.h"
#include "cg.h"

using namespace std;

namespace affck {

	class CKRun {
	public:
		static string ck_name;
		static string ck_symmetry;
		static string ck_dimension;
		static string ck_output_directory;
		static string ck_surface_file;
		static string ck_output_xyz_files;
		static int ck_output_start_number;
		static int ck_number;
		static int ck_fails_limit;
	};

	class FFRun {
	public:
		static string ff_output_file;
		static string ff_input_file;

		static string standard_list_file;
		static string standard_xyz_files;
		static int standard_list_file_id_column;
		static int standard_list_file_energy_column;
		static int standard_list_sort;
		static int standard_list_start_number;
		static int standard_list_accept_number;
		static double standard_list_accept_ratio;
		static string standard_fitted_output;
		static int standard_fitted_output_sort;

		static string test_list_file;
		static string test_xyz_files;
		static int test_list_file_id_column;
		static int test_list_file_energy_column;
		static int test_list_sort;
		static int test_list_start_number;
		static int test_list_accept_number;
		static double test_list_accept_ratio;
		static string test_fitted_output;
		static int test_fitted_output_sort;

		static string opt_list_file;
		static string opt_xyz_files;
		static int opt_list_file_id_column;
		static int opt_list_sort;
		static int opt_list_start_number;
		static int opt_list_accept_number;
		static double opt_list_accept_ratio;
		static string opt_output_directory;
		static string opt_output_xyz_files;
		static string opt_output_list_file;
		static int opt_output_list_sort;

	};

	class ParamTerm {
	public:
		void *ptr;
		string namea, nameb;
		string type, value;
		ParamTerm(void *p, string t, string na, string nb, string v) {
			ptr = p; type = t; namea = na; nameb = nb; value = v;
		}
	};

	extern vector<string> Execute; // exe list
	extern vector<ParamTerm> ParamSet; // type and initial values of parameters

	struct MapSV { // short type name for a pair vector
		vector<pair<string, vector<string> > > mapsv;
		MapSV() : mapsv(vector<pair<string, vector<string> > >()) {}
		MapSV(const vector<pair<string, vector<string> > >& mapsv) : mapsv(mapsv) {}
		MapSV& operator=(const MapSV& sv) {
			if (this != &sv) this->mapsv = sv.mapsv;
			return *this;
		}
	};

	// parsing methods

	const string IntToString(const int d);
	const string ToLower(const string& s);
	const string Capital(const string& s);
	double ParseUnit(const string& val, bool has_unit = true);
	int StringToInt(const string& x);
	vector<int> ParseIntList(const string& val);
	const string ReadStream(istream *input);
	const vector<pair<string, MapSV> > RawParse(const string& x);
	const vector<vector<pair<string, vector<pair<string, string> > > > > MultiSolve(const vector<pair<string, MapSV> >& vpm);
	const vector<vector<pair<string, MapSV> > > GroupExecuteSolve(const vector<pair<string, MapSV> >& vpm);
	const string MultiPrint(const vector<pair<string, vector<pair<string, string> > > >& vps);
	void SetParametersSingle(const pair<string, vector<pair<string, string> > >& pv);
	void ResetParametersSingle(pair<string, vector<pair<string, string> > >& pv);
	const pair<vector<int>, vector<double> > ReadEnergyList(const string& cont, int idc, int enc = -1);
	const string WriteEnergyList(const vector<pair<int, double> >& vpi);
	const string WriteEnergyList(const vector<pair<int, pair<double, double> > >& vpip);
	const string ReadFileContent(const string& filename);
	const vector<string> ExpandFileName(const string& a, const vector<int>& vi);
	const string ExpandFileName(const string& a, int i);
	Cluster* ReadCluster(const string& cont, const double energy);
	void PairSort(pair<vector<int>, vector<double> >& list, bool cp(const pair<int, double>&, const pair<int, double>&));
	string Join(const set<string>& a, const string& x);
	string Join(const vector<int>& a, const string& x);

	// comparing methods for sorting

	bool CmpPairF(const pair<int, double>& a, const pair<int, double>& b);
	bool CmpPairS(const pair<int, double>& a, const pair<int, double>& b);
	bool CmpTPairF(const pair<int, pair<double, double> >& a, const pair<int, pair<double, double> >& b);
	bool CmpTPairS(const pair<int, pair<double, double> >& a, const pair<int, pair<double, double> >& b);

	void InitParams();
}