
#include "affck.h"

namespace affck {

	ParamFitting *PF = NULL;
	vector<double> *Final = NULL;
	vector<vector<double> > CKSim;

	double ha_to_ev = 27.21138505;
	string default_settings = "default.txt";

	static const QF pcf = &ParamFitting::CF;
	static const QFp pcfp = &ParamFitting::CFp;
	static const QFpp pcfpp = &ParamFitting::CFpp;

	const string TranslateInteraction(int num) {
		int lo = num % AMP;
		int hi = num / AMP;
		string iname = "";
		if (lo >= LJ_BASE && lo < LJ_BASE + BASE_AMP)
			iname = "van der Waals ~" + IntToString(lo - LJ_BASE) + ": " + Cluster::smlj.Find(hi);
		else if (lo >= STR_BASE && lo < STR_BASE + BASE_AMP)
			iname = "stretching ~" + IntToString(lo - STR_BASE) + ": " + Cluster::smstr.Find(hi), lo -= STR_BASE;
		else if (lo >= BA_BASE && lo < BA_BASE + BASE_AMP)
			iname = "bending ~" + IntToString(lo - BA_BASE) + ": " + Cluster::smba.Find(hi), lo -= BA_BASE;
		else if (lo == -3)
			iname = "constant";
		else
			iname = "unknown";
		return iname;
	}

	pair<vector<Cluster*>, vector<int> > FFInner(const string& list_file, const string& xyz_files, int id_column, int energy_column,
		int list_sort, int accept_number, double accept_ratio, int start_number) {
		pair<vector<int>, vector<double> > el = ReadEnergyList(ReadFileContent(list_file),
			id_column, energy_column);
		int lines = el.first.size();
		cout << "   List lines: " << lines << endl;
		if (list_sort != 0) {
			PairSort(el, list_sort == 1 ? CmpPairF : CmpPairS);
			cout << "   List sorted." << endl;
		}
		int aclines = min(lines, int(accept_ratio*lines));
		aclines = min(aclines, accept_number);
		if (aclines - start_number < 1)
			throw ParseException("FF Error: Accepted energy list lines < 1: " + aclines - start_number);
		aclines -= start_number;
		if (start_number != 0) {
			el.first.erase(el.first.begin(), el.first.begin() + start_number);
			if (energy_column != -1) el.second.erase(el.second.begin(), el.second.begin() + start_number);
		}
		el.first.resize(aclines);
		if (energy_column != -1) el.second.resize(aclines);
		cout << "   Accepted energy list lines: " << aclines << endl;
		const vector<string>& fn = ExpandFileName(xyz_files, el.first);
		vector<Cluster*> vc = vector<Cluster*>();
		for (size_t i = 0; i < fn.size(); i++) {
			Cluster *clu = ReadCluster(ReadFileContent(fn[i]), energy_column != -1 ? el.second[i] : 0.0);
			vc.push_back(clu);
		}
		const set<string>& eles = vc[0]->ElementSet();
		cout << "   Cluster elements: " << Join(eles, " ") << endl;
		for (size_t i = 0; i < vc.size(); i++) {
			const set<string>& elev = vc[i]->ElementSet();
			if (eles.size() != elev.size() || !equal(eles.begin(), eles.end(), elev.begin()))
				throw ParseException("FF Error: Cluster elements mismatch #" +
				IntToString(el.first[i]) + ": " + Join(elev, " "));
		}
		cout << "   Cluster elements checking passed." << endl;
		InitSymMap(eles);
		for (size_t i = 0; i < vc.size(); i++)
			vc[i]->GenBonds();
		return pair<vector<Cluster*>, vector<int> >(vc, el.first);
	}

	void FFOpt() {
		string list_file = FFRun::opt_list_file;
		string xyz_files = FFRun::opt_xyz_files;
		int id_column = FFRun::opt_list_file_id_column;
		int energy_column = -1;
		int list_sort = FFRun::opt_list_sort;
		int accept_number = FFRun::opt_list_accept_number;
		double accept_ratio = FFRun::opt_list_accept_ratio;
		int start_number = FFRun::opt_list_start_number;

		cout << ">> AFFCK " << "Optimizing" << " Program <<" << endl << endl;

		const pair<vector<Cluster*>, vector<int> >& vcvi = FFInner(list_file, xyz_files, id_column,
			energy_column, list_sort, accept_number, accept_ratio, start_number);
		vector<Cluster*> vc = vcvi.first;
		vector<int> elf = vcvi.second;

		if (Final == NULL || PF == NULL)
			throw ParseException("Optimizing Error: Must run FF Fitting before optimizing.");

		string dir = GetUseableName(FFRun::opt_output_directory, Global::path_override);
		if (dir != "") {
			if (!FileOrPathExists(dir) && !MakeDirectory(dir))
				throw ParseException("Optimizing Error: Cannot make directory: \n" + dir);
		}

		const Cluster* vcf = NULL;
		vector<pair<int, pair<double, double> > > vpip(vc.size());
		for (size_t i = 0; i < vc.size(); i++) {
			cout << "   opt " << i + 1 << " of " << vc.size();
			PF->SetParam(*vc[i], *Final);
			PF->RefCluster = vc[i];
			PF->RefSolution = Final;
			int step = 0;
			vector<double> FinalGeo = CG::SolveNonLinear(pcf, pcfp, pcfpp, vc[i]->RefX(), &step, PF);
			vcf = PF->TransCluster(FinalGeo);
			if (dir != "" && FFRun::opt_output_xyz_files != "")
				WriteFileContent(dir + "/" + ExpandFileName(FFRun::opt_output_xyz_files, elf[i]), vcf->ToXYZ());
			cout << " (" << step << ") " <<
				vc[i]->FittedEnergy() << " " << vcf->XF() << endl;
			vpip[i] = pair<int, pair<double, double> >(elf[i], pair<double, double>(vc[i]->FittedEnergy(), vcf->XF()));
			delete vc[i];
			delete vcf;
		}

		if (FFRun::opt_output_list_file != "") {
			if (FFRun::opt_output_list_sort != 0) {
				sort(vpip.begin(), vpip.end(), FFRun::opt_output_list_sort == 1 ? CmpTPairF : CmpTPairS);
				cout << "   List sorted." << endl;
			}
			WriteFileContent(FFRun::opt_output_list_file, WriteEnergyList(vpip));
			cout << "   Optimized energiy list have been written. " << endl;
		}

		cout << endl << "-- End of AFFCK " << "Optimizing" << " Program --" << endl << endl;
	}

	string SolveVI(const vector<int>& vi) {
		if (accumulate(vi.begin(), vi.end(), 0) != 0)
			return " [" + IntToString(accumulate(vi.begin(), vi.end(), 0)) + ": "
			+ (vi[InitMinimum] != 0 ? "I=" + IntToString(vi[InitMinimum]) + " " : "")
			+ (vi[Fragment] != 0 ? "Fr=" + IntToString(vi[Fragment]) + " " : "")
			+ (vi[GhostMinimum] != 0 ? "G=" + IntToString(vi[GhostMinimum]) + " " : "")
			+ (vi[PassedSurface] != 0 ? "P=" + IntToString(vi[PassedSurface]) + " " : "")
			+ (vi[FinalMinimum] != 0 ? "Fi=" + IntToString(vi[FinalMinimum]) + " " : "")
			+ (vi[Similar] != 0 ? "S=" + IntToString(vi[Similar]) + " " : "") + "]";
		else return "";
	}

	void FFInput() {
		cout << ">> FF Input Program <<" << endl;
		if (FFRun::ff_input_file != "") {
			if (PF != NULL) delete PF;
			if (Final != NULL) delete Final;
			Final = new vector<double>();
			ifstream ifs(FFRun::ff_input_file.c_str());
			if (!ifs)
				throw ParseException("FF I/O Error: Cannot open file: " + FFRun::ff_input_file);
			const string& vps = ReadStream(&ifs);
			ifs.close();
			PF = new ParamFitting(vps, *Final);
		} else
			cout << "   Nothing did since parameter ff_input_file is empty." << endl;
		cout << "-- End of FF Input Program --" << endl << endl;
	}

	void FFOutput() {
		cout << ">> FF Output Program <<" << endl;
		if (FFRun::ff_output_file != "")
			WriteFileContent(FFRun::ff_output_file, PF->PFOutput(*Final));
		else
			cout << "   Nothing did since parameter ff_output_file is empty." << endl;
		cout << "-- End of FF Output Program --" << endl << endl;
	}

	void CKSimClear() {
		cout << ">> CKCS Similarity Clear Program <<" << endl;
		CKSim.clear();
		cout << "   Structure pool cleared." << endl;
		cout << "-- End of CKCS Similarity Clear Program --" << endl << endl;
	}

	void CKCS() {
		cout << ">> Coalescence Kick for Cluster with Surface (CKCS) Program <<" << endl << endl;
		if (Global::random_seed != -2) {
			if (Global::random_seed < -2)
				throw ParseException("CKCS Error: Invalid random seed: \n" + Global::random_seed);
			RandomInit(Global::random_seed);
			cout << "   Random seed set to: " << IntToString(Global::random_seed) << endl;
			Global::random_seed = -2;
		}
		CK ckx;
		ckx.InitName(CKRun::ck_name);
		cout << "   System name: " << CKRun::ck_name << endl;
		Cluster *clu = NULL; vector<Atom> clula;
		if (CKRun::ck_surface_file != "") {
			clu = ReadCluster(ReadFileContent(CKRun::ck_surface_file), 0.0);
			clula = clu->GetLA();
			ckx.InitSurface(&clula);
			cout << "   Surface file loaded." << endl;
		}

		string dir = GetUseableName(CKRun::ck_output_directory, Global::path_override);
		if (dir != "") {
			if (!FileOrPathExists(dir) && !MakeDirectory(dir))
				throw ParseException("CKCS Error: Cannot make directory: \n" + dir);
			cout << "   output directory: " << dir << endl << endl;
		}

		CKPointGroup *ckpg = ckx.GetCKPG();
		for (int i = 0; i < CKRun::ck_number; i++) {
			vector<int> vi(6);
			cout << "   generating " << i + 1 << " of " << CKRun::ck_number;
			bool suc = false;
			for (int j = 0; j < CKRun::ck_fails_limit; j++) {
				vector<Atom> va = ckpg->PG(CKRun::ck_symmetry, CKRun::ck_dimension);
				CKResults ckr = ckx.Generate(va);
				if (ckr != Success) { vi[ckr]++; continue; }
				const vector<double>& idis = CK::InnerDistance(va);
				bool similar = false;
				for (size_t k = 0; k < CKSim.size(); k++) {
					if (!ckx.TestDifferent(CKSim[k], idis)) {
						similar = true; break;
					}
				}
				if (similar) { vi[Similar]++; continue; }
				CKSim.push_back(idis);
				cout << SolveVI(vi);
				if (dir != "" && CKRun::ck_output_xyz_files != "") {
					Cluster clu(ckx.GenerateResults(va), 0.0);
					const string& q = ExpandFileName(CKRun::ck_output_xyz_files, CKRun::ck_output_start_number + i);
					cout << " -> " << GetUseableFileName(dir + "/" + q, Global::file_override, true) << endl;
					WriteFileContent(dir + "/" + q, clu.ToXYZ());
				} else
					cout << endl;
				suc = true;
				break;
			}
			if (!suc) {
				cout << SolveVI(vi);
				cout << " Failed. Meets the fails limit. Will not continue." << endl;
				cout << "!! WARNING: Please change the symmetry requirement, dmax, drel, or increase the fails limit!" << endl;
				break;
			}
		}
		cout << endl << "-- End of Coalescence Kick for Cluster with Surface (CKCS) Program --" << endl << endl;
		if (clu != NULL) delete clu;
		if (ckpg != NULL) delete ckpg;
	}

	void FFFit(bool test = false) {
		string list_file = test ? FFRun::test_list_file : FFRun::standard_list_file;
		string xyz_files = test ? FFRun::test_xyz_files : FFRun::standard_xyz_files;
		string fitted_output = test ? FFRun::test_fitted_output : FFRun::standard_fitted_output;
		int id_column = test ? FFRun::test_list_file_id_column : FFRun::standard_list_file_id_column;
		int energy_column = test ? FFRun::test_list_file_energy_column : FFRun::standard_list_file_energy_column;
		int list_sort = test ? FFRun::test_list_sort : FFRun::standard_list_sort;
		int accept_number = test ? FFRun::test_list_accept_number : FFRun::standard_list_accept_number;
		double accept_ratio = test ? FFRun::test_list_accept_ratio : FFRun::standard_list_accept_ratio;
		int start_number = test ? FFRun::test_list_start_number : FFRun::standard_list_start_number;
		int fitted_output_sort = test ? FFRun::test_fitted_output_sort : FFRun::standard_fitted_output_sort;

		cout << ">> AFFCK " << (test ? "Testing" : "Fitting") << " Program <<" << endl << endl;

		const pair<vector<Cluster*>, vector<int> >& vcvi = FFInner(list_file, xyz_files, id_column,
			energy_column, list_sort, accept_number, accept_ratio, start_number);
		vector<Cluster*> vc = vcvi.first;
		vector<int> elf = vcvi.second;

		if (PF != NULL) {
			ParamFitting *pfc = new ParamFitting(vc, PF);
			delete PF;
			PF = pfc;
		} else
			PF = new ParamFitting(vc);
		if (!test) {
			vector<double> *init = NULL;
			if (Final != NULL) {
				init = Final;
				cout << "   Start from last solution." << endl;
			} else {
				init = new vector<double>();
				*init = PF->Init();
			}
			Final = new vector<double>();
			*Final = CG::SolveLinearHH(PF->LinearA(*init), PF->LinearB());
			delete init;
		} else {
			vector<double> *init = new vector<double>();
			*init = PF->Init();
			if (Final == NULL)
				throw ParseException("Testing Error: Must run FF Fitting before Test.");
			if (Final->size() != init->size())
				throw ParseException("Testing Error: FF Fitting results mismatch with the data to be tested.");
			delete init;
		}
		PF->SetParam(*Final);
		double sqres = PF->F(*Final);
		double sig = sqrt(sqres / vc.size());
		cout << "   Sum of squared residual: " << sqres << endl;
		cout << "   Standard deviation: " << sig << endl;
		cout << "   Standard deviation (* htr->ev): " << sig * ha_to_ev << endl;
		if (fitted_output != "") {
			vector<pair<int, double> > vpi(vc.size());
			for (size_t i = 0; i < vc.size(); i++) {
				PF->SetParam(*vc[i], *Final);
				vpi[i] = pair<int, double>(elf[i], vc[i]->FittedEnergy());
			}
			if (fitted_output_sort != 0) {
				sort(vpi.begin(), vpi.end(), fitted_output_sort == 1 ? CmpPairF : CmpPairS);
				cout << "   List sorted." << endl;
			}
			WriteFileContent(fitted_output, WriteEnergyList(vpi));
			cout << "   Fitted energies have been written. " << endl;
		}
		if (!test) {
			cout << endl << "   Fitted functions:" << endl;
			const map<int, pair<int, GFunction*> >& mpg = PF->GFunctionStat();
			for (map<int, pair<int, GFunction*> >::const_iterator i = mpg.begin(); i != mpg.end(); i++) {
				cout << "   " << TranslateInteraction(i->first) << ":" << endl;
				cout << "     = $" << i->first << " (# " << i->second.first << " ) " <<
					i->second.second->Formula(i->second.second->Param) << endl << endl;
			}
		}
		cout << endl << "-- End of AFFCK " << (test ? "Testing" : "Fitting") << " Program --" << endl << endl;
	}

	void Run(istream *input) {
		time_t now_time;
		now_time = time(NULL);
		try {
			cout << "#################################################" << endl;
			cout << "#              AFFCK Program  2.1               #" << endl;
			cout << "#             Author: Huanchen ZHAI             #" << endl;
			cout << "#               stczhc@gmail.com                #" << endl;
			cout << "#                Nov. 21, 2015                  #" << endl;
			cout << "#################################################" << endl << endl;
			ifstream ifs(default_settings.c_str());
			if (!ifs)
				throw ParseException("Run Error: Cannot open default settings file: " + default_settings);
			const vector<pair<string, MapSV> >& vps = RawParse(ReadStream(&ifs));
			ifs.close();
			const vector<vector<pair<string, vector<pair<string, string> > > > >& vvp = MultiSolve(vps);
			if (vvp.size() != 1)
				throw ParseException("Run Error: Multi-value in default settings file");
			for (size_t i = 0; i < vvp[0].size(); i++)
				SetParametersSingle(vvp[0][i]);
			if (!*input)
				throw ParseException("Run Error: Cannot open input file.");
			const vector<pair<string, MapSV> >& vps2 = RawParse(ReadStream(input));
			const vector<vector<pair<string, MapSV> > >& vcs = GroupExecuteSolve(vps2);
			if (vcs.size() == 0)
				throw ParseException("Run Error: No execute instructions.");
			cout << "Total " << vcs.size() << " execution group(s) generated." << endl << endl << endl;
			for (size_t t = 0; t < vcs.size(); t++) {
				cout << " ** Exe group " << t + 1 << " of " << vcs.size() << endl << endl;
				vector<vector<pair<string, vector<pair<string, string> > > > > vvp2 = MultiSolve(vcs[t]);
				int n = vvp2.size();
				cout << "Total " << n << " set(s) of parameters generated." << endl << endl << endl;
				for (int k = 0; k < n; k++) {
					cout << " ** Run " << k + 1 << " of " << n << endl << endl;
					Execute.clear();
					for (size_t i = 0; i < vvp2[k].size(); i++) {
						SetParametersSingle(vvp2[k][i]);
						ResetParametersSingle(vvp2[k][i]);
					}
					cout << " Parameters set:" << endl << endl;
					cout << MultiPrint(vvp2[k]) << endl;
					for (size_t i = 0; i < Execute.size(); i++) {
						if (Execute[i] == "ff" || Execute[i] == "fitting" || Execute[i] == "fit" || Execute[i] == "fffit")
							FFFit();
						else if (Execute[i] == "test" || Execute[i] == "testing" || Execute[i] == "fftest")
							FFFit(true);
						else if (Execute[i] == "opt" || Execute[i] == "optimize" || Execute[i] == "relax"
							|| Execute[i] == "optimization" || Execute[i] == "ffopt" || Execute[i] == "optimizing")
							FFOpt();
						else if (Execute[i] == "ck" || Execute[i] == "ckcs")
							CKCS();
						else if (Execute[i] == "ckclear" || Execute[i] == "ckcsclear" ||Execute[i] == "clear")
							CKSimClear();
						else if (Execute[i] == "ffip" || Execute[i] == "ffinput")
							FFInput();
						else if (Execute[i] == "ffop" || Execute[i] == "ffoutput")
							FFOutput();
						else
							throw ParseException("Run Error: Cannot recognize program name: " + Execute[i]);
					}
				}
			}
			if (PF != NULL) delete PF;
			if (Final != NULL) delete Final;
		} catch (const ParseException& pe) {
			cout << "!! Error exit" << endl << pe.x << endl;
		} catch (const exception& e) {
			cout << "!! Fetal Error. Please contact the program author." << endl << e.what() << endl;
		}
		cout << endl << "Time used: " << time(NULL) - now_time << " sec" << endl;
		cout << endl << "## END ##" << endl;
	}

}
