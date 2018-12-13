
#include "parse.h"

namespace affck {

	string CKRun::ck_name;
	string CKRun::ck_symmetry;
	string CKRun::ck_dimension;
	string CKRun::ck_output_directory;
	string CKRun::ck_surface_file;
	string CKRun::ck_output_xyz_files;
	int CKRun::ck_output_start_number;
	int CKRun::ck_number;
	int CKRun::ck_fails_limit;

	string FFRun::ff_output_file;
	string FFRun::ff_input_file;

	string FFRun::standard_list_file;
	string FFRun::standard_xyz_files;
	int FFRun::standard_list_file_id_column;
	int FFRun::standard_list_file_energy_column;
	int FFRun::standard_list_sort;
	int FFRun::standard_list_start_number;
	int FFRun::standard_list_accept_number;
	double FFRun::standard_list_accept_ratio;
	string FFRun::standard_fitted_output;
	int FFRun::standard_fitted_output_sort;

	string FFRun::test_list_file;
	string FFRun::test_xyz_files;
	int FFRun::test_list_file_id_column;
	int FFRun::test_list_file_energy_column;
	int FFRun::test_list_sort;
	int FFRun::test_list_start_number;
	int FFRun::test_list_accept_number;
	double FFRun::test_list_accept_ratio;
	string FFRun::test_fitted_output;
	int FFRun::test_fitted_output_sort;

	string FFRun::opt_list_file;
	string FFRun::opt_xyz_files;
	int FFRun::opt_list_file_id_column;
	int FFRun::opt_list_sort;
	int FFRun::opt_list_start_number;
	int FFRun::opt_list_accept_number;
	double FFRun::opt_list_accept_ratio;
	string FFRun::opt_output_directory;
	string FFRun::opt_output_xyz_files;
	string FFRun::opt_output_list_file;
	int FFRun::opt_output_list_sort;

	double CK::ShiftLength = 0.2; // shift in each step
	double CK::RotateAngle = 8;
	double CK::MinFactor = 0.67; // the distance below minfactor * cov-radius will be discarded
	double CK::BSizeFactor = 4; // box size factor
	bool CK::WriteTrajectory = false; // write trajectory to traj.xyz
	double CK::EPSF = 1e-6; // the error bound of arr_len
	double CK::DMax = 0.7, CK::DRel = 0.03;
	map<string, double> CK::AtomWeights;
	vector<int> CK::SurfUFix, CK::SurfFix, CK::SurfUOut;
	
	double CK::ZGap = 1.0, CK::ZDown = 0.04;
	double CK::XShift = 0.0, CK::YShift = 0.0;
	bool CK::XYCOM = true; // xy use center of mass
	double CK::EPS = 1e-5; // the error bound of z coordinates of top layer
	string CK::TrajOutDir;
	bool CK::OutSurface = true;

	int CG::cg_imax = 200;
	double CG::cg_epsilon = 1e-7;
	int CG::cg_jmax = 10;
	double CG::cg_epsilon_n = 1e-6;

	int GLJ::lj_power = 9;
	int GLJ::lj_distance = 1;

}