#pragma once
#include <ctime>
#include "cluster.h"

namespace affck {

	void RandomInit(int i = -1);
	double RandomNumber();
	class DisjointSet { // disjoint-set for describe ck fragments
		vector<int> rank;
		vector<int> parent;
		int m;
	public:
		DisjointSet(size_t n, size_t m);
		int FindX(int x);
		void Clear();
		void Union(int x, int y);
	};

	class CKPointGroup { // ck point group symmetry generation
		const double BoxSize;
		const vector<pair<string, int> >& M;
		const vector<vector<int> > Combination(int *times, int *maxs, int tn, int an) const;
		int ArrayFind(const vector<int>& x, int xin, int v) const;
		double& ATI(Atom& a, char dir) const;
		char Move(char dir) const;
		char Extra(char dir1, char dir2) const;
		Atom CoI(const Atom& x) const; // reflection
		Atom CoMirror(const Atom& x, char dir, double eq) const; // mirror respect to the plane dir = eq
		Atom CoRotate(const Atom& x, char dir, double ang) const; // rotate along (bs/2, bs/2, *) (if dir = z), ang is in rad (counterclockwise)
		Atom CoMirrRot(const Atom& x, char dir, double ang) const; // rotary reflection
	public:
		CKPointGroup(double bs, const vector<pair<string, int> >& m);
		vector<Atom> PG(string sym, string dim) const;
	};

	enum CKResults { Success = -1, InitMinimum = 0, Fragment = 1, GhostMinimum = 2, PassedSurface = 3, FinalMinimum = 4, Similar = 5 };

	class CK { // ck main algorithm
	public:
		static double ShiftLength; // shift in each step
		static double RotateAngle; // rotate angle in each step
		static double MinFactor; // the distance below minfactor * cov-radius will be discarded
		static double BSizeFactor; // box size factor
		static bool WriteTrajectory; // write trajectory to traj.xyz
		static double EPSF; // the error bound of arr_len
		static double DMax, DRel;
		static map<string, double> AtomWeights;

		static double ZGap, ZDown;
		static double XShift, YShift;
		static bool XYCOM; // xy use center of mass
		static double EPS; // the error bound of z coordinates of top layer
		static string TrajOutDir;
		static bool OutSurface;
		static vector<int> SurfUFix, SurfFix, SurfUOut;
	private:
		vector<pair<string, int> > M; // number of each type of atoms in cluster
		vector<Atom> *Surface;
		Atom SShifted, SCenter;
		double BoxSize;
		int TotC;
		Atom GeometryCenter(const vector<Atom>& la, size_t len = 0) const;
		Atom MassCenter(const vector<Atom>& vt, int i = -1, DisjointSet *djs = NULL) const;
		bool FragCheck(DisjointSet *djs) const;
		bool MinDistanceCheck(const vector<Atom>& x, size_t m) const;
		double BondLength(const Atom *x, const Atom *y) const;
		bool Shifting(vector<Atom>& vt, DisjointSet *djs, const Atom mc, vector<vector<Atom> >& traj) const;
		void AnyRotate(vector<Atom>& vt, size_t m, const Atom& mc, double rr, double rtheta, double rphi) const;
	public:
		CKResults Generate(vector<Atom>& vc);
		const vector<Atom> GenerateResults(const vector<Atom>& vc) const;
		static vector<double> InnerDistance(const vector<Atom>& va);
		bool TestDifferent(const vector<double>& a, const vector<double>& b) const;
		void InitSurface(vector<Atom> *surf);
		void InitName(const string& sys);
		CKPointGroup *GetCKPG() const;
	};

}