#pragma once
#include "base.h"
#include "funcs.h"
#include "derivs.h"

namespace affck {

	extern const int AMP; // const used to denote FF functions
	extern const int LJ_BASE;
	extern const int STR_BASE;
	extern const int BA_BASE;
	extern const int BASE_AMP;

	struct Point { // 2 integers pair
		int X, Y;
		Point(const int x, const int y) : X(x), Y(y) {}
		Point() : X(0), Y(0) {}
		bool operator==(const Point& b) const {
			return this->X == b.X && this->Y == b.Y;
		}
	};

	struct Point3D :Point { // 3 integers
		int Z;
		Point3D(const int x, const int y, const int z) : Point(x, y), Z(z) {}
		Point3D() : Point(0, 0), Z(0) {}
		bool operator==(const Point3D& b) const {
			return this->X == b.X && this->Y == b.Y && this->Z == b.Z;
		}
	};

	struct Point4D :Point3D { // 4 integers
		int W;
		Point4D(const int x, const int y, const int z, const int w) : Point3D(x, y, z), W(w) {}
		Point4D() : Point3D(0, 0, 0), W(0) {}
		bool operator==(const Point4D& b) const {
			return this->X == b.X && this->Y == b.Y && this->Z == b.Z && this->W == b.W;
		}
	};

	class SymMap { // Helping class for dividing FF func types by elem's types
		vector<string> vec;
		int InnerFind(const string& x);
	public:
		SymMap();
		const string Find(size_t i) const;
		int Find(const string& a, const string& b);
		int Find(const string& a, const string& b, const string& c);
		int Find(const string& a, const string& b, const string& c, const string& d);
	};

	struct Atom { // atom name and coordinates data
		string Name;
		double X, Y, Z, GX, GY, GZ; int G;
		Atom();
		Atom(const string& name, const double x, const double y, const double z);
		Atom(const string& name, double x, double y, double z, double gx, double gy, double gz);
		double Length() const;
		const Atom operator*(const double k) const;
		double operator^(const Atom& a) const; // Dot product
		const Atom operator-(const Atom& a) const;
		const Atom operator-() const;
		double DistanceTo(const Atom& b) const;
		double AngleBy(const Atom& a, const Atom& b) const; // a - this - b
		double AngleTors(const Atom& a, const Atom& b, const Atom& c) const; // a - this - b - c
	};

	class Cluster { // cluster algorithms
	public:
		static vector<int> OpCode;
		static map<string, double> BLCovalent;
		static double BLTolerance;
		static double BLMinimum;
		static SymMap smlj, smstr, smba;
	private:
		const static Derivs::DL DL;
		const static Derivs::DTheta DTheta;
		const vector<Atom> LA;
		vector<Point> Bonds;
		vector<Point> ExBonds;
		vector<Point3D> BondAngles;
		vector<Point4D> TorsAngles;
		vector<int> BondsNum;
		vector<GFunction*> Funcs;
	private:
		Cluster& operator=(const Cluster& cluster);
	public:
		const double Energy;
		Cluster(const vector<Atom>& la, const double energy);
		~Cluster();
		const vector<GFunction*> GetFuncs() const;
		const string ToXYZ() const;
		double XF() const;
		const vector<double> XFp() const;
		const vector<vector<double> > XFpp() const;
		double FittedEnergy() const;
	private:
		template<class GFunc> void AFunc(const vector<int>& lix, const vector<int>& spi,
			const int num, int h, const double r, int base);
		vector<int> LClassify(int base);
		void GenFuncs();
	public:
		void GenBonds();
		const vector<double> RefX() const;
		static const Cluster* TransCluster(const vector<double>& x, const Cluster& refclu);
		const set<string> ElementSet() const;
		const vector<Atom>& GetLA() const;
	};

	class ParamFitting { // fitting algorithms
		const vector<Cluster*> LC;
		map<int, vector<int> > Indices; // from func's index to real parameters index list
		int PN;
	private:
		ParamFitting& operator=(const ParamFitting& pf);
	public:
		Cluster* RefCluster;
		vector<double>* RefSolution;
		ParamFitting(const vector<Cluster*>& lc, ParamFitting *pf);
		ParamFitting(const vector<Cluster*>& lc);
		ParamFitting(const string& input, vector<double>& f);
		~ParamFitting();
	private:
		const map<int, vector<double> > Trans(const vector<double>& x) const;
	public:
		const string PFOutput(const vector<double>& f) const;
		double F(const vector<double>& x) const;
		const vector<double> LinearB() const;
		const vector<vector<double> > LinearA(const vector<double>& x) const;
		void SetParam(const vector<double>& x) const;
		void SetParam(const Cluster& c, const vector<double>& x) const;
		const Cluster* TransCluster(const vector<double>& x) const;
		double CF(const vector<double>& x) const;
		const vector<double> CFp(const vector<double>& x) const;
		const vector<vector<double> > CFpp(const vector<double>& x) const;
		const vector<double> Init() const;
		const map<int, pair<int, GFunction*> > GFunctionStat() const;
	};

	void InitSymMap(const set<string>& ss);

}