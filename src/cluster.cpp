
#include "cluster.h"

namespace affck {

	const int AMP = 1000;
	const int LJ_BASE = 50;
	const int STR_BASE = 10;
	const int BA_BASE = 30;
	const int BASE_AMP = 15;

	int SymMap::InnerFind(const string& x) {
		for (size_t i = 0; i < vec.size(); i++) {
			if (vec[i] == x) return i;
		}
		vec.push_back(x);
		return vec.size() - 1;
	}

	SymMap::SymMap() { vec = vector<string>(); }

	int SymMap::Find(const string& a, const string& b) {
		string x = a + "-" + b;
		string y = b + "-" + a;
		return InnerFind(x < y ? x : y);
	}

	int SymMap::Find(const string& a, const string& b, const string& c) {
		string x = a + "-" + b + "-" + c;
		string y = c + "-" + b + "-" + a;
		return InnerFind(x < y ? x : y);
	}

	int SymMap::Find(const string& a, const string& b, const string& c, const string& d) {
		string x = a + "-" + b + "-" + c + "-" + d;
		string y = d + "-" + c + "-" + b + "-" + a;
		return InnerFind(x < y ? x : y);
	}

	const string SymMap::Find(size_t i) const {
		if (i < vec.size()) return vec[i];
		else return "";
	}

	Atom::Atom() : Name(""), X(0), Y(0), Z(0), G(0), GX(0), GY(0), GZ(0) {}

	Atom::Atom(const string& name, const double x, const double y, const double z)
		: Name(ToLower(name)), X(x), Y(y), Z(z), G(0), GX(0), GY(0), GZ(0) {
	}

	Atom::Atom(const string& name, double x, double y, double z, double gx, double gy, double gz)
		: Name(ToLower(name)), X(x), Y(y), Z(z), G(1), GX(gx), GY(gy), GZ(gz) {
	}

	double Atom::Length() const {
		return sqrt(X*X + Y*Y + Z*Z);
	}

	const Atom Atom::operator*(const double k) const {
		return Atom(this->Name, this->X*k, this->Y*k, this->Z*k);
	}

	// Dot product
	double Atom::operator^(const Atom& a) const {
		return a.X*this->X + a.Y*this->Y + a.Z*this->Z;
	}

	const Atom Atom::operator-(const Atom& a) const {
		return Atom(this->Name, this->X - a.X, this->Y - a.Y, this->Z - a.Z);
	}

	const Atom Atom::operator-() const {
		return Atom(this->Name, -this->X, -this->Y, -this->Z);
	}

	double Atom::DistanceTo(const Atom& b) const {
		double r = (this->X - b.X)*(this->X - b.X)
			+ (this->Y - b.Y)*(this->Y - b.Y)
			+ (this->Z - b.Z)*(this->Z - b.Z);
		if (r < 0) r = 0;
		return sqrt(r);
	}

	// a - this - b
	double Atom::AngleBy(const Atom& a, const Atom& b) const {
		const double lc = a.DistanceTo(b);
		const double la = this->DistanceTo(a);
		const double lb = this->DistanceTo(b);
		double g = (la*la + lb*lb - lc*lc) / (2 * la*lb);
		if (g>1)g = 1; else if (g < -1)g = -1;
		return acos(g);
	}

	// a - this - b - c
	double Atom::AngleTors(const Atom& a, const Atom& b, const Atom& c) const {
		const Atom tb = b - *this;
		const Atom bc = c - b;
		const Atom ta = a - *this;
		const Atom bh = tb * ((tb ^ bc) / tb.Length() / tb.Length());
		const Atom cc = bc - bh;
		const Atom th = (-tb) * (((-tb) ^ ta) / tb.Length() / tb.Length());
		const Atom aa = ta - th;
		double g = (aa^cc) / aa.Length() / cc.Length();
		if (g > 1)g = 1; else if (g < -1)g = -1;
		return acos(g);
	}

	Cluster& Cluster::operator = (const Cluster& cluster) { return *this; }
	Cluster::Cluster(const vector<Atom>& la, const double energy) : LA(la), Energy(energy) {}
	Cluster::~Cluster() {
		for (size_t i = 0; i < Funcs.size(); i++)
			delete Funcs[i];
	}

	const vector<GFunction*> Cluster::GetFuncs() const { return Funcs; }

	const string Cluster::ToXYZ() const {
		stringstream ss;
		ss.precision(10);
		ss.setf(ios::fixed);
		ss << LA.size() << endl << "ATOMIC_POSITIONS" << endl;
		for (size_t i = 0; i < LA.size(); i++) {
			ss.width(4);
			ss << left << Capital(LA[i].Name) << right;
			ss.width(18); ss << LA[i].X;
			ss.width(18); ss << LA[i].Y;
			ss.width(18); ss << LA[i].Z;
			if (LA[i].G) {
				ss.unsetf(ios::fixed);
				ss.width(10); ss << LA[i].GX;
				ss.width(10); ss << LA[i].GY;
				ss.width(10); ss << LA[i].GZ;
				ss.setf(ios::fixed);
			}
			ss << endl;
		}
		return ss.str();
	}

	double Cluster::XF() const {
		double r = 0;
		for (size_t i = 0; i < Funcs.size(); i++)
			r += Funcs[i]->Fl();
		return r;
	}

	const vector<double> Cluster::XFp() const {
		vector<double> r(LA.size() * 3);
		for (size_t i = 0; i < Funcs.size(); i++) {
			const GXL& fi = *Funcs[i];
			const double p = fi.Fpl();
			const vector<int>& fgi = fi.GetPointIndex();
			vector<double> x(fgi.size() * 3);
			for (size_t j = 0; j < fgi.size(); j++) {
				x[3 * j + 0] = LA[fgi[j]].X;
				x[3 * j + 1] = LA[fgi[j]].Y;
				x[3 * j + 2] = LA[fgi[j]].Z;
			}
			const Derivs::DFuntion *df = 0;
			if (fgi.size() == 2) df = &DL;
			else if (fgi.size() == 3) df = &DTheta;

			if (fgi.size() == 2 || fgi.size() == 3) {
				const vector<double>& lpx = df->Fp(x);
				for (size_t j = 0; j < fgi.size(); j++) {
					for (size_t k = 0; k < 3; k++)
						r[fgi[j] * 3 + k] += p * lpx[j * 3 + k];
				}
			}
		}
		return r;
	}

	const vector<vector<double> > Cluster::XFpp() const {
		vector<vector<double> > r(LA.size() * 3);
		for (size_t i = 0; i < LA.size() * 3; i++)
			r[i] = vector<double>(LA.size() * 3);
		for (size_t i = 0; i < Funcs.size(); i++) {
			const GXL& fi = *Funcs[i];
			const double p = fi.Fpl();
			const double pp = fi.Fppl();
			const vector<int>& fgi = fi.GetPointIndex();
			vector<double> x(fgi.size() * 3);
			for (size_t j = 0; j < fgi.size(); j++) {
				x[3 * j + 0] = LA[fgi[j]].X;
				x[3 * j + 1] = LA[fgi[j]].Y;
				x[3 * j + 2] = LA[fgi[j]].Z;
			}
			const Derivs::DFuntion *df = 0;
			if (fgi.size() == 2) df = &DL;
			else if (fgi.size() == 3) df = &DTheta;

			if (fgi.size() == 2 || fgi.size() == 3) {
				const vector<double>& lpx = df->Fp(x);
				const vector<vector<double> >& lppx = df->Fpp(x);
				for (size_t j = 0; j < fgi.size(); j++) {
					for (size_t k = 0; k < 3; k++) {
						for (size_t l = 0; l < fgi.size(); l++) {
							for (size_t m = 0; m < 3; m++) {
								r[fgi[j] * 3 + k][fgi[l] * 3 + m]
									+= p * lppx[j * 3 + k][l * 3 + m];
								r[fgi[j] * 3 + k][fgi[l] * 3 + m]
									+= pp * lpx[j * 3 + k] * lpx[l * 3 + m];
							}
						}
					}
				}
			}
		}
		return r;
	}

	double Cluster::FittedEnergy() const {
		double r = 0;
		for (size_t i = 0; i < Funcs.size(); i++)
			r += Funcs[i]->F(Funcs[i]->Param);
		return r;
	}

	const vector<Atom>& Cluster::GetLA() const { return LA; }

	template<class GFunc> void Cluster::AFunc(const vector<int>& lix, const vector<int>& spi,
		const int num, int h, const double r, int base) {
		if (lix.size() == 1)
			Funcs.push_back(new GFunc(r, h + lix[0], spi));
		else if (lix.size() >= 2) {
			bool k = false;
			if (num <= lix[0] - base) {
				Funcs.push_back(new GFunc(r, h + lix[0], spi));
				k = true;
			}
			for (size_t kkk = 1; kkk < lix.size(); kkk++) {
				if (num <= lix[kkk] - base && !k) {
					Funcs.push_back(new GFunc(r, h + lix[kkk], spi));
					k = true;
				}
			}
			if (!k) Funcs.push_back(new GFunc(r, h + lix[lix.size() - 1], spi));
		}
	}

	vector<int> Cluster::LClassify(int base) {
		vector<int> lix = vector<int>();
		for (vector<int>::const_iterator ll = OpCode.begin(); ll != OpCode.end(); ll++) {
			if (*ll >= base && *ll < base + BASE_AMP)
				lix.push_back(*ll);
		}
		sort(lix.begin(), lix.end());
		return lix;
	}

	void Cluster::GenFuncs() {

		Funcs = vector<GFunction*>();
		vector<int> ljd = LClassify(LJ_BASE);
		vector<int> lstr = LClassify(STR_BASE);
		vector<int> lba = LClassify(BA_BASE);

		for (size_t j = 0; j < ExBonds.size(); j++) {
			const Point& b = ExBonds[j];
			const double r = LA[b.X].DistanceTo(LA[b.Y]);
			vector<int> spi(2); spi[0] = b.X, spi[1] = b.Y;
			const int h = AMP* smlj.Find(LA[b.X].Name, LA[b.Y].Name);
			const vector<int>& lix = ljd;
			AFunc<GLJ>(lix, spi, BondsNum[b.X] + BondsNum[b.Y], h, r, LJ_BASE);
		}

		for (size_t j = 0; j < Bonds.size(); j++) {
			const Point& b = Bonds[j];
			const double r = LA[b.X].DistanceTo(LA[b.Y]);
			vector<int> spi(2); spi[0] = b.X, spi[1] = b.Y;
			const int h = AMP* smstr.Find(LA[b.X].Name, LA[b.Y].Name);
			const vector<int>& lix = lstr;
			AFunc<GP2C>(lix, spi, BondsNum[b.X] + BondsNum[b.Y], h, r, STR_BASE);
		}

		for (size_t j = 0; j < BondAngles.size(); j++) {
			const Point3D& a = BondAngles[j];
			const double r = LA[a.X].AngleBy(LA[a.Y], LA[a.Z]);
			vector<int> spi(3);
			spi[0] = a.Y, spi[1] = a.X, spi[2] = a.Z;
			const int h = AMP* smba.Find(LA[a.Y].Name, LA[a.X].Name, LA[a.Z].Name);
			const vector<int>& lix = lba;
			AFunc<GP2C>(lix, spi, BondsNum[a.X] + BondsNum[a.Y] + BondsNum[a.Z], h, r, BA_BASE);
		}

		Funcs.push_back(new GConst(-3));
	}

	void Cluster::GenBonds() {
		Bonds = vector<Point>();
		ExBonds = vector<Point>();
		BondAngles = vector<Point3D>();
		TorsAngles = vector<Point4D>();
		BondsNum = vector<int>();
		for (size_t i = 0; i < LA.size(); i++) BondsNum.push_back(0);
		for (size_t i = 0; i < LA.size(); i++) {
			for (size_t j = i + 1; j < LA.size(); j++) {
				const double dist = LA[i].DistanceTo(LA[j]);
				const double max = BLCovalent[LA[i].Name] + BLCovalent[LA[j].Name] + BLTolerance;
				if (dist <= max && dist >= BLMinimum) {
					Bonds.push_back(Point(i, j));
					BondsNum[i]++;
					BondsNum[j]++;
				}
			}
		}

		for (size_t i = 0; i < Bonds.size(); i++)
		for (size_t j = i + 1; j < Bonds.size(); j++) {
			vector<int> indices(4);
			indices[0] = Bonds[i].X, indices[1] = Bonds[i].Y, indices[2] = Bonds[j].X, indices[3] = Bonds[j].Y;
			sort(indices.begin(), indices.end());

			if (indices[0] == indices[1])
				BondAngles.push_back(Point3D(indices[0], indices[2], indices[3]));
			else if (indices[1] == indices[2])
				BondAngles.push_back(Point3D(indices[1], indices[0], indices[3]));
			else if (indices[2] == indices[3])
				BondAngles.push_back(Point3D(indices[2], indices[0], indices[1]));
			else {
				if (find(Bonds.begin(), Bonds.end(), Point(Bonds[i].X, Bonds[j].X)) != Bonds.end())
					TorsAngles.push_back(Point4D(Bonds[i].Y, Bonds[i].X, Bonds[j].X, Bonds[j].Y));
				if (find(Bonds.begin(), Bonds.end(), Point(Bonds[i].X, Bonds[j].Y)) != Bonds.end())
					TorsAngles.push_back(Point4D(Bonds[i].Y, Bonds[i].X, Bonds[j].Y, Bonds[j].X));
				if (find(Bonds.begin(), Bonds.end(), Point(Bonds[i].Y, Bonds[j].X)) != Bonds.end())
					TorsAngles.push_back(Point4D(Bonds[i].X, Bonds[i].Y, Bonds[j].X, Bonds[j].Y));
				if (find(Bonds.begin(), Bonds.end(), Point(Bonds[i].Y, Bonds[j].Y)) != Bonds.end())
					TorsAngles.push_back(Point4D(Bonds[i].X, Bonds[i].Y, Bonds[j].Y, Bonds[j].X));
			}
		}
		vector<vector<int> > dis(LA.size());
		for (size_t i = 0; i < LA.size(); i++) {
			dis[i] = vector<int>(LA.size());
			for (size_t j = 0; j < LA.size(); j++)
				dis[i][j] = i == j ? 0 : LA.size() * 2;
		}
		for (size_t i = 0; i < Bonds.size(); i++)
			dis[Bonds[i].X][Bonds[i].Y] = dis[Bonds[i].Y][Bonds[i].X] = 1;
		for (size_t k = 0; k < LA.size(); k++) {
			for (size_t i = 0; i < LA.size(); i++) {
				for (size_t j = 0; j < LA.size(); j++) {
					if (dis[i][k] + dis[k][j] < dis[i][j])
						dis[i][j] = dis[i][k] + dis[k][j];
				}
			}
		}
		for (size_t i = 0; i < LA.size(); i++) {
			for (size_t j = i + 1; j < LA.size(); j++) {
				if (dis[i][j] >= GLJ::lj_distance)
					ExBonds.push_back(Point(i, j));
			}
		}
		GenFuncs();
	}

	const vector<double> Cluster::RefX() const {
		vector<double> r(LA.size() * 3);
		for (size_t i = 0; i < LA.size(); i++) {
			r[i * 3 + 0] = LA[i].X;
			r[i * 3 + 1] = LA[i].Y;
			r[i * 3 + 2] = LA[i].Z;
		}
		return r;
	}

	const Cluster* Cluster::TransCluster(const vector<double>& x,
		const Cluster& refclu) {
		vector<Atom> la;
		for (size_t i = 0; i < refclu.LA.size(); i++)
			la.push_back(Atom(refclu.LA[i].Name, x[i * 3 + 0], x[i * 3 + 1], x[i * 3 + 2]));
		Cluster *clu = new Cluster(la, 0);
		clu->GenBonds();
		return clu;
	}

	const set<string> Cluster::ElementSet() const {
		set<string> ss;
		for (size_t i = 0; i < this->LA.size(); i++)
			ss.insert(this->LA[i].Name);
		return ss;
	}

	const Derivs::DL Cluster::DL = Derivs::DL();
	const Derivs::DTheta Cluster::DTheta = Derivs::DTheta();
	vector<int> Cluster::OpCode = vector<int>();
	map<string, double> Cluster::BLCovalent = map<string, double>();
	double Cluster::BLTolerance(0);
	double Cluster::BLMinimum(0);

	SymMap Cluster::smlj = SymMap(); // init directly according to type of elem's
	SymMap Cluster::smstr = SymMap();
	SymMap Cluster::smba = SymMap();

	void InitSymMap(const set<string>& ss) {
		Cluster::smlj = SymMap();
		Cluster::smba = SymMap();
		for (set<string>::const_iterator i = ss.begin(); i != ss.end(); i++) {
			for (set<string>::const_iterator j = ss.begin(); j != ss.end(); j++) {
				Cluster::smlj.Find(*i, *j);
				for (set<string>::const_iterator k = ss.begin(); k != ss.end(); k++)
					Cluster::smba.Find(*i, *j, *k);
			}
		}
		Cluster::smstr = Cluster::smlj;
	}

	ParamFitting& ParamFitting::operator=(const ParamFitting& pf) { return *this; }

	ParamFitting::ParamFitting(const vector<Cluster*>& lc) : LC(lc), Indices() {
		PN = 0;
		for (size_t i = 0; i < LC.size(); i++) {
			const vector<GFunction*>& vgf = LC[i]->GetFuncs();
			for (size_t j = 0; j < vgf.size(); j++) {
				const GFunction *gf = vgf[j];
				if (Indices.find(gf->Index) == Indices.end()) {
					vector<int> kk(gf->Num);
					for (int k = 0; k < gf->Num; k++)
						kk[k] = PN + k;
					PN += gf->Num;
					Indices[gf->Index] = kk;
				}
			}
		}
		if (int(LC.size()) < PN)
			throw ParseException("Fitting Error: Need at least " + IntToString(PN) + " structures for fitting.");
	}

	ParamFitting::ParamFitting(const vector<Cluster*>& lc, ParamFitting *pf) : LC(lc) {
		PN = pf->PN;
		Indices = pf->Indices;
	}

	ParamFitting::~ParamFitting() {
		for (size_t i = 0; i < LC.size(); i++)
			delete LC[i];
	}

	const string ParamFitting::PFOutput(const vector<double>& f) const {
		stringstream ss;
		ss.precision(15);
		ss << "! Automatic generated file for FF Fitting, please do not change." << endl << endl;
		ss << " { final: " << endl << "   ";
		for (size_t i = 0; i < f.size(); i++) ss << f[i] << "; ";
		ss << endl << " }" << endl;
		ss << " { pn: " << PN << " }" << endl;
		ss << " { indices: " << endl;
		for (map<int, vector<int> >::const_iterator i = Indices.begin(); i != Indices.end(); i++)
			ss << "   " << i->first << " = " << Join(i->second, " ") << ";" << endl;
		ss << " }" << endl;
		return ss.str();
	}

	const map<int, vector<double> > ParamFitting::Trans(const vector<double>& x) const {
		map<int, vector<double> > id;
		for (map<int, vector<int> >::const_iterator t = Indices.begin(); t != Indices.end(); ++t) {
			vector<double> g(t->second.size());
			for (size_t i = 0; i < g.size(); i++)
				g[i] = x[t->second[i]];
			id[t->first] = g;
		}
		return id;
	}

	double ParamFitting::F(const vector<double>& x) const {
		const map<int, vector<double> >& id = Trans(x);
		double r = 0;
		for (size_t i = 0; i < LC.size(); i++) {
			double rr = 0;
			const vector<GFunction*>& vgf = LC[i]->GetFuncs();
			for (size_t j = 0; j < vgf.size(); j++)
				rr += vgf[j]->F(id.at(vgf[j]->Index));
			r += (LC[i]->Energy - rr) * (LC[i]->Energy - rr);
		}
		return r;
	}

	const vector<double> ParamFitting::LinearB() const {
		vector<double> b(LC.size());
		for (size_t i = 0; i < LC.size(); i++)
			b[i] = LC[i]->Energy;
		return b;
	}

	const vector<vector<double> > ParamFitting::LinearA(const vector<double>& x) const {
		const map<int, vector<double> >& id = Trans(x);
		vector<vector<double> > a(LC.size());
		for (size_t i = 0; i < LC.size(); i++) {
			vector<double> g(PN);
			const vector<GFunction*>& vgf = LC[i]->GetFuncs();
			for (size_t j = 0; j < vgf.size(); j++) {
				const vector<double>& rr = vgf[j]->Fp(id.at(vgf[j]->Index));
				const vector<int>& t = Indices.at(vgf[j]->Index);
				for (size_t k = 0; k < t.size(); k++)
					g[t[k]] += rr[k];
			}
			a[i] = g;
		}
		return a;
	}

	void ParamFitting::SetParam(const vector<double>& x) const {
		const map<int, vector<double> >& id = Trans(x);
		for (size_t i = 0; i < LC.size(); i++) {
			const vector<GFunction*>& vgf = LC[i]->GetFuncs();
			for (size_t j = 0; j < vgf.size(); j++)
				vgf[j]->Param = id.at(vgf[j]->Index);
		}
	}

	void ParamFitting::SetParam(const Cluster& c, const vector<double>& x) const {
		const map<int, vector<double> >& id = Trans(x);
		const vector<GFunction*>& vgf = c.GetFuncs();
		for (size_t j = 0; j < vgf.size(); j++)
			vgf[j]->Param = id.at(vgf[j]->Index);
	}

	const Cluster* ParamFitting::TransCluster(const vector<double>& x) const {
		const Cluster *clu = Cluster::TransCluster(x, *RefCluster);
		SetParam(*clu, *RefSolution);
		return clu;
	}

	double ParamFitting::CF(const vector<double>& x) const { 
		const Cluster *c = TransCluster(x);
		const double d = c->XF();
		delete c;
		return d;
	}
	const vector<double> ParamFitting::CFp(const vector<double>& x) const { 
		const Cluster *c = TransCluster(x);
		const vector<double> d = c->XFp();
		delete c;
		return d;
	}
	const vector<vector<double> > ParamFitting::CFpp(const vector<double>& x) const { 
		const Cluster *c = TransCluster(x);
		const vector<vector<double> > d = c->XFpp();
		delete c;
		return d;
	}

	const vector<double> ParamFitting::Init() const {
		vector<double> r(PN);
		for (int i = 0; i < PN; i++)
			r[i] = 1;
		return r;
	}

	const map<int, pair<int, GFunction*> > ParamFitting::GFunctionStat() const {
		map<int, pair<int, GFunction*> > mpg;
		for (size_t i = 0; i < LC.size(); i++) {
			const vector<GFunction*>& vgf = LC[i]->GetFuncs();
			for (size_t j = 0; j < vgf.size(); j++)
			if (mpg.find(vgf[j]->Index) == mpg.end())
				mpg[vgf[j]->Index] = pair<int, GFunction*>(1, vgf[j]);
			else
				mpg[vgf[j]->Index].first++;
		}
		return mpg;
	}

}