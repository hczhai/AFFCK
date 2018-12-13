
#include "ck.h"

namespace affck {

	void RandomInit(int i) { 
		if (i < 0) srand((unsigned int)(time(NULL)));
		else srand((unsigned int)(i));
	}

	double RandomNumber() {
		return rand()*1.0 / RAND_MAX;
	}

	DisjointSet::DisjointSet(size_t n, size_t m) : rank(vector<int>(n)), parent(vector<int>(n)), m(m) {
		Clear();
	}

	void DisjointSet::Clear() {
		for (size_t i = 0; i < rank.size(); i++) {
			rank[i] = 0; parent[i] = i;
		}
	}

	int DisjointSet::FindX(int x) {
		if (parent[x] != x)
			parent[x] = FindX(parent[x]);
		return parent[x];
	}

	void DisjointSet::Union(int x, int y) {
		x = FindX(x), y = FindX(y);
		if (x == y) return;
		else if (x >= m && y < m)
			parent[y] = x;
		else if (y >= m && x < m)
			parent[x] = y;
		else if (rank[x] < rank[y])
			parent[x] = y;
		else if (rank[x]>rank[y])
			parent[y] = x;
		else {
			parent[y] = x;
			rank[x]++;
		}
	}

	const vector<vector<int> > CKPointGroup::Combination(int *times, int *maxs, int tn, int an) const {
		vector<int> ind(tn);
		vector<vector<int> > gens;
		for (int i = 0; i < tn; i++) ind[i] = 0;
		int l = 0;
		for (;;) {
			int j = 0;
			for (int k = 0; k < tn; k++) j += times[k] * ind[k];
			if (j == an) // an is the total number of atoms to generate 
				gens.push_back(ind);

			bool over = false; // generate next occupation
			for (int k = tn - 1; k >= 0; k--) {
				if (ind[k] < maxs[k]) {
					ind[k]++;
					for (int j = k + 1; j < tn; j++)
						ind[j] = 0;
					over = true;
					break;
				}
			}
			if (!over) break; // no more occupation
		}
		return gens;
	}

	int CKPointGroup::ArrayFind(const vector<int>& x, int xin, int v) const {
		for (int i = 0; i < xin; i++) {
			if (x[i] == v) return i;
		}
		return -1;
	}

	double& CKPointGroup::ATI(Atom& a, char dir) const {
		if (dir == 'x') return a.X;
		else if (dir == 'y') return a.Y;
		else if (dir == 'z') return a.Z;
		else
			throw ParseException("CK Point Group Error : Invalid direction:\n" + dir);
	}

	char CKPointGroup::Move(char dir)const {
		dir++;
		if (dir > 'z') dir = 'x';
		return dir;
	}

	char CKPointGroup::Extra(char dir1, char dir2) const { return 'x' + 'y' + 'z' - dir1 - dir2; }

	// reflection
	Atom CKPointGroup::CoI(const Atom& x) const { return Atom("", BoxSize - x.X, BoxSize - x.Y, BoxSize - x.Z); }

	// mirror respect to the plane dir = eq
	Atom CKPointGroup::CoMirror(const Atom& x, char dir, double eq) const {
		Atom cm = x;
		ATI(cm, dir) = eq * 2 - ATI(cm, dir);
		return cm;
	}

	// rotate along (bs/2, bs/2, *) (if dir = z), ang is in rad (counterclockwise)
	Atom CKPointGroup::CoRotate(const Atom& x, char dir, double ang)const {
		Atom ax("", BoxSize / 2, BoxSize / 2, BoxSize / 2), ay;
		ATI(ax, dir) = 0;
		Atom cor = x - ax;
		char idx = Move(dir), idy = Move(idx);
		ATI(ay, idx) = ATI(cor, idx) *cos(ang) - ATI(cor, idy) *sin(ang);
		ATI(ay, idy) = ATI(cor, idx) *sin(ang) + ATI(cor, idy) *cos(ang);
		ATI(ay, dir) = ATI(cor, dir);
		return ay - -ax;
	}

	// rotary reflection
	Atom CKPointGroup::CoMirrRot(const Atom& x, char dir, double ang)const {
		const Atom& tmp = CoRotate(x, dir, ang);
		return CoMirror(tmp, dir, BoxSize / 2);
	}

	CKPointGroup::CKPointGroup(double bs, const vector<pair<string, int> >& m) :BoxSize(bs), M(m) {}

	vector<Atom> CKPointGroup::PG(string sym, string dim) const {
		size_t m = M.size();
		vector<Atom> cd;
		Atom ac("", BoxSize / 2, BoxSize / 2, BoxSize / 2);

		for (size_t r = 0; r < m; r++) {
			int nct = M[r].second;
			vector<Atom> cdx(nct);
			for (int i = 0; i < nct; i++) {
				cdx[i].X = RandomNumber() * BoxSize;
				cdx[i].Y = RandomNumber() * BoxSize;
				cdx[i].Z = RandomNumber() * BoxSize;
			}
			if (sym == "C1");
			else if (sym.substr(0, 2) == "Cs") {
				if (sym.length()>3)
					throw ParseException("CK Point Group Error : Need Cs axis, Csx, Csy or Csz:\n" + sym);
				char idir = sym[2];
				if (nct % 2 == 1)
					ATI(cdx[nct / 2], idir) = BoxSize / 2;
				for (int i = 0; i < nct / 2; i++)
					cdx[nct - 1 - i] = CoMirror(cdx[i], idir, BoxSize / 2);
			} else if (sym.substr(0, 2) == "Ci") {
				if (nct % 2 == 1)
					cdx[nct / 2] = ac;
				for (int i = 0; i < nct / 2; i++)
					cdx[nct - 1 - i] = CoI(cdx[i]);
			} else if (sym.substr(0, 1) == "C") {
				size_t k = 1;
				double pi = acos(-1.0);
				for (; k < sym.length() && sym[k] >= '0'&&sym[k] <= '9'; k++);
				int g = atoi(sym.substr(1, k - 1).c_str());
				if (g <= 1 || k >= sym.length())
					throw ParseException("CK Point Group Error : Invalid n-fold axis:\n" + sym);
				else if (k == sym.length() - 1) { // Cn
					char idir = sym[k];
					for (int i = nct / g*g; i < nct; i++) {
						double dran = ATI(cdx[i], idir);
						cdx[i] = ac;
						ATI(cdx[i], idir) = dran;
					}
					for (int i = 1; i < g; i++) {
						for (int j = 0; j < nct / g; j++)
							cdx[i*(nct / g) + j] = CoRotate(cdx[j], idir, 2 * pi / g*i);
					}
				} else if (k + 1 == sym.length() - 1) { // Cnh
					char idir = sym[k];
					char idir2 = sym[k + 1];
					if (idir2 == 'h') {
						int dtimes[] = { 2 * g, g, 2, 1 };
						int dmaxs[] = { nct / (2 * g), nct / g, nct / 2, 1 };
						const vector<vector<int> >& dgens = Combination(dtimes, dmaxs, 4, nct);
						int dg = int(RandomNumber()*dgens.size());
						if (dg == dgens.size()) dg = 0;
						int j = 0;
						for (int i = 0; i < dgens[dg][3]; i++)
							cdx[j + i] = ac;
						j += dgens[dg][3];
						for (int i = 0; i < dgens[dg][2]; i++) {
							double dran = ATI(cdx[i + j], idir);
							cdx[i + j] = ac;
							ATI(cdx[i + j], idir) = dran;
							cdx[i + j + dgens[dg][2]] = CoMirror(cdx[i + j], idir, BoxSize / 2);
						}
						j += dgens[dg][2] * 2;
						for (int i = 0; i < dgens[dg][1]; i++) {
							ATI(cdx[i + j], idir) = BoxSize / 2;
							for (int l = 1; l < g; l++)
								cdx[i + j + dgens[dg][1] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
						}
						j += dgens[dg][1] * g;
						for (int i = 0; i < dgens[dg][0]; i++) {
							for (int l = 0; l < g; l++) {
								cdx[i + j + dgens[dg][0] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
								cdx[i + j + dgens[dg][0] * (l + g)] = CoMirror(cdx[i + j + dgens[dg][0] * l], idir, BoxSize / 2);
							}
						}
					} else
						throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
				} else if (k + 2 == sym.length() - 1) { // Cnv
					char idir = sym[k];
					char idir2 = sym[k + 1];
					char idir3 = sym[k + 2];
					if (idir2 == 'v') {
						if (idir == idir3)
							throw ParseException("CK Point Group Error : For Cn*v** point group, * and ** cannot be the same:\n" + sym);
						int dtimes[] = { 2 * g, g, 1 };
						int dmaxs[] = { nct / (2 * g), nct / g, nct };
						const vector<vector<int> >& dgens = Combination(dtimes, dmaxs, 3, nct);
						int dg = int(RandomNumber()*dgens.size());
						if (dg == dgens.size()) dg = 0;
						int j = 0;
						for (int i = 0; i < dgens[dg][2]; i++) {
							double dran = ATI(cdx[i + j], idir);
							cdx[i + j] = ac;
							ATI(cdx[i + j], idir) = dran;
						}
						j += dgens[dg][2];
						char idirm = Extra(idir, idir3);
						for (int i = 0; i < dgens[dg][1]; i++) {
							ATI(cdx[i + j], idirm) = BoxSize / 2;
							for (int l = 1; l < g; l++)
								cdx[i + j + dgens[dg][1] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
						}
						j += dgens[dg][1] * g;
						for (int i = 0; i < dgens[dg][0]; i++) {
							cdx[i + j + dgens[dg][0] * g] = CoMirror(cdx[i + j], idirm, BoxSize / 2);
							for (int l = 1; l < g; l++) {
								cdx[i + j + dgens[dg][0] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
								cdx[i + j + dgens[dg][0] * (l + g)] = CoRotate(cdx[i + j + dgens[dg][0] * g], idir, 2 * pi / g*l);
							}
						}
					} else
						throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
				} else
					throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
			} else if (sym.substr(0, 1) == "S") {
				size_t k = 1;
				double pi = acos(-1.0);
				for (; k < sym.length() && sym[k] >= '0'&&sym[k] <= '9'; k++);
				int g = atoi(sym.substr(1, k - 1).c_str());
				if (g <= 1 || k >= sym.length())
					throw ParseException("CK Point Group Error : Invalid n-fold axis:\n" + sym);
				else if (k == sym.length() - 1) { // Sn
					char idir = sym[k];
					for (int i = nct / g*g; i < nct; i++) {
						double dran = ATI(cdx[i], idir);
						cdx[i] = ac;
						ATI(cdx[i], idir) = dran;
					}
					for (int i = 1; i < g; i++) {
						for (int j = 0; j < nct / g; j++)
							cdx[i*(nct / g) + j] = CoMirrRot(cdx[(i - 1)*(nct / g) + j], idir, 2 * pi / g);
					}
				} else
					throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
			} else if (sym.substr(0, 1) == "D") {
				size_t k = 1;
				double pi = acos(-1.0);
				for (; k < sym.length() && sym[k] >= '0'&&sym[k] <= '9'; k++);
				int g = atoi(sym.substr(1, k - 1).c_str());
				if (g <= 1 || k >= sym.length())
					throw ParseException("CK Point Group Error : Invalid n-fold axis:\n" + sym);
				else if (k + 1 == sym.length() - 1 || (k + 2 == sym.length() - 1 && (sym[k + 2] == 'h' || sym[k + 2] == 'd'))) { // Dn & Dnh & Dnd
					char idir = sym[k];
					char idir2 = sym[k + 1];
					int ll = 0;
					if (k + 2 == sym.length() - 1 && sym[k + 2] == 'h') ll = 1;
					else if (k + 2 == sym.length() - 1 && sym[k + 2] == 'd') ll = 2;
					if (idir == idir2)
						throw ParseException("CK Point Group Error : For Dn* ** & Dn* **h & Dn* "
						"**d point group, * and ** cannot be the same:\n" + sym);
					int dtimes[] = { 2 * g, g, 2, 1 };
					int dmaxs[] = { nct / (2 * g), nct / g, nct / 2, 1 };
					const vector<vector<int> >& dgens = Combination(dtimes, dmaxs, 4, nct);
					int dg = int(RandomNumber()*dgens.size());
					if (dg == dgens.size()) dg = 0;
					int j = 0;
					for (int i = 0; i < dgens[dg][3]; i++)
						cdx[j + i] = ac;
					j += dgens[dg][3];
					for (int i = 0; i < dgens[dg][2]; i++) {
						double dran = ATI(cdx[i + j], idir);
						cdx[i + j] = ac;
						ATI(cdx[i + j], idir) = dran;
						cdx[i + j + dgens[dg][2]] = CoRotate(cdx[i + j], idir2, pi);
					}
					j += dgens[dg][2] * 2;
					for (int i = 0; i < dgens[dg][1]; i++) {
						double dran = ATI(cdx[i + j], idir2);
						cdx[i + j] = ac;
						ATI(cdx[i + j], idir2) = dran;
						for (int l = 1; l < g; l++)
							cdx[i + j + dgens[dg][1] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
					}
					j += dgens[dg][1] * g;
					char idirm = Extra(idir, idir2);
					for (int i = 0; i < dgens[dg][0]; i++) {
						if (ll == 1)
							ATI(cdx[i + j], idirm) = BoxSize / 2;
						else if (ll == 2) {
							ATI(cdx[i + j], idirm) = BoxSize / 2;
							cdx[i + j] = CoRotate(cdx[i + j], idir, 2 * pi / g / 4);
						}
						cdx[i + j + dgens[dg][0] * g] = CoRotate(cdx[i + j], idir2, pi);
						for (int l = 1; l < g; l++) {
							cdx[i + j + dgens[dg][0] * l] = CoRotate(cdx[i + j], idir, 2 * pi / g*l);
							cdx[i + j + dgens[dg][0] * (l + g)] = CoRotate(cdx[i + j + dgens[dg][0] * g], idir, 2 * pi / g*l);
						}
					}
				} else
					throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
			} else
				throw ParseException("CK Point Group Error : Unimplemented point group:\n" + sym);
			for (size_t i = 0; i < cdx.size(); i++) {
				cdx[i].Name = M[r].first;
				cd.push_back(cdx[i]);
			}
		}

		if (dim.length() == 1) {
			int g = atoi(dim.c_str());
			if (g<1 || g>3)
				throw ParseException("CK Point Group Error : Invalid dimension:\n" + dim);
			for (size_t i = 0; i < cd.size(); i++) {
				if (g <= 2) {
					cd[i].Z = BoxSize / 2;
					if (g <= 1) cd[i].Y = BoxSize / 2;
				}
			}
		} else if (dim.length() == 2) {
			int g = atoi(dim.c_str());
			if (g != 1)
				throw ParseException("CK Point Group Error : Invalid dimension:\n" + dim);
			char idir = dim[1];
			for (size_t i = 0; i < cd.size(); i++) {
				double dran = ATI(cd[i], idir);
				cd[i] = ac;
				ATI(cd[i], idir) = dran;
			}
		} else if (dim.length() == 3) {
			int g = atoi(dim.c_str());
			if (g != 2)
				throw ParseException("CK Point Group Error : Invalid dimension:\n" + dim);
			char idir = dim[1];
			char idir2 = dim[2];
			char idirm = Extra(idir, idir2);
			for (size_t i = 0; i < cd.size(); i++)
				ATI(cd[i], idirm) = BoxSize / 2;
		} else
			throw ParseException("CK Point Group Error : Invalid dimension:\n" + dim);

		return cd;
	}

	Atom CK::GeometryCenter(const vector<Atom>& la, size_t len) const {
		Atom c;
		if (len == 0) len = la.size();
		for (size_t j = 0; j < len; j++)
			c = c - -la[j];
		return c*(1.0 / len);
	}

	Atom CK::MassCenter(const vector<Atom>& vt, int i, DisjointSet *djs) const { // -2: only cluster; -1: all
		Atom c;
		double n = 0.0;
		for (size_t j = 0; j < (i == -2 ? TotC : vt.size()); j++) {
			if (i < 0 || djs->FindX(j) == i) {
				c = c - -vt[j] * AtomWeights[vt[j].Name];
				n += AtomWeights[vt[j].Name];
			}
		}
		return c*(1 / n);
	}

	bool CK::FragCheck(DisjointSet *djs) const {
		for (int i = 0; i < TotC; i++) {
			if (djs->FindX(i) != djs->FindX(0))
				return false;
		}
		return true;
	}

	bool CK::MinDistanceCheck(const vector<Atom>& x, size_t m) const { // only check cluster atoms
		for (size_t i = 0; i < m; i++) {
			for (size_t j = i + 1; j < x.size(); j++) {
				if (x[i].DistanceTo(x[j]) <= BondLength(&x[i], &x[j])*MinFactor)
					return false;
			}
		}
		return true;
	}

	double CK::BondLength(const Atom *x, const Atom *y) const {
		return Cluster::BLCovalent[x->Name] + Cluster::BLCovalent[y->Name];
	}
	
	void CK::AnyRotate(vector<Atom>& vt, size_t m, const Atom& mc, double rr, double rtheta, 
		double rphi) const {
		double alpha = rr * sin(rtheta) * cos(rphi);
		double beta = rr * sin(rtheta) * sin(rphi);
		double phi = rr * cos(rtheta);
		for (size_t i = 0; i < m; i++) {
			Atom cor = vt[i] - mc, ay;
			ay.X = cos(beta)*cos(phi) * cor.X - sin(phi)*cos(beta) * cor.Y
				+ sin(beta) * cor.Z;
			ay.Y = (sin(alpha)*sin(beta)*cos(phi) + cos(alpha)*sin(phi)) * cor.X 
				+ (-sin(alpha)*sin(beta)*sin(phi) + cos(alpha)*cos(phi)) * cor.Y
				- sin(alpha)*cos(beta)* cor.Z;
			ay.Z = (-sin(beta)*cos(alpha)*cos(phi) + sin(alpha)*sin(phi)) * cor.X
				+ (sin(beta)*cos(alpha)*sin(phi) + sin(alpha)*cos(phi)) * cor.Y
				+ cos(alpha)*cos(beta) * cor.Z;
			ay.Name = vt[i].Name;
			vt[i] = ay - -mc;
		}
	}
	
	bool CK::Shifting(vector<Atom>& vt, DisjointSet *djs, const Atom mc, vector<vector<Atom> >& traj) const {

		djs->Clear();

		for (int i = 0; i < TotC; i++) {
			for (size_t j = i + 1; j < vt.size(); j++) {
				if (djs->FindX(i) != djs->FindX(j)) {
					if (vt[i].DistanceTo(vt[j]) <= BondLength(&vt[i], &vt[j]))
						djs->Union(i, j);
				}
			}
		}

		if (Surface != NULL) {
			int fragz = -1;
			for (int i = 0; i < TotC; i++) {
				if (djs->FindX(i) < TotC)
					fragz = djs->FindX(i);
			}
			if (fragz == -1) return true; // all connected to surface
		} else {
			bool ok = true;
			for (int i = 1; i < TotC; i++) {
				if (djs->FindX(i) != djs->FindX(0)) {
					ok = false; break;
				}
			}
			if (ok) return true;
		}

		bool all_passed = true;
		bool all_zero = true;
		Atom zer;
		for (int i = 0; i < TotC; i++) {
			if (djs->FindX(i) == i) {
				const Atom& mcx = MassCenter(vt, i, djs);
				Atom arr = mc - mcx;
				double arr_len = arr.Length();
				if (arr_len > EPSF && arr_len > ShiftLength)
					arr = arr * (1 / arr_len*ShiftLength);
				if (arr_len > EPSF) all_zero = false;
				int nk = 0;
				for (int j = 0; j < TotC; j++) {
					if (djs->FindX(j) == i)
						vt[j] = vt[j] - -arr, nk++;
				}

				all_passed = false;
			}
		}

		static const double pi = acos(-1.0);
		AnyRotate(vt, TotC, MassCenter(vt, -2), CK::RotateAngle / 180.0 * pi, RandomNumber()* pi, RandomNumber() * 2 * pi);
		
		vector<Atom> vax;
		if (CK::WriteTrajectory) {
			for (int j = 0; j < TotC; j++) {
				vax.push_back(Surface != NULL ? vt[j] - SShifted : vt[j]);
			}
		}
		traj.push_back(vax);

		return all_passed || all_zero;
	}

	CKPointGroup *CK::GetCKPG() const {
		return new CKPointGroup(BoxSize, M);
	}

	CKResults CK::Generate(vector<Atom>& vc) {

		Atom zg = Atom("", 0, 0, ZGap);
		vector<Atom> vt(vc.size());
		for (size_t i = 0; i < vc.size(); i++) vt[i] = Surface != NULL ? vc[i] - -zg : vc[i];
		if (Surface != NULL) {
			vt.resize(vc.size() + Surface->size());
			for (size_t i = 0; i < Surface->size(); i++)
				vt[i + vc.size()] = (*Surface)[i] - -SShifted;
		}

		DisjointSet djs(vt.size(), vc.size());

		if (!MinDistanceCheck(vt, vc.size())) return InitMinimum;

		bool passed = false;
		
		vector<vector<Atom> > traj;

		Atom mc = MassCenter(vt, -2);
		if (Surface != NULL)
			if (!XYCOM) mc.X = SCenter.X, mc.Y = SCenter.Y;
		double zcenter = mc.Z;
		while (!Shifting(vt, &djs, mc, traj)) {
			mc = MassCenter(vt, -2);
			if (mc.Z < 0) { passed = true; break; } // all cluster.Z < 0 means one atom go extremely far
			if (Surface != NULL) {
				if (!XYCOM) mc.X = SCenter.X, mc.Y = SCenter.Y;
				zcenter -= ZDown;
				mc.Z = zcenter > 0 ? zcenter : 0;
			}
		}

		if (!FragCheck(&djs)) 
			return Fragment;
		if (!MinDistanceCheck(vt, vc.size())) 
			return FinalMinimum;
		if (passed) return PassedSurface;
		
		for (size_t i = 0; i < vc.size(); i++)
			vc[i] = Surface != NULL ? vt[i] - SShifted : vt[i];

		if (WriteTrajectory) {
			string trajx = "";
			for (size_t i = 0; i < traj.size(); i++)
				trajx += Cluster(GenerateResults(traj[i]), 0.0).ToXYZ();
			WriteFileContent(TrajOutDir + "/traj.xyz", trajx);
		}

		return Success;
	}

	const vector<Atom> CK::GenerateResults(const vector<Atom>& vc) const {
		bool ws = Surface != NULL && OutSurface;
		size_t y = ws ? vc.size() + Surface->size() : vc.size();
		vector<Atom> vs(y);
		for (size_t i = 0; i < vc.size(); i++) vs[i] = vc[i];
		if (ws) {
			for (size_t i = 0; i < Surface->size(); i++)
				vs[i + vc.size()] = (*Surface)[i];
		} else if (Surface == NULL) {
			double xmin = vs[0].X, ymin = vs[0].Y, zmin = vs[0].Z;
			for (size_t i = 1; i < vs.size(); i++) {
				if (vs[i].X < xmin) xmin = vs[i].X;
				if (vs[i].Y < ymin) ymin = vs[i].Y;
				if (vs[i].Z < zmin) zmin = vs[i].Z;
			}
			Atom mm("", xmin, ymin, zmin);
			for (size_t i = 0; i < vs.size(); i++)
				vs[i] = vs[i] - mm;
		}
		return vs;
	}

	vector<double> CK::InnerDistance(const vector<Atom>& va) {
		vector<double> vd;
		for (size_t i = 0; i < va.size(); i++) {
			for (size_t j = 0; j < i; j++)
				vd.push_back(va[i].DistanceTo(va[j]));
		}
		sort(vd.begin(), vd.end());
		return vd;
	}

	bool CK::TestDifferent(const vector<double>& a, const vector<double>& b) const {
		bool diff = false;
		double y = 0, z = 0;
		for (size_t j = 0; j < a.size(); j++) {
			y += abs(a[j] - b[j]);
			z += b[j];
			if (abs(a[j] - b[j])>DMax)
				return true;
		}
		if (y / z > DRel)
			return true;
		else return false;
	}

	void CK::InitSurface(vector<Atom> *surf) {

		Surface = surf;
		if (surf == NULL) return;

		if (surf->size() == 0)
			throw ParseException("CK Error : No atoms found on surface.");

		// count the number of top layer of surface and its z coordinate

		vector<vector<int> > vvi;
		for (size_t i = 0; i < Surface->size(); i++) {
			bool handled = false;
			for (size_t j = 0; j < vvi.size(); j++)
				if (abs((*Surface)[i].Z - (*Surface)[vvi[j][0]].Z) < EPS) {
					vvi[j].push_back(i);
					handled = true; break;
				}
			if (handled) continue;
			else {
				size_t g = 0;
				while (g < vvi.size() && (*Surface)[i].Z < (*Surface)[vvi[g][0]].Z) g++;
				vector<int> vi;
				vi.push_back(i); vvi.insert(vvi.begin() + g, vi);
			}
		}

		for (size_t i = 0; i < SurfFix.size(); i++)
			if (vvi.size() > (size_t)SurfFix[i])
				for (size_t j = 0; j < vvi[SurfFix[i]].size(); j++) {
					(*Surface)[vvi[SurfFix[i]][j]].G = 1;
					(*Surface)[vvi[SurfFix[i]][j]].GX = (*Surface)[vvi[SurfFix[i]][j]].GY = (*Surface)[vvi[SurfFix[i]][j]].GZ = 0;
				}

		for (size_t i = 0; i < SurfUFix.size(); i++)
			if (vvi.size() >(size_t)SurfUFix[i])
				for (size_t j = 0; j < vvi[SurfUFix[i]].size(); j++) {
					(*Surface)[vvi[SurfUFix[i]][j]].G = 0;
					(*Surface)[vvi[SurfUFix[i]][j]].GX = (*Surface)[vvi[SurfUFix[i]][j]].GY = (*Surface)[vvi[SurfUFix[i]][j]].GZ = 0;
				}

		for (size_t i = 0; i < SurfUOut.size(); i++)
			if (vvi.size() >(size_t)SurfUOut[i])
				for (size_t j = 0; j < vvi[SurfUOut[i]].size(); j++) {
					(*Surface)[vvi[SurfUOut[i]][j]].G = -1;
					(*Surface)[vvi[SurfUOut[i]][j]].GX = (*Surface)[vvi[SurfUOut[i]][j]].GY = (*Surface)[vvi[SurfUOut[i]][j]].GZ = 0;
				}

		vector<Atom> va;
		for (size_t i = 0; i < Surface->size(); i++)
			if ((*Surface)[i].G != -1) va.push_back((*Surface)[i]);

		Surface->clear();
		for (size_t i = 0; i < va.size(); i++)
			Surface->push_back(va[i]);

		Atom sc = GeometryCenter(*surf);
		sc.Z = (*Surface)[vvi[0][0]].Z;
		Atom st = Atom("", BoxSize / 2, BoxSize / 2, 0);
		st.X += (RandomNumber() - 0.5) * XShift;
		st.Y += (RandomNumber() - 0.5) * YShift;
		SCenter = st;
		SShifted = st - sc;
	}

	void CK::InitName(const string& sys) {
		vector<pair<string, string> > vps;
		bool need_new = true;

		if (sys.size() == 0)
			throw ParseException("CK Error : Cluster name is empty!");

		for (size_t i = 0; i < sys.size(); i++) {
			if (sys[i] >= '0' && sys[i] <= '9') {
				if (vps.size() == 0)
					throw ParseException("CK Error : Cluster name starts with digit: \n" + sys);
				else
					vps[vps.size() - 1].second += sys[i];
				need_new = true;
			} else if (need_new) {
				vps.push_back(pair<string, string>(string("") + (char)sys[i], ""));
				need_new = false;
			} else {
				vps[vps.size() - 1].first += sys[i];
			}
		}
		M = vector<pair<string, int> >(vps.size());
		for (size_t i = 0; i < vps.size(); i++) {
			string lo = ToLower(vps[i].first);
			if (AtomWeights.find(lo) == AtomWeights.end()
				|| Cluster::BLCovalent.find(lo) == Cluster::BLCovalent.end())
				throw ParseException("CK Error : Bad element name in cluster: \n" + vps[i].first);
			int hi = 1;
			if (vps[i].second != "")
				hi = StringToInt(vps[i].second);
			if (hi < 1)
				throw ParseException("CK Error : Bad element number: \n" + vps[i].second);
			M[i] = pair<string, int>(lo, hi);
		}

		// calculate box size and total number of atoms for cluster

		double box_size = 0.0;
		int totc = 0;
		for (size_t i = 0; i < M.size(); i++) {
			box_size += BSizeFactor * M[i].second * Cluster::BLCovalent[M[i].first];
			totc += M[i].second;
		}

		BoxSize = box_size;
		TotC = totc;

		Surface = NULL;
	}

}