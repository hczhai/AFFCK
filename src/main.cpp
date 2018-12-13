
#include "affck.h"

int main(int argc, char *argv[]) {
	affck::InitParams();
  // set default
  char *pvar;
  pvar = getenv("AFFCKHOME");
  if (string(pvar).length() != 0)
    affck::default_settings = string(pvar) + string("/res/default.txt");
	// find input
  if (argc == 1)
		affck::Run(&cin);
	else {
		if (argc == 3)
			affck::default_settings = string(argv[2]);
		if (!ifstream(argv[1]))
			cout << "Run Error: Cannot open input file." << endl;
		else {
			ifstream ifs(argv[1]);
			affck::Run(&ifs);
		}
	}
	return 0;
}
