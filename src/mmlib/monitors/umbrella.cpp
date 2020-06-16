#include "global.h"

#include "umbrella.h"

namespace Monitors{

	void printUmbrellaProfile(UmbrellaMonitor &umbrella, double binwidth){

		Maths::Histogram1D hist = umbrella.getHistogram(binwidth);
		std::deque<double> probability;
		hist.getProb(probability);

		for(size_t x = 0; x < hist.data_size(); x++) {
			double Q = hist.x1() + (double) x * (double) hist.xgridsize();

			printf("%shst%10.5lf\t%8d\t%10.8lf\t%10.8lf\t%10.8lf\n",
				umbrella.name.c_str(),
				Q,
				int(hist.data(x)),
				probability[x],
				-0.592*log(probability[x]),
				-0.592*log(probability[x]) - umbrella.ff->calcEnergyAtQ(Q) );
		}
	}
}

