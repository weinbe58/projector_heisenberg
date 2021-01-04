
#include "projector_heisenberg/heis.h"
#include <vector>

int main(int argc, char const *argv[])
{
	const int M = 1000;
	const int L = 10;
	std::vector<double> J_list(L,1.0),h(L*L);
	std::vector<int> bond_list(2*J_list.size()), phi(L);


	for(int i=0;i<J_list.size();i++){
		J_list[i] = 1.0;
		bond_list[2*i] = std::min(i,(i+1)%L);
		bond_list[2*i+1] = std::max(i,(i+1)%L);;
		phi[i] = (i%2 ? -1 : 1);


	}
	// setting up amplitudes for initial vbs wavefunction, amplitude product state.
	// this is a square matrix in which the i,j correspond to the amplitude between sites i and j.
	// Note: that in order for this code to work properly the amplitudes for i,j on the 
	// same sublattice must be set to 0. 
	int lin_ind = 0;
	for(int i=0;i<J_list.size();i++){
		for(int j=0;j<L;j++){
			if(phi[i] != phi[j]){ // amplitude decays as a power-law of distance between sites w/ PBC
				const double d = std::abs(i-j);
				h[lin_ind++] = 1.0/std::pow(std::min(d,L-d),1.5);
			}
			else{
				h[lin_ind++] = 0.0;
			}
		}
	}


	proj_heis qmc(L,J_list.size(),M,&J_list[0],&bond_list[0],&h[0]);

	for(int i=0;i<100;i++){
		qmc.MCsweep();
	}


	const int nbin = 10;
	std::vector<double> bins;
	std::cout << std::setprecision(7) << std::scientific;
	for(int i=0;i<nbin;i++){
		std::cout << std::setw(5) << i;
		double E = 0;
		size_t count = 0;
		for(int i=0;i<10000;i++){
			// qmc.print_opstr();
			qmc.MCsweep();
			double e = 0; 
			// qmc.clusters is a vector that labels each site with a number. each number corresponds to a unique 
			// loop that has passed through the site the center of the projector string. This gives the loop graph
			// for the overlap between the left and right VBS states after evolving in imaginary time. You can use 
			// this data to calculate various quantities in the VBS basis. 

			// two point correlator:
			// <L|S_i \cdot S_j|R> = 0.75*phi[i]*phi[j] if cluster[i]==cluster[j] else 0
			// phi[i] : the staggared signs for each sublattice.

			// the code below calculates the expectation value of the energy:
			for(int b=0;b<J_list.size();b++){
				const int i = bond_list[2*b];
				const int j = bond_list[2*b+1];
				if(qmc.clusters[i]==qmc.clusters[j]) 
					e += 0.75*phi[i]*phi[j]*J_list[b];
			}

			// dynamic update of mean value, E, which is energy density
			count++; 
			double delta = e/L - E;
			E += delta/count;
		}
		bins.push_back(E);

		std::cout << std::setw(20) << E << std::endl;

	}

	double init = 0;
	double E = std::accumulate(bins.begin(), bins.end(), init)/bins.size();
	double var = std::accumulate(bins.begin(), bins.end(), init,[E](double a , double b){return a + (b-E)*(b-E);})/(bins.size()-1);
	std::cout << std::setw(20) << E << std::setw(20) << std::sqrt(var) << std::endl;


	return 0;
}