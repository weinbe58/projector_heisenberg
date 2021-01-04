#ifndef _heis_h_
#define _heis_h_

#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <cmath>
#include <vector>
#include "random.h"
#include "utils.h"


// implementation for proejctor QMC on Heisenberg model.
class proj_heis
{

  std::unique_ptr<select_random<double>> select_bond;

  std::vector<std::unique_ptr<select_random<double>>> select_h;

  rand_uniform ran;

  std::vector<int> X,opstr,end_caps,first,last;
  std::vector<signed char> spins,l_spins,r_spins;
  int *Vl,*Vr;
  const int * bond_list;
  const double * J_list, * h;
  const int N,Nb,M,H;

  	// select_bond : random number generator that generates integers based on weights.
    // ran : random number between 0 and 1;

	// N : number of sites on the lattice (must be even)
	// Nb : number of bonds connecting two sites on lattice
	// M : number of projector steps: H^M.

	// X : (2*M,) array for linked list for traversing loops
	// opstr : (2*M,) array which stores the type and bond operators at each time slice.
	// Vl/Vr : (N,) array, vbs configutation of left/right ket
	// l_spins/r_spins : (N,) array, spin configutation of left/right ket
	// spins : (N,) array, current spin-configutation during during imaginary time sweep
	// first,last : (N,) array,  list used to construct linked list
	// bond_list : (Nb,2) array which contains the sites associated with a particular bond
	// J_list : (Nb,) array which contains the coupling for each bond in bond_list
	// h : (N,N) array that contains the weight >=0 for the bond connecting sites i and j

	void init(); // initialize
	void diag_update(); // diagonal update for operator string
	void loop_update(); // loop update for operator string
	void linked_list(); // set up linked list

	template<class callback>
	void diag_update(callback&); // diagonal update for operator string
	void ends_update(); // update vbs configurations   

public:

	std::vector<int> clusters;

	proj_heis(const int N_,const int Nb_,const int M_,
		const double * J_list_,const int * bond_list_,const double * h_) : 
	N(N_), Nb(Nb_), M(M_), H(8*M_), J_list(J_list_), bond_list(bond_list_), h(h_) {this->init();}
	~proj_heis(){}



	template<class callback>
	void MCsweep(callback&);

	void MCsweep();
	void print_opstr();
	
};



void proj_heis::init() {

	assert(N%2==0);

	for(int b=0;b<2*Nb;b++){assert(bond_list[b] < N);}


	std::vector<double> J_abs(Nb);

	for(int i=0;i<Nb;i++){J_abs[i] = std::abs(J_list[i]);}

	select_bond = std::make_unique<select_random<double>>(&J_abs[0],Nb);
	ran = rand_uniform();

	clusters.resize(N);
	end_caps.resize(2*N);

	spins.resize(N);
	l_spins.resize(N);
	r_spins.resize(N);
	first.insert(first.end(),N,-1);
	last.insert(last.end(),N,-1);
	X.resize(H+4*N);
	opstr.resize(2*M);

	Vl = &end_caps[0];
	Vr = &end_caps[N];

	for(int i=0;i<N;i++){Vl[i] = Vr[i] = -1;}


	// create initial vbs state using greedy algorithm
	for(int i=0;i<N;i++){

		select_h.emplace_back(std::make_unique<select_random<double>>(&h[i*N],N));
		
		if(Vl[i]>=0){continue;} // bond already set

		std::vector<double> h_row(&h[i * N],&h[(i+1) * N]);
		std::vector<int> sites(N);
		std::iota(sites.begin(), sites.end(), 0);

		while(true){
			auto ele = std::max_element(h_row.begin(),h_row.end());
			auto site = (size_t)(ele - h_row.begin()) + sites.begin();
			const int j = *site;

			if(*ele > 0 && Vl[j] == -1){
			
				Vl[i] = Vr[i] = j;
				Vl[j] = Vr[j] = i;

				// spins on bond must be opposite signs
				l_spins[i] = r_spins[i] =  1;
				l_spins[j] = r_spins[j] = -1;

				break;
			}

			h_row.erase(ele);
			sites.erase(site);


		}
	}

	for(int p=0;p<2*M;p++){

		while(true){ // do not stop until bond is set
			const int ib = (*select_bond)();
			const int J_sign = (J_list[ib] > 0 ? 1 : -1);
			const int i = bond_list[2*ib];
			const int j = bond_list[2*ib+1];
			if(r_spins[i] * r_spins[j] * J_sign < 0){
				opstr[p] = 2*ib; 
				break;
			}
		}

	}
}


void proj_heis::print_opstr(){
	std::copy(r_spins.begin(),r_spins.end(),spins.begin());
	std::cout << "-----------config-----------" << std::endl;

	std::vector<char> bond_label(N);
	// print bonds by labeling pairs of spins
	char b = 'a';
	std::fill(bond_label.begin(),bond_label.end(),-1);

	for(int i=0;i<N;i++){
		const int j = Vr[i];
		if(bond_label[i]==-1 && bond_label[j]==-1){
			bond_label[i] = bond_label[j] = b++;
		}
	}

	std::cout << "     ";
	for(auto b : bond_label){
		std::cout << " " << b << " ";
	}
	std::cout << std::endl;
	// print out spins + operator configurations.
	for(int p=0;p<2*M;p++){
		std::cout << "     ";
		for(auto s: spins){
			std::cout << (s>0? " + " : " - ");
		}
		std::cout << std::endl;

		const int ib = opstr[p]/2;
		const int i = bond_list[2*ib];
		const int j = bond_list[2*ib+1];
		std::cout << "     ";
		if(opstr[p]%2){
			std::swap(spins[i],spins[j]);
			for(int s=0;s<N;s++){
				std::cout << ((s<i || s>j)? "   ": "===");
			}
		}			
		else{
			for(int s=0;s<N;s++){
				std::cout << ((s<i || s>j)? "   ": "---");
			}
		}
		std::cout << std::endl;

	}
	std::cout << "     ";
	for(auto s: spins){
		std::cout << (s>0? " + " : " - ");
	}
	std::cout << std::endl;


	// print bonds by labeling pairs of spins
	b = 'a';
	std::fill(bond_label.begin(),bond_label.end(),-1);

	for(int i=0;i<N;i++){
		const int j = Vl[i];
		if(bond_label[i]==-1 && bond_label[j]==-1){
			bond_label[i] = bond_label[j] = b++;
		}
	}
	std::cout << "     ";
	for(auto b : bond_label){
		std::cout << " " << b << " ";
	}
	std::cout << std::endl;
	std::cout << "-----------config-----------" << std::endl;
}

template<class callback>
void proj_heis::MCsweep(callback &cb){
	this->diag_update(cb);
	this->loop_update();
	this->ends_update();
}


void proj_heis::MCsweep(){
	auto cb = [](int p,signed char * s) {};
	this->diag_update(cb);
	this->loop_update();
	this->ends_update();
}

template<class callback>
void proj_heis::diag_update(callback &cb) {
	std::copy(r_spins.begin(),r_spins.end(),spins.begin());

	for(int p=0;p<2*M;p++){
		cb(p,&spins[0]);

		if(opstr[p]%2){
			const int ib = opstr[p]/2;
			const int i = bond_list[2*ib];
			const int j = bond_list[2*ib+1];

			std::swap(spins[i],spins[j]);
		}
		else{
			while(true){ // do not stop until bond is set
				const int ib = (*select_bond)();
				const int J_sign = (J_list[ib] > 0 ? 1 : -1);
				const int i = bond_list[2*ib];
				const int j = bond_list[2*ib+1];
				if(spins[i] * spins[j] != J_sign){
					opstr[p] = 2*ib; 
					break;
				}
			}			
		}
	}

	for(int i=0;i<N;i++){assert((spins[i]==l_spins[i]));}
	cb(2*M-1,&spins[0]);
}



void proj_heis::linked_list(){

	std::fill(last.begin(),last.end(),-1);
	std::fill(first.begin(),first.end(),-1);
	std::fill(X.begin(),X.end(),-1);

	for(int p=0;p<2*M;p++){
		const int ib = opstr[p]/2;
		const int v0 = 4*p;
		const int i1 = bond_list[2*ib];
		const int i2 = bond_list[2*ib+1];
		const int v1 = last[i1];
		const int v2 = last[i2];

		if(v1!=-1){X[v0]=v1;X[v1]=v0;}
		else{ first[i1] = v0;}

		if(v2!=-1){X[v0+1]=v2;X[v2]=v0+1;}
		else{first[i2] = v0+1;}

		last[i1] = v0+2;
		last[i2] = v0+3;
	}

	for(int i=0;i<N;i++){
		const int q0 = H+2*i;
		const int q1 = H+2*(N+i);
		if(first[i] == -1 && last[i] == -1){
			X[q0] = q1;	X[q1] = q0;
		}
		else{
			X[first[i]] = q0; X[q0] = first[i];
			X[last[i]] = q1; X[q1] = last[i];
		}
	}

	for(int i=0;i<N;i++){
		const int q0 = H+2*i+1;
		const int q1 = H+2*Vr[i]+1;
		const int q2 = H+2*(N+i)+1;
		const int q3 = H+2*(N+Vl[i])+1;

		X[q0] = q1; X[q1] = q0;
		X[q2] = q3; X[q3] = q2;
	}

}


void proj_heis::loop_update(){
	this->linked_list();

	std::fill(clusters.begin(),clusters.end(),-1);
	int cluster = 0;

	for(int v0=0;v0<X.size();v0+=2){

		if(X[v0]<0){continue;}
		int v = v0;
		bool found_mid = false;
		if(ran()<0.5){ // visit loop

			while(true){
				X[v] = -1;
				int p,p1,site;
				if(v < H){
					p = v/4;
					const int ib = opstr[p] / 2;
					
					if(J_list[ib] < 0){v = (v^2)^1;} // ferromagnetic: 0 <-> 3, 1 <-> 2
					else{v ^= 1;} // antiferromagnetic: 0 <-> 1, 2 <-> 3

					site = bond_list[2*ib+v%2];
				}
				else{
					v ^= 1;
					const int q = v-H;
					p = (q < 2*N ? -1 : 2*M);
					site = (q/2)%N;
				}

				const int v1 = X[v];
				if(v1 < H){
					p1 = v1/4;
				}
				else{
					p1 = ((v1-H) < 2*N ? -1 : 2*M);
				}

				if((p < M && p1 >= M ) || (p1 < M && p >= M)){
					clusters[site] = cluster;
					found_mid = true;
				}
				X[v] = -1; v = v1;
				if(v==v0) break;
			}
		}
		else{ // flip loop
			while(true){
				X[v] = -2;
				int p,p1,site;
				if(v < H){
					p = v/4;
					opstr[p] ^= 1;
					const int ib = opstr[p] / 2;
					
					if(J_list[ib] < 0){v = (v^2)^1;} // ferromagnetic: 0 <-> 3, 1 <-> 2
					else{v ^= 1;} // antiferromagnetic: 0 <-> 1, 2 <-> 3

					site = bond_list[2*ib+v%2];
				}
				else{
					v ^= 1;
					const int q = v-H;
					p = (q < 2*N ? -1 : 2*M);
					site = (q/2)%N;
				}

				const int v1 = X[v];
				if(v1 < H){
					p1 = v1/4;
				}
				else{
					p1 = ((v1-H) < 2*N ? -1 : 2*M);
				}

				if((p < M && p1 >= M ) || (p1 < M && p >= M)){
					clusters[site] = cluster;
					found_mid = true;
				}
				X[v] = -2; v = v1;
				if(v==v0) break;
			}
		}

		if(found_mid) cluster++;

	}



	for(int i=0;i<N;i++){
		const int q0 = H+2*i;
		const int q1 = H+2*(N+i);

		if(X[q0] == -2){r_spins[i] *= -1;}
		if(X[q1] == -2){l_spins[i] *= -1;}

	}


}




void proj_heis::ends_update(){
	//right side
	int count = 0;
	while(2*count < N){

		const int j0 = std::floor(N*ran());
		const int i0 = Vr[j0];
		const int j = (*select_h[i0])();
		if(r_spins[i0] == r_spins[j] || j==j0) continue;

		const int i = Vr[j];
		const int j1 = (*select_h[i])();

		if(j1!=j0) continue;

		Vr[i0] = j;
		Vr[j]  = i0;
		Vr[j0] = i;
		Vr[i]  = j0;

		count++;

	}

	count = 0;
	// left side
	while(2*count < N){

		const int j0 = std::floor(N*ran());
		const int i0 = Vl[j0];
		const int j = (*select_h[i0])();
		if(l_spins[i0] == l_spins[j] || j==j0) continue;

		const int i = Vl[j];
		const int j1 = (*select_h[i])();
		if(j1!=j0) continue;

		Vl[i0] = j;
		Vl[j]  = i0;
		Vl[j0] = i;
		Vl[i]  = j0;
		count++;
	}


}

#endif