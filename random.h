#ifndef _random_h_
#define _random_h_

#include <random>
#include <vector>
#include <stack>


class rand_uniform
{
	std::mt19937_64 gen;
	std::uniform_real_distribution<double> dist;
	unsigned int s;

	public:


		rand_uniform(){
			//seeding random number generator
			unsigned int lo,hi;
			__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
			s=((unsigned long long)hi << 32) | lo;

			gen.seed(s);
			dist = std::uniform_real_distribution<double>(0.0,1.0);
		};

		inline double operator()(void){
			return dist(gen);
		}

		inline
		unsigned int seed(void) const {return s;}

};



template<class T>
class select_random
{
	const size_t N;
	rand_uniform ran;

	std::vector<T> P;
	std::vector<size_t> A;

public:
	select_random(const T * weights,const size_t N_): N(N_){
		ran = rand_uniform();
		P.reserve(N_);
		A.reserve(N_);
		A.insert(A.end(),N,-1);

		std::stack<size_t> under,over;


		T norm = 0;
		for(size_t i=0;i<N_;i++){
			norm += weights[i];
		}

		norm /= N_;

		for(size_t i=0;i<N_;i++){
			const T p = weights[i] / norm;
			if(p > 1){
				over.push(i);
			}
			else if (p < 1){
				under.push(i);
			}

			P.push_back(p);
		}


		while(!(over.empty() || under.empty())){
			const size_t i = over.top(); over.pop();
			const size_t j = under.top(); under.pop();
			
			const T res = P[i] + P[j] - 1.0;

			P[i] = res;
			A[j] = i;

			if(res > 1){
				over.push(i);
			}
			else if (res < 1){
				under.push(i);
			}
		}
	}

	inline
	size_t operator()(){
		const double q = N * ran();
		const size_t i = (size_t)std::floor(q);
		if((q-i) < P[i]){
			return i;
		}
		else{
			return A[i];
		}

	}


};


#endif