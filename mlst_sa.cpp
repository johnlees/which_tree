#include <Rcpp.h>
#include <RInside.h>
#include <string>
#include <vector>
#include <random>

// Namespaces
using namespace std;
using namespace Rcpp;

// Helper R script
string r_script_name("mlst_init.R");

int total_genes = 2175;
int mlst_scheme_genes = 7;

// Simulated annealing parameters
// Currently use linear temp decrease
// TODO consider implementing temp as class
long int states = 1000000;
long int sample_rate = 1000;
long int final_temp = 0;

default_random_engine generator;
uniform_int_distribution<int> mlst_index(1,mlst_scheme_genes);
uniform_int_distribution<int> rand_gene(1,total_genes);

// MLST class contains current vector and distance. Call update to propose
// a new value, and accept it based on the current temperature
class MLST
{
   public:
      // Initialisation
      MLST(const vector<int>& initial);

      // Non-modifying
      double distance() const { return _distance; }
      NumericVector genes() const { return _mlst; }

      // Modifying
      int get_distance(string) { return _R.parseEval("tree_metric(" + string + ")"); }

      // Propose and evaluate function value
      int update(double temperature);

   private:
      vector<int> _mlst;
      double _distance;
      RInside _R;
};

void printInfo(const MLST&, double temperature, int acceptances, int iterations);

MLST::MLST(const vector<int>& initial)
{
   // Set up R environment
   RInside R(argc, argv);
   R.parseEvalQ("source(" + r_script_name + ")");

   MLST::_R = R;

   // Initialise mlst vec

   // Calculate distance for this value
}

// Propose, evaluate and accept based on temperature
// Returns 1 if accepted, 0 if rejected
int MLST::update(double temperature)
{
   int accepted = 0;

   // Proposal
   int replace_idx = mlst_index(generator);
   int replaced_gene = _mlst[replace_idx];

   int new_gene = rand_gene(generator);
   while (new_gene == replaced_gene)
   {
      new_gene = rand_gene(generator);
   }

   _mlst[replace_idx] = new_gene;
   std::sort(_mlst.begin(), _mlst.end());

   // Calculate distance
   string mlst_genes
   double new_dist = _R.parseEval("tree_metric()");

   // Decide whether to accept
   if (new_dist < _distance)
   {
      _distance = new_dist;
      accepted = 1;
   }

   return accepted;
}

void printInfo(const MLST&, double temperature, int acceptances, int iterations)
{
   double acceptance_rate = acceptances / (double) iterations;
   cout << acceptance_rate << "\t" << temperature << endl;
   cout << state.genes() << "\t" << state.distance() << endl;
}

int main(int argc, char *argv[])
{
   // Set an initial value
   vector start = 100, 200, 300, 400, 500, 600, 700;
   MLST state(start);

   long int acceptances, iterations = 0;
   long int temperature = states;

   // Main loop
   while (temperature > final_temp)
   {
      // Update state
      if (state.update())
      {
         acceptances++;
      }

      // Update temperature
      temperature--;
      iterations++;

      if (temperature % sample_rate)
      {
         printInfo(state, temperature, acceptances, iterations);
      }
   }

   // Diagnostics
   printInfo(state, temperature, acceptances, iterations);
}
