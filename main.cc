#include <iostream>
#include <math.h>
#include <cstdlib>
#include <time.h>

using namespace std;

// board size NxN
int N = 8;

// initial population size
int PopSize;

double SurvivalRate = 0.05;

// determines to percentage of random queen positions in a board to alter
double MutationRate;

class Board {
	public:
		int* queens;

//------------------------------------Board---------------------------------------
//
//                  Constructor : creates array with queens
//--------------------------------------------------------------------------------
		Board() {
			queens = new int[N];
            srand((unsigned)time(0));
				for(int j=0; j<N; j++) {
				queens[j]=(rand()%N);
			}
		}

//------------------------------------Print---------------------------------------
//
//                      This function prints board.
//--------------------------------------------------------------------------------
		void Print() {
			for(int i=0; i<N; i++) {
				for(int j=0; j<queens[i]; j++) {
					cout << "| ";
				}
				cout << "|Q";
				for(int j=0; j<N-queens[i]-1; j++) {
					cout << "| ";
				}
				cout << "|\n";
			}
			cout << "Fitness: " << Fitness() << "\n\n";
		}

//---------------------------------- Fitness --------------------------------------
//
//                        This function counts fitness.
//----------------------------------------------------------------------------------
		int Fitness() {
			int colCount[N];
			int upperDiagCount[(2*N)-1];
			int lowerDiagCount[(2*N)-1];

            for(int i=0; i<N; i++)
				colCount[i] = 0;
            for(int i=0; i<(2*N-1); i++){
				upperDiagCount[i] = 0;
				lowerDiagCount[i] = 0;
            }

			for(int i=0; i<N; i++) {
				colCount[queens[i]] += 1;
				upperDiagCount[queens[i]+i]+=1;
				lowerDiagCount[(N-queens[i])+i-1]+=1;
			}
			int fitness = 0;
			for(int i=0; i<2*N-1; i++) {
				if(i<N) {
					fitness += ((colCount[i]-1)*colCount[i])/2;
				}
				fitness += ((upperDiagCount[i]-1)*upperDiagCount[i])/2;
				fitness += ((lowerDiagCount[i]-1)*lowerDiagCount[i])/2;
			}
			return fitness;
		}

//------------------------------------Mutate---------------------------------------
//
//      This function mutates a population's bits dependent on the mutation rate
//----------------------------------------------------------------------------------
		void Mutate() {
            int mutationCount = max(min((int)floor(MutationRate*N), N), 1);
			for(int i=0; i<mutationCount; i++) {
				//random position
				int j=(rand()%N);
				//random value
				char v =(rand()%N);
				queens[j]=v;
			}
}

//---------------------------------- Mate ---------------------------------------
//
//            This function selects a random point along the length
//         of the population and swaps all the  bits after that point.
//                Here a pair of parents produces one children.
//-------------------------------------------------------------------------------
        Board Mate(Board Parent) {
			Board Child;
			// random position between 1 and N-1
			int r=(rand()%(N-1))+1;
			for(int i=0; i<r; i++) {
				Child.queens[i]=queens[i];
			}
			for(int i=r; i<N; i++) {
				Child.queens[i]=Parent.queens[i];
			}
			return Child;
        }

};

// the population of boards
Board* Population;

//---------------------------------- Initialization -------------------------------------
//
//                             This function creates populations
//----------------------------------------------------------------------------------------
void Initialization() {
	Population = new Board[PopSize];
	for(int i=0; i<PopSize; i++) {
		Population[i] = Board();
	}
}

//--------------------------------- sortPopulation ---------------------------------------
//
//            This function sorts populations by fitness in decreasing order
//                            using modified quicksort.
//----------------------------------------------------------------------------------------
void sortPopulation(int Left, int Right) {
	int i = Left, j = Right;
	Board tmp;
	int pivot = Population[(Left+Right)/2].Fitness();
	while(i <= j) {
		while(Population[i].Fitness() < pivot)
			i++;
		while(Population[j].Fitness() > pivot)
			j--;
		if(i <= j) {
			tmp = Population[i];
			Population[i] = Population[j];
			Population[j] = tmp;
			i++;
			j--;
		}
	};

	if(Left < j)
		sortPopulation(Left, j);
	if(i < Right)
		sortPopulation(i, Right);
}

//--------------------------------Elitist-------------------------------------------
//
//       If the best individual from the new population is better than
//          the best individual from the previous population, then
//           copy the best from the new population; else replace the
//            worst individual from the current population with the
//              best one from the previous generation.
//----------------------------------------------------------------------------------
void Elitist(int *bestCon, int best[])
{
  if ( Population[0].Fitness() >= *bestCon ){
    for(int i=0; i<N;i++)
        Population[PopSize-1].queens[i] = best[i];
  }
  else{
    *bestCon=Population[0].Fitness();
    for(int i=0; i<N;i++)
        best[i] = Population[0].queens[i];
  }
}

//--------------------------------------- GA ---------------------------------------
//
//           This function calls all genetic operations in order one by one:
//                    selection, crossover and mutation.
//              Checks fitness after each round of operations.
//------------------------------------------------------------------------------------
int GA() {
	int numberOfsteps=0;
    int best[N], best1[N];
	int lastFitness = -1,bestCon=-1, bestCon1=-1;
	while(true) {
        //values are needed for elitist selection
        if(numberOfsteps%2!=0 and bestCon1>lastFitness)
                bestCon1=lastFitness;
        if(numberOfsteps%2==0 and bestCon>lastFitness)
                bestCon=lastFitness;

		// sort population by fitness level
		sortPopulation(0, PopSize-1);
		int smallestFitness = Population[0].Fitness();
		if(lastFitness != smallestFitness) {
			lastFitness = smallestFitness;
			cout << "Fitness: " << smallestFitness << endl;
		}
		if(smallestFitness == 0) {
			break;
		}

		// crossover
		int cut = min(PopSize - 1, max((int)floor(PopSize*SurvivalRate), 1));
		int numberOfKids = floor(1/SurvivalRate)-1;
		for(int i=cut, j=0; i<PopSize; i+=numberOfKids, j++) {
			for(int k=0; k<numberOfKids; k++) {
				Population[i+k]=Population[j].Mate(Population[(j+k)%cut]);
			}
		}

        //values are needed for elitist selection
        if(numberOfsteps==0){
            bestCon=lastFitness;
            bestCon1=lastFitness;
            for(int i=0; i<N;i++){
                best1[i]=Population[0].queens[i];
                best[i]=Population[0].queens[i];
            }
        }
        if(numberOfsteps%2==0 and bestCon1==lastFitness){
            for(int i=0; i<N;i++)
                best[i]=Population[0].queens[i];
		}
		if(numberOfsteps%2!=0 and bestCon==lastFitness){
            for(int i=0; i<N;i++)
                best1[i]=Population[0].queens[i];
        }

        //Elitist selection
        if(numberOfsteps>0){
            if(numberOfsteps%2!=0)
                Elitist(&bestCon1, best);
            else
                Elitist(&bestCon, best1);
		}

        // mutate every board with the mutation rate
        for(int i=0; i<PopSize; i++) {
			Population[i].Mutate();
    }
		// now that we have a new generation increment the number of moves
		numberOfsteps++;
	}

return numberOfsteps;
}

//-------------------------------------------------------------------------------------
//                                      Main
//-------------------------------------------------------------------------------------
int main() {

	PopSize = 10000;
	MutationRate = 0.01;

	cout << "Number of queens: " << N << endl;
	cout << "Population size: " << PopSize << endl;
	cout << "Mutation rate: " << MutationRate << endl;

	Initialization();

	cout << "Number of Moves: " << GA() << endl;

	Population[0].Print();

	return 0;
}
