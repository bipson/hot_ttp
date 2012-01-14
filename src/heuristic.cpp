#include <iostream>
#include <string>
#include <fstream>

#include <vector>
#include <map>

#include <time.h>

using namespace std;


#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)
#define ABS(a)    (a < 0 ? -a : a)

#define ALREADY_SET -1 // exact same spot has been taken
#define IS_DOUBLE -2   // same game has already been played
#define HOME_REPEATED_BEFORE -3 
#define HOME_REPEATED_AFTER -4
#define AWAY_REPEATED_BEFORE -5
#define AWAY_REPEATED_AFTER -6

#define CONSECUTIVE_HOME -7
#define CONSECUTIVE_AWAY -8

#define OK 0

#define IGNORE_CONSTRAINTS 1

//#define NEIGHBOURHOOD_OUTPUT_ITERATION


class Solution;

typedef vector<Solution> SolutionVector;

typedef SolutionVector * (Solution::*NeighbourStructure) (int k) ;
typedef Solution       * (* StepFunction)                (SolutionVector * solution_space, Solution * initial) ;

typedef pair<int, int> IntPair;

map<IntPair, double> pheromones;
map<IntPair, double>::const_iterator pheromones_iter;

double pheromones_factor = 100;

double pheromone_power  = 3;
double local_info_power = 1;

double pheromone_quotient_update_power = 6.0;

double t_min = 1.0;
double t_max = 30.0;

double pheromone_evaporation_rate = 1.05;
double pheromone_increase_rate    = 0.90;

int num_ants = 80;

class Solution
{
	int * plan;
	// plan[i] is an array representing list of teams playing on (i+1)th week
	// plan[i][j] therefore represents that team |j| (index |j-1|) plays on (i+1)th week against team |plan[i][j]| (index |plan[i][j]| - 1)
	// if plan[i][j] > 0, game is home, for negative, it is away. 0 means no game for team
	// in this array, teams are always referred to thru nonzero integers (1,2,...), but such integers actaually represent teams with indices (0, 1, ...)


public:

	static int * weights;
	static int dim;

	int is_valid; 

	const static int check_doubles				   = 1;
	const static int check_repeaters			   = 1;
	const static int check_consecutive_constraints = 1;

	const static int consecutive_games = 3;

	// set game on week, m_home against m_away (both m_ 1-indexed). 

	inline int CanSetGame(int week, int m_home, int m_away) const
	{
		if(GetGame(week, m_home) || GetGame(week, m_away)) 
			return ALREADY_SET; // game already set error

		if(check_doubles)
		{
			for(int i = 1; i <= 2 * (dim - 1); i++)
				if(GetGame(i, m_home) == m_away) return IS_DOUBLE;
		}

		if(check_repeaters)
		{
			if(week > 0 &&             (GetGame(week - 1, m_home) == -m_away)) return HOME_REPEATED_BEFORE; // home team, repeated before
			if(week < 2 * (dim - 1) && (GetGame(week + 1, m_home) == -m_away)) return HOME_REPEATED_AFTER; // home team, repeated after

			if(week > 0 &&             (GetGame(week - 1, m_away) ==  m_home)) return AWAY_REPEATED_BEFORE; // away team, repeated before
			if(week < 2 * (dim - 1) && (GetGame(week + 1, m_away) ==  m_home)) return AWAY_REPEATED_AFTER; // away team, repeated after
		}

		if(check_consecutive_constraints)
		{
			int sp = MAX(week - consecutive_games, 1);
			int ep = MIN(week + consecutive_games, 2 * (dim - 1) );
			
			int cons_t1 = 0;
			int cons_t2 = 0;

			for(int i = sp; i <= ep; i++)
			{
				if(i == week) 
				{
					cons_t1++;
					cons_t2++;
				}
				else 
				{
					int hg = GetGame(i, m_home);
					int ag = GetGame(i, m_away);

					cons_t1 = hg &&   IsHome(hg) ? cons_t1 + 1 : 0; // we want 0 if game is not yet set in both cases
					cons_t2 = ag &&  !IsHome(ag) ? cons_t2 + 1 : 0;

				}

				if     (cons_t1 > consecutive_games) return CONSECUTIVE_HOME; // fail, home team
				else if(cons_t2 > consecutive_games) return CONSECUTIVE_AWAY; // fail, away team
			}
		}

		return OK;
	}

	inline int SetGame(int week, int m_home, int m_away, int ignore_constraints = 0)
	{
		int csg = OK;
		
		if(ignore_constraints || (csg = CanSetGame(week, m_home, m_away)) == OK)
		{
			if(ignore_constraints) // if a game was set with ignore constraints, immediatelly assume its invalidity
				is_valid = 0;

			plan[dim * (week - 1) + m_home - 1] =   m_away;
			plan[dim * (week - 1) + m_away - 1] = - m_home;	
		}
	
		return csg;
	};
	inline int GetGame(int week, int team) const
	{
		return plan[dim * (week - 1) + team - 1];
	}

	inline int UnsetGame(int week, int m_home)
	{
		int m_away = ABS(plan[dim * (week - 1) + m_home - 1]);
		if(m_away != 0)
		{
			plan[dim * (week - 1) + m_home - 1] = 0;
			plan[dim * (week - 1) + m_away - 1] = 0;
			return 1;
		}
		else return 0;			
	}

	inline static int GetDistance(int m1, int m2) // can accept negative values
	{
		//return weights[ABS(m1) + ABS(m2) -2];
		return weights[(ABS(m1) - 1) * dim + ABS(m2) - 1];
	}

	inline int IsHome(int game) const
	{
		return game > 0;
	}
	inline int IsHome(int week, int team) const // answer for week/team pair directly rather than for getgame result
	{
		//getgame = return plan[dim * (week - 1) + team - 1];
		// ishome = getgame > 0
		//return plan[dim * (week - 1) + team - 1] > 0;
		return IsHome(GetGame(week, team));
	}

	inline int RespectsConstraints(int * error = NULL) // does this solution respect constraints? possible that it does not if plan was built up manually or with ignore constraints
	{
		Solution s;

		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			for(int m = 1; m <= dim; m++)
			{
				int g = GetGame(w, m); // try to set all games of this to s, if any screws up a constraint return 0
				if(IsHome(g)) // set only for home teams, this will be exactly half of the games, away will thus be set implicitly
				{
					int res = s.SetGame(w, m, g);

					if(res != OK)
					{
						is_valid = 0;
						if(error) *error = res;
						return 0;
					}
				}
			}
		}
		is_valid = 1;
		return 1; // everything was settable, return 1
	}

	// answer whether a plan is full (there is no empty spot)
	inline int IsFull()
	{
		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			for(int m = 1; m <= dim; m++)
			{
				if(GetGame(w, m) == 0) // one empty spot -> plan is not full
					return 0;
			}
		}
		return 1;
	}

	Solution()
	{
		is_valid = 1;
		plan = new int[dim * ( 2 * (dim - 1) )];
		memset(plan, 0, sizeof(int) * dim * ( 2 * (dim - 1) ));
	};
	Solution(const Solution &s)
	{
		this->plan = new int[dim * ( 2 * (dim - 1) )];
		memcpy(this->plan, s.plan, sizeof(int) *  dim * ( 2 * (dim - 1) ));
		this->is_valid = s.is_valid;
		//memset(plan, 0, sizeof(int) *  dim * ( 2 * (dim - 1) ));
	};

	int Evaluate() const
	{
		int dist = 0;

		for(int m = 1; m <= dim; m++)
		{
			int g = GetGame(1, m); // look at first week
			if(!IsHome(g)) // and for every away team
				dist += GetDistance(m, g); // count distance from its home to the location of the game

			    g = GetGame(2 * (dim - 1), m); // look at last week
			if(!IsHome(g))  // and for every away team
				dist += GetDistance(m, g); // count distance from its away location of the game to home
		}


		for(int w = 2; w <= 2 * (dim - 1); w++) // begin from 2nd week
		{
			for(int m = 1; m <= dim; m++) // for every team
			{
				int g = GetGame(w, m); // look at the game m plays on week w
				int game_before = GetGame(w-1, m); // look at the game it has played before

				if(!IsHome(g)) // we count only away team s(by convention)
				{
					int initial_pos = IsHome(game_before) ? m : abs(game_before); // team m's initial position is, if previous game was at home for m, at home (m),
																				  // otherwise, at team against which it has previously played
					dist += GetDistance(initial_pos, g);// g is away game, therefore g is negative, abs(g) == -g
				}
				else if(!IsHome(game_before))// they are home now, we must check whether they were home week before, if not,w e must add trip to home
				{
					dist += GetDistance(game_before, m); // from previous location to home location
				}
			}
		}

		return dist;
	}


	static int Greedy(Solution &s, int depth = 0)
	{
		int w = 0; // first week with blank spot (at least one)
		for(int wj = 1; wj <= 2 * (s.dim - 1); wj++)
		{
			int m = 0;
			for(int mj = 1; mj <= s.dim; mj++)
			{
				if(s.GetGame(wj, mj) == 0) 
				{
					m = mj;
					break;
				}
			}

			if(m) // a blank spot was found
			{
				w = wj;
				break;
			}
		}

		if(w != 0) // week with empty spot found, proceed with recursion
		{			
			vector<char> ban_list_home; // save space, use char
			vector<char> ban_list_away;

			while(true)
			{
				int home_candidate = 0; // home and away candidates
				int away_candidate = 0; // or both 0 if no applicable found

				double min_dist = numeric_limits<double>::max();

				// search for optimal (home, away) pair that will together travel shortest
				for(int m_home = 1; m_home <= s.dim; m_home++) // assume outer team is home team, assumption ok since another team will also get value of m_home
				{

					// get previous location for home team (assume it was at home)
					int m_home_previous_loc = m_home;
					if(w != 1) // on 0thweek  every team is at home
					{
						int pg = s.GetGame(w - 1, m_home); // get what happened week before
						if(!s.IsHome(pg))  // if team was away, calculate distance to home
							m_home_previous_loc = pg;
					}
					
					// calculate home distance
					int dist_home = GetDistance(m_home_previous_loc, m_home); // distance that must be travelled by home team


					// search for optimal away team that will travel the shortest
					for(int m_away = 1; m_away <= s.dim; m_away++) // for every other team
					{
						if(m_away == m_home) 
							continue;

						if(s.GetGame(w, m_home) || s.GetGame(w, m_away)) // if any of the teams already has a game ignore the pair
							continue;

						// check for ban
						int is_on_ban_list = 0;
						for(int bi = 0; bi < ban_list_home.size(); bi++)
						{
							if(ban_list_home[bi] == m_home && ban_list_away[bi] == m_away)
							{
								is_on_ban_list = 1;
								break;
							}
						}
						if(is_on_ban_list) 
							continue;
						// end ban check


						// get previous location for away team (assume it was away)
						int m_away_previous_loc = m_away;
						if(w != 1) // on 0thweek  every team is at home
						{
							int pg = s.GetGame(w - 1, m_away); // get what happened week before
							if(!s.IsHome(pg))  // if team was away, calculate distance to home
								m_away_previous_loc = pg;
						}

						// calculate away distance
						int dist_away = GetDistance(m_away_previous_loc, m_home); // distance that must be travelled by away team


						double local_info = dist_home + dist_away;

						IntPair pheromone_key;
						pheromone_key.first  = m_home;
						pheromone_key.second = m_away;

						double dist = pow(local_info, local_info_power) + pheromones_factor * pow(pheromones[pheromone_key], pheromone_power); // consists of local info + pheromone info				
						
						//dist = pow(local_info, local_info_power) ; // consists of local info + pheromone info
						//cout << pow(local_info, local_info_power) << ", " << pow(pheromones[pheromone_key], pheromone_power) << endl;
						// if this is minimum, and this game can be set
						if(dist < min_dist && s.CanSetGame(w, m_home, m_away) == OK)
						{	
							min_dist = dist;
							home_candidate = m_home;
							away_candidate = m_away;
						}
					}
				}

				if(home_candidate && away_candidate) // if a suitable candidate is found
				{
					s.SetGame(w, home_candidate, away_candidate);
				
					// recursive call
					if(Greedy(s, depth + 1)) // recursive call success
					{
						return 1;
					}
					else // failure, ban current candidate
					{
						s.UnsetGame(w, home_candidate);
						//ban_list_home[ban_list_size  ] = home_candidate;
						//ban_list_away[ban_list_size++] = away_candidate;
						ban_list_home.push_back(home_candidate);
						ban_list_away.push_back(away_candidate);
					}
				}
				else // when no suitable candidate found, return failure
					return 0;
			}
		}
		else // no empty spot, correct solution found
		{
			return 1; 
		}
	};

	static Solution Flip_Game(const Solution *s, int m1, int m2, int ignore_constraints = 0)
	{
		int w1 = 0, w2 = 0; 
		for(int w = 1; w <= 2 * (dim - 1); w++) // optimization idea - construct, for a solution, a lookup table that would, for every team combo, put weeks in which tey play
		{
			if(s->GetGame(w, m1) == m2) w1 = w; // find the week of the game in which m1 plays against m2, m1 home
			if(s->GetGame(w, m2) == m1) w2 = w; // find the week of the game in which m1 plays against m2, m2 home
		}	

		Solution new_solution (*s); // copy this solution
						
		// unset the games
		new_solution.UnsetGame(w1, m1);
		new_solution.UnsetGame(w2, m2);

		//set flipped, if allowed
		//if(new_solution.SetGame(w1, m2, m1, ignore_constraints) == OK && new_solution.SetGame(w2, m1, m2, ignore_constraints) == OK)
		if(new_solution.SetGame(w1, m2, m1, ignore_constraints) != OK || new_solution.SetGame(w2, m1, m2, ignore_constraints) != OK)
			new_solution.is_valid = 0;
		
		return new_solution;
	}


	SolutionVector * Neighbourhood_k_Flip(int k)
	{
		SolutionVector *ss = new SolutionVector(0); // TODO set proper limit

		int a_limit = -1; // first element for k - 1
		int b_limit = -1; // last element for k- 1 th iteration

		for(int i = 1; i <= k; i++)
		{
			#ifdef NEIGHBOURHOOD_OUTPUT_ITERATION
			cout << "k-Flip, " << i << "th iteration" << endl;
			#endif

			for(int el = a_limit; el <= b_limit; el++)
			{
				for(int m1 = 1; m1 <= dim; m1++)// select 2 teams, find weeks they play in
				{
					for(int m2 = m1 + 1; m2 <= dim; m2++) // only look over/right the diagonal (from 1,1 to n,n)
					{
						// first time, both limits -1, s = this
						Solution *s = (a_limit == -1 ? this : &(*ss)[el]); // must be inside these loops, since it points to location in ss, which can be moved as ss grows						

						Solution ns = Flip_Game(s, m1, m2); // try to flip
						if(ns.is_valid)
						{
							// unique filter
							int already_exists = 0;
							for(int j = 0; j < ss->size(); j++) // check whether it exists
							{
								if(ns == (*ss)[j])
								{
									already_exists = 1;
									break;
								}
							}
							if(!already_exists)
								ss->push_back(ns);
						}
					}
				}
			}

			a_limit = b_limit    + 1;
			b_limit = ss->size() - 1;
		}
		
		return ss;
	};

	SolutionVector * Neighbourhood_k_Opt(int k)
	{
		SolutionVector ss;

		ss.push_back(*this); // add self to neighbourhood

		int a_limit = 0; // first element for k - 1
		int b_limit = 0; // last  element for k - 1 th iteration

		for(int i = 1; i <= k; i++)
		{

			#ifdef NEIGHBOURHOOD_OUTPUT_ITERATION
			cout << "k-Opt, " << i << "th iteration" << endl;
			#endif

			for(int el = a_limit; el <= b_limit; el++)
			{	
				for(int w1 = 1; w1 < 2 * (dim - 1); w1++)
				{
					for(int w2 = w1 + 1; w2 <= 2 * (dim - 1); w2++) // pick another, different week
					{
						Solution s(ss[el]);

						int * w2_games = new int[dim];
						memcpy(w2_games,		        s.plan + dim * (w2 - 1), sizeof(int) * dim); // save w2th row
						memcpy(s.plan + dim * (w2 - 1), s.plan + dim * (w1 - 1), sizeof(int) * dim); // w1th row -> w2th row
						memcpy(s.plan + dim * (w1 - 1), w2_games,				 sizeof(int) * dim); // w2th row (saved) -> w1th row
						delete w2_games;

						// unique filter
						int already_exists = 0;
						for(int j = 0; j < ss.size(); j++) // check whether it exists
						{
							if(s == ss[j])
							{
								already_exists = 1;
								break;
							}
						}
						if(!already_exists)
							ss.push_back(s);
					}
				}
			}

			a_limit = b_limit   + 1;
			b_limit = ss.size() - 1;
		}
		
		SolutionVector * proper_solution_vector = new SolutionVector;
		
		for(int i = 0; i < ss.size(); i++) // set i=1 instead of i=0 to ignore self 
		{
			if(ss[i].RespectsConstraints()) // if an element within a list respects constraints, put it in proper solutions list
			{
				//if(ss[i] > *this) // ignore worse solutions?
					proper_solution_vector->push_back(ss[i]);
			}
		}
		
		return proper_solution_vector;
	};

	static void DisplayWeights()
	{
		for(int i = 1; i <= dim; i++)
		{
			for(int j = 1; j <= dim; j++)
			{
				//cout << weights[dim * i + j] << "\t ";
				cout << GetDistance(i, j) << "\t";
			}
			cout << endl;
		}
	};

	void DisplaySolutionMatrix(int show_symmetry = 1, int show_team_names = 0) const
	{
		// initialize team names
		char * team_names[12];
		{
			char tni = 0;
			team_names[tni++] = "ATL";
			team_names[tni++] = "NYM";
			team_names[tni++] = "PHI";
			team_names[tni++] = "MON";
			team_names[tni++] = "FLA";
			team_names[tni++] = "PIT";
			team_names[tni++] = "CIN";
			team_names[tni++] = "CHI";
			team_names[tni++] = "STL";
			team_names[tni++] = "MIL";	
			team_names[tni++] = "HOU";
			team_names[tni++] = "COL";
		}


		cout << "\nSolution matrix: \n W/M\t";
		// mannschaft nummer
		for(int m = 1; m <= dim; m++)
		{
			if(show_team_names)
				cout << team_names[m-1] << "\t";
			else
				cout << m << "\t";
		}
		cout << "\n";

		if(show_team_names)
		{
			for(int m = 1; m <= dim + 1; m++)
				cout << "---\t";
		}
		cout << endl;

		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			cout << w << "\t";
			for(int m = 1; m <= dim; m++)
			{
				int el = GetGame(w, m);

				if (el == 0) cout << "-";
				else
				{
					if(el < 0 && show_symmetry) cout << "@";

					if(el > 0 || show_symmetry)		  
					{
						if(show_team_names)
							cout << team_names[ABS(el) - 1];
						else cout << el;
					}
				}

				cout << "\t";
			}
			cout << "\n";
		}
		cout << endl;
	}


	void UpdatePheromones(double amount)
	{
		IntPair key;

		for(int m = 1; m <= dim; m++)
		{
			key.first = m;

			int * schedule_for_m = new int[2 * (dim  - 1)];
			
			for(int w = 1; w <= 2 * (dim - 1); w++)
				schedule_for_m[w - 1] = GetGame(w, m);
						
			for(int w = 1; w <= 2 * (dim - 1); w++)
			{
				key.second = schedule_for_m[w - 1];

				IntPair k2;
				k2.first = key.first;
				k2.second = key.second;

				double val = pheromones[key] * amount * pheromone_increase_rate;				
				 // val = pheromones[key];
				pheromones[key] = val;
				
				if(pheromones[key] < t_min)
					pheromones[key] = t_min;

				key.first = key.second; 
			}

			delete schedule_for_m;
		}
	}

	bool operator==(const Solution &other)
	{
		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			for(int m = 1; m <= dim; m++) 
			{
				if(this->IsHome(w, m) && this->GetGame(w, m) != other.GetGame(w, m))
					return 0;
			}
		}
		return 1;
	}
	bool operator!=(const Solution &other)
	{
		return !(*this == other);
	}
	// take into account validity?
	bool operator>(const Solution &other)
	{
		return this->Evaluate() < other.Evaluate();
	}
	bool operator>=(const Solution &other)
	{
		return this->Evaluate() <= other.Evaluate();
	}
	bool operator<(const Solution &other)
	{
		return !(*this >= other);
	}
};

int Solution::dim;
int * Solution::weights;

class Loader
{
	string str;
	int dim;

public:
	int * data;

	Loader(int _dim)
	{
		dim = _dim;
		data = new int[dim * dim];
	};

	int Load(string filename)
	{
		char * buffer = 0;
		long length;
		FILE * f = fopen (filename.c_str(), "rb");

		if (f)
		{
		  fseek (f, 0, SEEK_END);
		  length = ftell (f);
		  fseek (f, 0, SEEK_SET);
		  buffer = (char *) malloc (length + 1);
		  if (buffer)
		  {
			fread (buffer, 1, length, f);
		  }
		  buffer[length] = 0;
		  fclose (f);
		}

		if (buffer)
		{			
			int i = 0, j = 0;

			char *p = strtok(buffer, " \t\n");

			while (p) 
			{
				if(i == dim) break;
				else if (j == dim) 
				{
					j = 0;
					i++;
				}	

				data[dim * i + j++] = atol(p);		

				p = strtok(NULL, " \t\n");	
			}
			
		}

		//int 
		
		Solution::dim = dim;
		Solution::weights = data;
		
		return 0;
	};
};


int local_search_ctr = 0;
int vnd_ctr = 0;
int gvns_ctr = 0;

Solution * LocalSearch(Solution *initial, 
						NeighbourStructure neighbourhood_structure,
						int neighbourhood_structure_parameter, 
						StepFunction step_function
						)
{
	// neighbourhood_structure = NeighbourhoodStructure(int k)
	// step_function		   = StepFunction          (SolutionVector *neighbourhood, Solution * initial)

	Solution *current = initial;

	int ctr = 0;

	int iterations_without_improvement = 0;

	SolutionVector * solution_space = NULL;

	while(iterations_without_improvement < 10)
	{
		ctr++;

		if(solution_space == NULL)
			solution_space = (current->*neighbourhood_structure)(neighbourhood_structure_parameter);
		
		Solution	   * next = step_function (solution_space, current);

		if(next != NULL && *next > *current)
		{
			current = next;

			//cout << "Local search (Iteration " << ctr << "), distance: " << current->Evaluate() << ". Solution space size: " << solution_space->size() << ".\n";
			//current->DisplaySolutionMatrix();

			solution_space = NULL; // we must calculate new solution space
			iterations_without_improvement = 0;
		}
		else 
		{
			iterations_without_improvement++;
			//cout << "No neighbour found. " << iterations_without_improvement << " such iteration. " << endl;
		}
	}

	local_search_ctr += ctr;

	return current;
}

Solution * NextImprovement(SolutionVector *neighbourhood, Solution * initial)
{
	for(int i = 0; i < neighbourhood->size(); i++)
	{
		if(neighbourhood->at(i) > *initial)
			return neighbourhood->data() + i;
	}
	//return initial;
	return NULL; // null pointer on failure
};

Solution * BestImprovement(SolutionVector *neighbourhood, Solution * initial)
{
	Solution *best = initial;
		
	for(int i = 0; i < neighbourhood->size(); i++)
	{
		if(neighbourhood->at(i) > *best)
			best = neighbourhood->data() + i;
	}
	return best == initial ? NULL : best;
};

Solution * RandomNeighbour(SolutionVector *ss, Solution * initial__unused = NULL)
{
	return ss->data() + rand() % ss->size();
};

Solution * VND(Solution *initial, NeighbourStructure * neighbourhood_structures, int *neighbour_structures_params, StepFunction step_function, int num_neighbourhood_structures = 2) // Variable Neighborhood Descent
{
	int l = 0;
	int ctr = 0;

	Solution current = *initial;

	while(l < num_neighbourhood_structures)
	{
		cout.flush();

		SolutionVector * solspace = (current.*neighbourhood_structures[l])(neighbour_structures_params[l]);

		Solution *p_candidate = step_function(solspace, &current);

		ctr++;

		if(p_candidate != NULL && *p_candidate > current)
		{
			current = *p_candidate;
			l = 0;

			//cout << "    VND (Iteration " << ctr << "), distance: " << current.Evaluate() << ". Solution space size: " << solspace->size() << "."  << endl;
			//current.DisplaySolutionMatrix();
		}
		else 
		{
			if(p_candidate != NULL)
				delete p_candidate;
			//cout << "    VND (Iteration " << ctr << "). No viable candidate found (considered one with distance " << current.Evaluate() << ". Solution space size: " << solspace->size() << "." << endl;
			l++;
		}
		if(solspace)
			delete solspace;
	}

	vnd_ctr += ctr;

	return new Solution(current);
};

// entscheiden zwischen folgendem
Solution * GRASP() // Greedy Randomized Adaptive Search Procedure.
{
	return NULL;
};

// Generalized Variable Neighborhood Search
Solution * GVNS(Solution *initial, NeighbourStructure * neighbourhood_structures, int *neighbour_structures_params, StepFunction step_function, int num_neighbourhood_structures = 2) 
{
	Solution current = *initial;

	int iterations = 0;

	int ctr = 0;

	while(++iterations <= 10)
	{
		int l = 0;

		while(l < num_neighbourhood_structures)
		{
			//cout << "\n\n\n\nGVNS (Iteration " << ++ctr << " Begin)\n" << endl;

			SolutionVector * solspace = (current.*neighbourhood_structures[l])(neighbour_structures_params[l]);

			Solution * p_candidate = VND(
														RandomNeighbour(solspace),  // shaking - generate x shaken
														neighbourhood_structures, 
														neighbour_structures_params, 
														step_function, 
														num_neighbourhood_structures);

			//cout << "Back in GVNS..\n";

			if(*p_candidate > current)
			{
				current = *p_candidate;
				l = 0;

				//cout << "\nGVNS (Iteration " << ctr << "), distance: " << current.Evaluate() << ". Solution space size: " << solspace->size() << "." << endl;
				//current.DisplaySolutionMatrix();
			}
			else
			{
				//cout << "\nGVNS (Iteration " << ctr << "). No viable candidate found (considered one with distance " << current.Evaluate() << ". Solution space size: " << solspace->size() << ".\n" << endl;
				l++;
			}

			if(solspace)
				delete solspace;
		}
	}

	gvns_ctr += ctr;

	return new Solution(current);
};

int ex1(int argc, char ** argv)
{
	/* // sample solution output
	{

	int dim = 8;
	char filename[25] = "";
	sprintf(filename, "data%d.txt", dim);

	Loader l(dim);
	
	l.Load(filename);

	Solution s;

	int w = 1;
	
	s.SetGame(w,   1, 5);
	s.SetGame(w,   4, 2);
	s.SetGame(w,   6, 8);
	s.SetGame(w++, 7, 3);
	
	s.SetGame(w,   2, 8);
	s.SetGame(w,   3, 5);
	s.SetGame(w,   4, 1);
	s.SetGame(w++, 6, 7);
	
	s.SetGame(w,   2, 1);
	s.SetGame(w,   3, 8);
	s.SetGame(w,   4, 7);
	s.SetGame(w++, 5, 6);
	
	s.SetGame(w,   1, 6);
	s.SetGame(w,   5, 3);
	s.SetGame(w,   7, 2);
	s.SetGame(w++, 8, 4);
	
	s.SetGame(w,   1, 3);
	s.SetGame(w,   5, 4);
	s.SetGame(w,   7, 6);
	s.SetGame(w++, 8, 2);
	
	s.SetGame(w,   1, 4);
	s.SetGame(w,   6, 2);
	s.SetGame(w,   7, 5);
	s.SetGame(w++, 8, 3);
	
	s.SetGame(w,   2, 5);
	s.SetGame(w,   3, 7);
	s.SetGame(w,   4, 8);
	s.SetGame(w++, 6, 1);
	
	s.SetGame(w,   2, 7);
	s.SetGame(w,   3, 1);
	s.SetGame(w,   4, 5);
	s.SetGame(w++, 8, 6);
	
	s.SetGame(w,   3, 2);
	s.SetGame(w,   5, 8);
	s.SetGame(w,   6, 4);
	s.SetGame(w++, 7, 1);
	
	s.SetGame(w,   1, 8);
	s.SetGame(w,   5, 2);
	s.SetGame(w,   6, 3);
	s.SetGame(w++, 7, 4);
	
	s.SetGame(w,   1, 2);
	s.SetGame(w,   4, 3);
	s.SetGame(w,   6, 5);
	s.SetGame(w++, 7, 8);
	
	s.SetGame(w,   1, 7);
	s.SetGame(w,   2, 4);
	s.SetGame(w,   3, 6);
	s.SetGame(w++, 8, 5);
	
	s.SetGame(w,   2, 6);
	s.SetGame(w,   3, 4);
	s.SetGame(w,   5, 7);
	s.SetGame(w++, 8, 1);

	s.SetGame(w,   2, 3);
	s.SetGame(w,   4, 6);
	s.SetGame(w,   5, 1);
	s.SetGame(w++, 8, 7);
	
	s.DisplaySolutionMatrix(1, 1);
	cout << s.Evaluate(); 
	system("pause");
	return 0;}
	*/
	// cmd line params
	// 1: dim
	// 2: k-flip k parameter
	// 3: k-opt k parameter
	// 4: 1 for only g, l for only local search benchmark

	int dim = 4; // default 8
	char output_filename[64] = ""; // blank string -> output to console

	int k_flip_k = 4;
	int k_opt_k  = 2;

	int make_gvns = 1;
	int bmark_ls  = 1;

	int start_time = time(NULL);

	switch(argc)
	{
		case 5:
			if(argv[4][0] == 'g')
				bmark_ls = 0;
			else if(argv[4][0] == 'l')
				make_gvns = 0;

		case 4:
			if(! (k_opt_k = atoi(argv[3]))) 
				k_opt_k = 2; // otherwise stick to default
			
		case 3:
			if(! (k_flip_k = atoi(argv[2]))) 
				k_flip_k = 4; // otherwise stick to default
			
		case 2: // only dim given
			if(! (dim = atoi(argv[1]))) 
				dim = 8; // otherwise stick to default
			break;
	}	
	

	streambuf *psbuf, *cout_backup;
	ofstream filestr;

	sprintf(output_filename, "%d_dim_%d_Flip_%d_Opt%s%s.txt", dim, k_flip_k, k_opt_k, make_gvns ? "_GVNS" : "", bmark_ls ? "_LS" : "");
	

	if(*output_filename)
	{
		filestr.open (output_filename);
		cout_backup = cout.rdbuf();     // back up cout's streambuf
		psbuf = filestr.rdbuf();   // get file's streambuf
		cout.rdbuf(psbuf);         // assign streambuf to cout
	}

	srand(start_time);

	char filename[25] = "";
	sprintf(filename, "data%d.txt", dim);

	Loader l(dim);
	
	l.Load(filename);

	Solution greedy_solution;

	
	int i = 0;

	// comments - set initial game
	{
		//s.SetGame(1, 1, 2);
		//s.SetGame(2, 2, 1);

		/*
		// example game load
		s.SetGame(++i, 1, 3);
		s.SetGame(  i, 2, 4);

		s.SetGame(++i, 1, 2);
		s.SetGame(  i, 3, 4);

		s.SetGame(++i, 1, 4);
		s.SetGame(  i, 3, 2);

	
		s.SetGame(++i, 3, 1);
		s.SetGame(  i, 4, 2);

		s.SetGame(++i, 2, 1);
		s.SetGame(  i, 4, 3);

		s.SetGame(++i, 4, 1);
		s.SetGame(  i, 2, 3);
		//*/
	}

	Solution::Greedy(greedy_solution);
	

	cout << "Input parameters: \nDim: " << dim << "\n" << k_flip_k << "-Flip\n" << k_opt_k << "-Opt\n" << "Output filename: " << output_filename << "\n\n";

	greedy_solution.DisplaySolutionMatrix(0);
	
	cout << "Above: initial solution. Distance: "<< greedy_solution.Evaluate() <<"\n\n";

	NeighbourStructure neighbourhood_structures[2] =
	{
		&Solution::Neighbourhood_k_Flip,
		&Solution::Neighbourhood_k_Opt,
	};
	int neighbour_structures_params[2] = {k_flip_k, k_opt_k};
	StepFunction step_function = BestImprovement;

	Solution * s;
	
	if(make_gvns)
	{
		s = GVNS(&greedy_solution, neighbourhood_structures, neighbour_structures_params, step_function);

		cout << "\n\n\nFinal solution. Distance: " << s->Evaluate() << "." << endl;
		s->DisplaySolutionMatrix(0);

		cout << "GVNS iterations: "			<< gvns_ctr			<< "\n"
			 << "VND iterations: "			<< vnd_ctr			<< "\n";
			 //<< "Local Search iterations: " << local_search_ctr << endl;
	
		cout << "\nRun Time: " << time(NULL) - start_time << " s" << endl;
	}

	if(!bmark_ls)
	{
		cout << endl;
		filestr.flush();
		filestr.close();
		return 0;
	}

	start_time = time(NULL); // reset ctr

	cout << "\n\n\n";

	cout << "Local Search Benchmarks (initial solution distance: " << greedy_solution.Evaluate() << ")\n";

	s = LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Flip, k_flip_k, NextImprovement);
	cout << "Next Improvement, " << k_flip_k << "-Flip. Distance: " << s->Evaluate() << ". Run Time: " << time(NULL) - start_time << " s\n";
	start_time = time(NULL); // reset ctr

	s = LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Opt,  k_opt_k,  NextImprovement);
	cout << "Next Improvement, " << k_opt_k << "-Opt. Distance: " << s->Evaluate()   << ". Run Time: " << time(NULL) - start_time << " s\n\n";
	start_time = time(NULL); // reset ctr

	s = LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Flip, k_flip_k, BestImprovement);
	cout << "Best Improvement, " << k_flip_k << "-Flip. Distance: " << s->Evaluate()  << ". Run Time: " << time(NULL) - start_time << " s\n";
	start_time = time(NULL); // reset ctr

	s = LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Opt,  k_opt_k,  BestImprovement);
	cout << "Best Improvement, " << k_opt_k << "-Opt. Distance: " << s->Evaluate()   << ". Run Time: " << time(NULL) - start_time << " s\n\n";
	start_time = time(NULL); // reset ctr
		

	const int num_iterations = 30;

	double avg_distance = 0;
	for(int i = 0; i < num_iterations; i++)
		avg_distance += LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Flip, k_flip_k, RandomNeighbour)->Evaluate() / num_iterations;
	cout << "Random Neighbour, " << k_flip_k << "-Flip. Distance: " << avg_distance << ". Run Time: " << ((time(NULL) - start_time)) / num_iterations << " s\n";
	start_time = time(NULL); // reset ctr

	avg_distance = 0;
	for(int i = 0; i < num_iterations; i++)
		avg_distance += LocalSearch(&greedy_solution, &Solution::Neighbourhood_k_Opt, k_opt_k, RandomNeighbour)->Evaluate() / num_iterations;
	cout << "Random Neighbour, " << k_opt_k << "-Opt. Distance: " << avg_distance << ". Run Time: " << (time(NULL) - start_time) / num_iterations << " s\n\n";
	start_time = time(NULL); // reset ctr

	s = VND(&greedy_solution, neighbourhood_structures, neighbour_structures_params, BestImprovement);
	cout << "VND, " << k_flip_k << "-Flip. " << k_opt_k << "-Opt. Distance: " << s->Evaluate() << ". Run Time: " << time(NULL) - start_time << " s\n\n";

	cout << endl;

	//*/

	filestr.close();
	if(!*output_filename)
		system("pause");

	return 0;
}

int main(int argc, char ** argv)
{
	int dim = 8;
	
	int start_time = time(NULL);
	//srand(start_time);

	char filename[25] = "";
	sprintf(filename, "data%d.txt", dim);

	Loader l(dim);
	
	l.Load(filename);

	
	pair<int, int> key;
	
	for(key.first = -dim; key.first <= dim; key.first++)
	{
		for(key.second = -dim; key.second <= dim; key.second++)
		{
			if(key.second == 0) continue;

			pheromones[key] = t_min * pheromone_evaporation_rate; // init dict, all poossible edges
		}
	}


	int k_flip_k = 2;
	int k_opt_k  = 1;
	NeighbourStructure neighbourhood_structures[2] =
	{
		&Solution::Neighbourhood_k_Flip,
		&Solution::Neighbourhood_k_Opt,
	};
	int neighbour_structures_params[2] = {k_flip_k, k_opt_k};
	StepFunction step_function = BestImprovement;


	Solution best_solution;
	Solution::Greedy(best_solution);
	
	
	int best_sol_val = best_solution.Evaluate();
	int best_ant_index = -1;

	double initial_sol_val = best_sol_val;

	for(int ant = 1; ant <= num_ants; ant++)
	{
		Solution s; // to every ant his own
		Solution::Greedy(s);
		
		Solution * s_after_gvns = GVNS(&s, neighbourhood_structures, neighbour_structures_params, step_function);

		//s.UpdatePheromones(pow( s.Evaluate() / initial_sol_val, pheromone_quotient_update_power));
		s_after_gvns->UpdatePheromones(pow( s.Evaluate() / initial_sol_val, pheromone_quotient_update_power));

		cout << ant << ": " << s_after_gvns->Evaluate() << endl;

		if(best_solution < *s_after_gvns)
		{			
			best_solution = *s_after_gvns;
			best_sol_val  = s_after_gvns->Evaluate();
			best_ant_index = ant;
		}

		delete s_after_gvns;

		//*
		for (pheromones_iter=pheromones.begin(); pheromones_iter != pheromones.end(); ++pheromones_iter) 
		{
			pheromones[pheromones_iter->first] *= pheromone_evaporation_rate; // evaporate pheromones
			
			if(pheromones[pheromones_iter->first] > t_max)
				pheromones[pheromones_iter->first]  = t_max;
		}
		//*/
	}


	best_solution.DisplaySolutionMatrix();

	cout << endl << endl << best_solution.Evaluate() << endl;
	
	cout << "Best ant index: " << best_ant_index << endl;

	cout << "Run Time: " << time(NULL) - start_time << " s" << endl;

	system("pause");

	return 0;
}
