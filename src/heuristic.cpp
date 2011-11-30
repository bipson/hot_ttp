#include <iostream>
#include <string>
#include <fstream>

#include <vector>

#ifdef linux
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <string.h>
#endif /*linux*/

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


class Solution;

typedef vector<Solution> SolutionVector;

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

	const static int check_doubles				   = 1;
	const static int check_repeaters			   = 1;
	const static int check_consecutive_constraints = 1;

	const static int consecutive_games = 3;

	// set game on week, m_home against m_away (both m_ 1-indexed). 

	inline int CanSetGame(int week, int m_home, int m_away)
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
			plan[dim * (week - 1) + m_home - 1] =   m_away;
			plan[dim * (week - 1) + m_away - 1] = - m_home;	
		}
	
		return csg;
	};
	inline int GetGame(int week, int team)
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

	inline int IsHome(int game)
	{
		return game > 0;
	}
	inline int IsHome(int week, int team) // answer for week/team pair directly rather than for getgame result
	{
		return IsHome(GetGame(week, team));
	}

	inline int RespectsConstraints() // does this solution respect constraints? possible that it does not if plan was built up manually or with ignore constraints
	{
		Solution s;

		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			for(int m = 1; m <= dim; m++)
			{
				int g = GetGame(w, m); // try to set all games of this to s, if any screws up a constraint return 0
				if(IsHome(g)) // set only for home teams, this will be exactly half of the games, away will thus be set implicitly
				{
					if(!s.SetGame(w, m, g))
					{
						return 0;
					}
				}
			}
		}

		return 1; // everything was settable, return 1
	}

	Solution()
	{
		plan = new int[dim * ( 2 * (dim - 1) )];
		memset(plan, 0, sizeof(int) * dim * ( 2 * (dim - 1) ));
	};
	Solution(const Solution &s)
	{
		this->plan = new int[dim * ( 2 * (dim - 1) )];
		memcpy(this->plan, s.plan, sizeof(int) *  dim * ( 2 * (dim - 1) ));
		//memset(plan, 0, sizeof(int) *  dim * ( 2 * (dim - 1) ));
	};

	int Evaluate()
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
			int   ban_list_size = 0;
			int * ban_list_home = new int[dim * dim];
			int * ban_list_away = new int[dim * dim];

			memset(ban_list_home, 0, sizeof(int) * dim * dim);
			memset(ban_list_away, 0, sizeof(int) * dim * dim);

			while(true)
			{
				int home_candidate = 0; // home and away candidates
				int away_candidate = 0; // or both 0 if no applicable found

				int min_dist = INT_MAX;

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
						for(int bi = 0; bi < ban_list_size; bi++)
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

						// if this is minimum, and this game can be set
						if((dist_home + dist_away) < min_dist && s.CanSetGame(w, m_home, m_away) == OK)
						{	
							min_dist = dist_home + dist_away;
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
						ban_list_home[ban_list_size  ] = home_candidate;
						ban_list_away[ban_list_size++] = away_candidate;
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

	static Solution * Flip_Game(Solution *s, int m1, int m2, int ignore_constraints = 0)
	{
		int w1 = 0; // find the week of the game in which m1 plays against m2, m1 home
		for(int w = 1; w <= 2 * (dim - 1); w++) // optimization idea - construct, for a solution, a lookup table that would, for every team combo, put weeks in which tey play
			if(s->GetGame(w, m1) == m2) w1 = w;

		int w2 = 0;
		for(int w = 1; w < 2 * (dim - 1); w++)    // find the week of the game in which m1 plays against m2, m2 home
			if(s->GetGame(w, m2) == m1) w2 = w;

		Solution * new_solution = new Solution(*s); // copy this solution
						
		// unset the games
		new_solution->UnsetGame(w1, m1);
		new_solution->UnsetGame(w2, m2);

		//set flipped, if allowed
		if(new_solution->SetGame(w1, m2, m1, ignore_constraints) == OK && new_solution->SetGame(w2, m1, m2, ignore_constraints) == OK)
			return new_solution;
		else
			return NULL;
	}


	SolutionVector * Neighbourhood_k_Flip(int k)
	{
		SolutionVector *ss = new SolutionVector(0); // TODO set proper limit

		int a_limit = -1; // first element for k - 1
		int b_limit = -1; // last element for k- 1 th iteration

		for(int i = 1; i <= k; i++)
		{
			for(int el = a_limit; el <= b_limit; el++)
			{
				for(int m1 = 1; m1 <= dim; m1++)// select 2 teams, find weeks they play in
				{
					for(int m2 = m1 + 1; m2 <= dim; m2++) // only look over/right the diagonal (from 1,1 to n,n)
					{
						// first time, both limits -1, s = this
						Solution *s = (a_limit == -1 ? this : &(*ss)[el]); // must be inside these loops, since it points to location in ss, which can be moved as ss grows						

						Solution * ns = Flip_Game(s, m1, m2); // try to flip
						if(ns != NULL)
							ss->push_back(*ns); // unique filter here? 
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
		SolutionVector ss; // TODO set proper limit

		ss.push_back(*this);

		int a_limit = 0; // first element for k - 1
		int b_limit = 0; // last  element for k - 1 th iteration

		for(int i = 1; i <= k; i++)
		{
			for(int el = a_limit; el <= b_limit; el++)
			{
				// pick a week
				for(int w = 1; w <= 2 * (dim - 1); w++)
				{
					for(int m1 = 1; m1 <= dim; m1++) // pick m1
					{
						if(!IsHome(w, m1)) 
							continue;

						for(int m2 = m1 + 1; m2 <= dim; m2++) // pick m2 diagonally
						{
							if(!IsHome(w, m2)) 
								continue;

							// pick 2 different teams, both of which play at home
							// this will give us 2 pairs of teams, m1 and its opponent, and m2 and its opponent
							// we know that m1 does not play against m2 in week w, since both play at home
							
							Solution s(ss[el]); // copy corresponding solution

							// get opponents of both teams, we have pairs (m1, m1_op) @ m1, and (m2, m2_op) @ m2
							int m1_op = s.GetGame(w, m1);
							int m2_op = s.GetGame(w, m2);
							
							s.UnsetGame(w, m1); // unset these two games
							s.UnsetGame(w, m2);
							
							s.SetGame(w, m1, m2_op, IGNORE_CONSTRAINTS); // set that each plays against the former opponent of the other
							s.SetGame(w, m2, m1_op, IGNORE_CONSTRAINTS);
							
							ss.push_back(s); // put the solution into solutions list

							/*
							Solution * flip_one = Flip_Game(&s, m1, m2_op, IGNORE_CONSTRAINTS); // flip one of these matches, and add it
							ss.push_back(*flip_one);

							ss.push_back(*Flip_Game(flip_one, m2, m1_op, IGNORE_CONSTRAINTS)); // flip another one along with the first

							ss.push_back(*Flip_Game(&s, m2, m1_op, IGNORE_CONSTRAINTS)); // now flip 2nd one, without the fisrt
							//*/
						}
					}
				}
			}

			a_limit = b_limit   + 1;
			b_limit = ss.size() - 1;
		}
		
		SolutionVector * proper_solution_vector = new SolutionVector;
		
		for(int i = 0; i < ss.size(); i++)
		{
			if(ss[i].RespectsConstraints()) // if an element within a list respects constraints, put it in proper solutions list
			{
				proper_solution_vector->push_back(ss[i]);
			}
		}

		return proper_solution_vector;
	};

	Solution * NextImprovement(SolutionVector ss)
	{
		return ss.data();
	};

	Solution * BestImprovement(SolutionVector ss)
	{
		int min = INT_MAX;
		int cand = 0;

		for(int i = 0; i < ss.size(); i++)
		{
			int val = ss[i].Evaluate();

			if(val < min)
			{
				min = val;
				cand = i;
			}
		}
		return &ss[cand];
	};

	Solution * RandomNeighbour(SolutionVector ss)
	{
		return ss.data() + rand() % ss.size();
	};

	// entscheiden zwischen folgendem
	Solution VND() // Variable Neighborhood Descent
	{
	};

	Solution GRASP() // Greedy Randomized Adaptive Search Procedure.
	{
	};

	Solution GVNS() // Generalized Variable Neighborhood Search
	{
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

	void DisplaySolutionMatrix(int show_symmetry = 1, int show_team_names = 0)
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


		cout << "\nSolution matrix: \n\n W/M\t";
		// mannschaft nummer
		for(int m = 1; m <= dim; m++)
			cout << m << "\t";
		cout << endl;

		for(int w = 1; w <= 2 * (dim - 1); w++)
		{
			cout << w << "\t";
			for(int m = 1; m <= dim; m++)
			{
				int el = GetGame(w, m);

				if(el > 0)		  
				{
					if(show_team_names)
						cout << team_names[el - 1];
					else cout << el;
				}
				else if (el == 0) cout << "-";
				else if(show_symmetry) cout << "@" << -el;

				cout << "\t";
			}
			cout << endl;
		}
		cout << endl;

	}
};

typedef vector<Solution> SolutionVector;

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

int main()
{
	int dim = 4;

	Loader l(dim);
	
	l.Load("data4.txt");

	//Solution::DisplayWeights();
	
	Solution s;

	
	int i = 0;



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

	Solution::Greedy(s);

	//SolutionVector *ss = s.Neighbourhood_k_Flip(1);
	SolutionVector *ss = s.Neighbourhood_k_Opt(1);
	
	for(int i = 0; i < ss->size(); i++)
	{
		(*ss)[i].DisplaySolutionMatrix(0, 1);
		cout << (*ss)[i].Evaluate() << endl;
	}
	cout << "Total number: " << ss->size();
	//s.DisplaySolutionMatrix(0);
	//cout << s.Evaluate() << endl;

	
	
	//*/

	system("pause");

	return 0;
}
