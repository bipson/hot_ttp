#Ant Colony Optimization: An important issue for an effective ACO is a meaningful pheromone model.
#Therefore, you have to find a compromise between a complex but strong model and an efficient, compact
#one (concerning memory usage and speed). Another important factor is the use of heuristic information
#during the construction process and the use of local search.
#Therefore, you should implement an ACO variant of your own choice for which you design at least two
#ant construction heuristics and incorporate one or more local search procedures from exercise one.
#
#Test your approach with the instances mentioned above. Do not only optimize your algorithm just for some
#special cases but also for general use. Therefore, you should test with other instances, too. Again, your program
#must create a protocol file, which contains the best solutions and their objective values. You must also log
#the average objective value of each iteration.
#Again, write an abstract of two pages and bring a pri_nted version to the interview. You should discuss your
#algorithm and the problems you had during the implementation. The abstract must include a table which
#documents your results. Compare the different operators/heuristics you used and document your results within
#the abstract. Keep in mind that you should take the average value of multiple test runs to get a meaningful
#result.

#Double round robin (A at B and B at A): n teams need 2*n-2 slots.
#No more than three consecutive home or three consecutive road games for any team
#No repeaters (A at B, followed immediately by B at A)
#Objective is to minimize distance traveled (assume teams begin in their home city and must return there after the tournament)
#City names are in the order ATL NYM PHI MON FLA PIT CIN CHI STL MIL HOU COL SF SD LA ARI where NLx takes the first x cities.
#Data files contain the distance matrix (square, symmetric).

import math
import random
import itertools
import copy
from collections import defaultdict

dim = 4

class Game:
    def __init__(self, m_home, m_away):
        self.m_home = m_home
        self.m_away = m_away
    def __repr__(self):
        return "(" + str(self.m_away) + " @ " + str(self.m_home) + ")"
    def copy(self):
        return Game(self.m_home, self.m_away)
    def invert(self): #inverts home and away team
        return Game(self.m_away, self.m_home)
    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.m_home == other.m_home and self.m_away == other.m_away
    def __hash__(self):
        return hash((self.m_home,self.m_away))

class Distances:
    def __init__(self):
        global dim
        data = []
        with open("data/data" + str(dim) + ".txt", "r") as f:
            for line in f.read().strip().splitlines():
                data.append(list(map(int, line.strip().split())))
        
        for i in range(dim): # quick validation
            if int(data[i][i]) != 0:
                raise Exception("Distance to self must be 0 - failure at " + str(i))
            for j in range(dim):
                if int(data[i][j]) != int(data[j][i]):
                    raise Exception("Distance matrix must be symmetrical - failure at " + str((i, j)))
        self.__data = data
    def get(self, m1, m2):
        return self.__data[m1 - 1][m2 - 1]
    
    def __repr__(self):
        return "\n".join(
                         ["\t".join(
                                  [self.get(i + 1, j + 1) for j in range(dim)]
                        )  for i in range(dim)]
                ) 
            

class Solution:
    distances = Distances()
    consecutive_games = 3
    
    class FailureType:
        ALREADY_SET = 0
        IS_DOUBLE = 1
        #HOME_REPEATED_BEFORE = 2
        #HOME_REPEATED_AFTER = 3
        #AWAY_REPEATED_BEFORE = 4
        #AWAY_REPEATED_AFTER = 5
        REPEATED_BEFORE = 4
        REPEATED_AFTER = 5
        CONSECUTIVE_HOME = 6
        CONSECUTIVE_AWAY = 7
    
    
    def __init__(self):        
        global dim
        self.dim = dim   
        self.plan = [set() for i in range(2 * (dim - 1))] # 2 * dim * (dim - 1) weeks, initialize them
        self.is_valid     = True  # are there any violations of constraints (blank schedule has no violations)
        self.is_complete  = False # are there any empty spots in the schedule
        self.ban_list     = defaultdict(lambda: set([])) # set of banned games for a given week - this is a dict - thus 1! indexed
        self.game_assign_sequence = [] # lists the sequence of assignments of games (week, game)
        
    def copy(self):
        s = Solution()
        s.plan = copy.deepcopy(self.plan)
        (s.is_valid, s.is_complete) = (self.is_valid, self.is_complete)
        s.ban_list     = copy.copy(self.ban_list)
        s.game_assign_sequence = self.game_assign_sequence[:]
        return s
    
    # get all games that were set, regardless of the week
    def get_all_games(self):
        all_games = set()
        for game_list in self.plan:
            all_games = all_games.union(game_list)
        return all_games
    
    def get_games(self, week): # week (and team) indices start at 1
        return self.plan[week - 1]
    
    def get_game(self, week, team):
        for game in self.get_games(week):
            if game.m_home == team or game.m_away == team:
                return game
        return None # if nothing found, return none
    
    def can_set_game(self, week, game = None): 
        if game == None:
            game = week
            week = self.get_first_empty_week_index()
               
        failure_types = set() # assume all is well
        
        # game already setta
        if self.get_game(week, game.m_home) or self.get_game(week, game.m_away):
            failure_types.add(self.FailureType.ALREADY_SET)
            
        # doubles - check if the exact game was already set
        if game in self.get_all_games():
            failure_types.add(self.FailureType.IS_DOUBLE) 
        
         
        # repeaters - fail if the same game, with different teams was played immediately afterwards or before
        if week > 0             and (game.invert() in self.get_games(week - 1)): failure_types.add(self.FailureType.REPEATED_BEFORE) # repeated before
        if week < 2 * (dim - 1) and (game.invert() in self.get_games(week + 1)): failure_types.add(self.FailureType.REPEATED_AFTER)  # repeated after

    
        # check_consecutive_constraints
        sp = max(week - self.consecutive_games, 1); # top and bottom border
        ep = min(week + self.consecutive_games, 2 * (dim - 1) );
        cons_t1 = 0
        cons_t2 = 0

        for i in range(sp, ep + 1): # go thru possible affected weeks
            if i == week:           # if at current week, increase both counters
                cons_t1 += 1
                cons_t2 += 1
            else:
                hg = self.get_game(i, game.m_home); # get the games for home and away team, played at ith week
                ag = self.get_game(i, game.m_away);

                cons_t1 = (hg and hg.m_home == game.m_home) and cons_t1 + 1 or 0
                cons_t2 = (ag and ag.m_away == game.m_away) and cons_t2 + 1 or 0

            
            if cons_t1 > self.consecutive_games:
                failure_types.add(self.FailureType.CONSECUTIVE_HOME)  
                return False #  fail, home team
            if cons_t2 > self.consecutive_games:
                failure_types.add(self.FailureType.CONSECUTIVE_AWAY) #  fail, away team
        
        # here, log failure types if required
        if len(failure_types) > 0:
            pass
            #rint(failure_types)
        
        return len(failure_types) == 0
    
    #tries to set game, if successful return true, otherwise false
    # if only one param given, it is a game - week is current week
    def set_game(self, week, game = None):
        if game == None:
            game = week
            week = self.get_first_empty_week_index()
            
        if self.can_set_game(week, game):
            self.get_games(week).add(game.copy())
            
            self.game_assign_sequence.append((week, game.copy())) # memorize the last decision - for backtracking
            
            self.is_complete = True # assume it is true, if a non-full week is found, set to false
            for week in self.plan:
                if len(week) < dim / 2: # a week is full if it has exactly dim/ 2 games
                    self.is_complete = False
                    break
            
            return True
        else:
            return False
        
    def unset_last_game(self):
        if len(self.game_assign_sequence) == 1:
            self.ban_list[0].add(self.game_assign_sequence[0][1])
            self.ban_list[1] = set(self.ban_list[0])
            #self.game_assign_sequence = []
            #self.is_complete = False # no more complete
            #return
        
        week, game = self.game_assign_sequence.pop()
        
        # if the week is full, we must delete all entries in the ban list for the following weeks
        if len(self.get_games(week)) == self.dim / 2:
            #for w in range(week + 1, 2 * (self.dim - 1)):
            # can only be the next week
            self.ban_list[week + 1] = set()
            
        self.get_games(week).remove(game) #simply remove the last decision
        self.ban_list[week].add(game) # and memorize that it was removed
        
        self.is_complete = False # no more complete
        # TODO !!! set is_valid ?
        
    
    def __add__(self, game):
        s = self.copy()
        success = s.set_game(game)
        return success and s or None
    
    def __iadd__(self, game):
        return s.set_game(game) and self or None
    
    # cost of adding a game (in terms of distance)
    def game_cost(self, week, game = None):  
        if game == None:
            game = week
            week = self.get_first_empty_week_index()
        
        if not self.can_set_game(week, game):
            return None
        
        if week == 1: # 1st week, special calculation
            return self.distances.get(game.m_away, game.m_home) # evaluate for the current game 
        
        # dist that must be travelled by home team + dist that must be travelled by the away team
        return self.distances.get(self.get_game(week - 1, game.m_home).m_home, game.m_home) + \
               self.distances.get(self.get_game(week - 1, game.m_away).m_home, game.m_home)  
        
    def disp(self, show_symmetry = True, show_team_names = False):        
        s = "Solution matrix \n W/M\t"
        
        if show_team_names:
            team_names = "ATL NYM PHI MON FLA PIT CIN CHI STL MIL HOU COL SF SD LA ARI".split(" ")
            s += "\t".join(team_names[:dim])
        else:
            s += "\t".join(map(str, range(1, dim + 1)))
        
        s += "\n" + "---\t" * (dim + 1)
        
        for w in range(2 * (dim - 1)):
            s += "\n" + str(w + 1) + "\t"    
            for m in range(dim):
                game = self.get_game(w + 1, m + 1) 
                                
                if not game: 
                    s += "-" # if tehre is no game output a -
                else:
                    if m + 1 == game.m_home: # whether the current team (m) is away
                        s += " " + (show_team_names and team_names[game.m_away - 1] or str(game.m_away))
                    elif show_symmetry:
                        s += "@" + (show_team_names and team_names[game.m_home - 1] or str(game.m_home))
                s += "\t"
        return s
    
    def __repr__(self):
        return self.disp()         
    
    def __eq__(self, other):
        return self.plan == other.plan
    
    def evaluate(self):
        # solution must be valid and complete in order to be evaluable
        if not self.is_valid or not self.is_complete:
            return None
        
        dist = 0
        # a row in a plan is a set of games

        # look at the first week, there, every away team has to travel from its home(m_away) to the location of the game (m_home) 
        # and, after the last game, every away team has to travel back home
        for week in [1, 2 * (self.dim - 1)]: # get first and last week
            dist += sum(map(lambda game: self.distances.get(game.m_away, game.m_home), self.get_games(week)))

        
        # for every week except first, observe from where the team has travelled (where it was the week before)
        
        for week in range(2, 2 * (self.dim - 1) + 1):
            for game in self.get_games(week):
                # for both home and away team, check where they were previously - they had to travel from there
                # if a home team was previously at home, distance will be 0
                # note: every match is played at home team's location
                dist += self.distances.get(self.get_game(week - 1, game.m_home).m_home, game.m_home)
                dist += self.distances.get(self.get_game(week - 1, game.m_away).m_home, game.m_home)

        return dist
    
    def __lt__(self, other): return self.evaluate() < other.evaluate()
    def __gt__(self, other): return self.evaluate() > other.evaluate()

    def get_first_empty_week_index(self):
        for week in range(1, 2 * (self.dim)):
            if len(self.get_games(week)) < dim / 2:
                break
        return week
    
    def get_possible_decisions(self):
        first_empty_week = self.get_first_empty_week_index()        
        # get all the teams that are not playing on first empty week

        remaining_teams = set(range(1, self.dim + 1)).difference(*list(map(lambda g: set([g.m_home, g.m_away]), self.get_games(first_empty_week))))
        
        # possible games are pairs of teams - use combinations of 2 
        
        # exclude from possible games those in the tabu list for the given week
        possibilities = []
        for tg in itertools.combinations(remaining_teams, 2):
            possibilities.append(Game(*tg))
            possibilities.append(Game(*tg).invert())
            
        possibilities = set(possibilities).difference(self.ban_list[first_empty_week]).difference(self.get_all_games())
        
        possibilities = set(filter(lambda x: self.can_set_game(x), possibilities))
        
        #if first_empty_week == 1 and len(possibilities) == 0: # no rescue here -
        #    return None  
        
        return possibilities
    
    def get_stench(self, pheromones, game): # return how much would adding a certain game stinks
        # assume possible to set game
        week_before = ((self.get_first_empty_week_index()- 2) % 6) + 1
        
        if self.get_games(week_before) == set():
            return t_max
        else:
            t1 = (game.m_away, game.m_home, self.get_game(week_before, game.m_away).m_home)
            t2 = (game.m_home, game.m_home, self.get_game(week_before, game.m_home).m_home)
            #t_both = (t1, t2)
            #return pheromones[t_both]
            return (pheromones[t1] + pheromones[t2] ) / 2
    
    def update_pheromones(self, pheromones):
        for week in range(1, 2 * (self.dim) - 1):
            week_before = ((week-2) % 6) + 1
                        
            for g in self.get_games(week):
                t1 = (g.m_away, g.m_home, self.get_game(week_before, g.m_away).m_home)
                t2 = (g.m_home, g.m_home, self.get_game(week_before, g.m_home).m_home)
                t_both = (t1, t2)
                pheromones[t1] *= 0.2
                pheromones[t2] *= 0.2
                #pheromones[t_both] *= 0.1
                
    def __hash__(self):
        return hash(tuple(map(tuple, self.plan)))
        
            

random.seed()

iterations = 0

stopping_criteria = lambda: iterations > 50


(t_min, t_max) = (1.0, 10.0)
(stench_power, local_info_power) = 8, 2
pheromones = defaultdict(lambda: t_max)

num_ants = 10

while not stopping_criteria():
    
    best_solution, best_solution_value = None, 99999999
     
    for i in range(num_ants):
        
        s = Solution() # to every ant his own
        while not s.is_complete:
            decisions_list = []
            
            #if (s.get_first_empty_week_index()) == 2:
                #print (s.get_possible_decisions())
            #print("1:",s.ban_list[1],"\n2:",s.ban_list[2])
            #print(s.plan)
            total_prob_value = 0
            possible_decisions = s.get_possible_decisions()
            
            if possible_decisions == None: # this ant really sucks, try another one (not able to proceed, not even with backtracking)
                break
                # if this happens, it will kill evaluate below
            elif len(possible_decisions) > 0:
                for g in possible_decisions:
                    probability = math.pow(s.get_stench(pheromones, g), stench_power) * math.pow(1000/s.game_cost(g), local_info_power)
                    
                    total_prob_value += probability 
                    decisions_list.append((probability, g))
            else: # reverse last decision, add it to tabu list (automatic), and try again
                s.unset_last_game()
                continue
            
            random.shuffle(decisions_list)
            random_value = random.random() * total_prob_value
            
            decision = None
            for random_decrementer, d in decisions_list:
                random_value -= random_decrementer
                if random_value <= 0:
                     decision = d
                     break
            
            s += decision
        
        val = s.evaluate() # if possible to evaluate solution
        if val != None and val < best_solution_value: # update bestw
            best_solution, best_solution_value = s, val
            print(best_solution)
            print(val)
    
    if best_solution.is_complete:        
        #pheromones[best_solution] += 1000/best_solution_value # update pheromones
        #pheromones[best_solution] += pheromones[best_solution]/10 # update pheromones
        best_solution.update_pheromones(pheromones)
    
    for i in pheromones: # evaporate pheromones
        pheromones[i] *= 0.9
    
    iterations += 1
    
print(pheromones)
     
#===============================================================================
# 
# s = Solution()
# #// example game load
# 
# 
# s.set_game(  1, Game(1, 3))
# s.set_game(  1, Game(2, 4))
# 
# s.set_game(  2, Game(1, 2))
# s.set_game(  2, Game(3, 4))
# 
# s.set_game(  3, Game(1, 4))
# s.set_game(  3, Game(3, 2))
# 
# 
# s.set_game(  4, Game(3, 1))
# s.set_game(  4, Game(4, 2))
# 
# s.set_game(  5, Game(2, 1))
# s.set_game(  5, Game(4, 3))
# 
# s.set_game(  6, Game(4, 1))
# s.set_game(  6, Game(2, 3))
# 
# print (s,s.evaluate())
#===============================================================================
