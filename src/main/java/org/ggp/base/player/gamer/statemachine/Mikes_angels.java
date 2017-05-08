package org.ggp.base.player.gamer.statemachine;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import org.ggp.base.apps.player.Player;
import org.ggp.base.player.gamer.exception.GamePreviewException;
import org.ggp.base.util.game.Game;
import org.ggp.base.util.statemachine.MachineState;
import org.ggp.base.util.statemachine.Move;
import org.ggp.base.util.statemachine.Role;
import org.ggp.base.util.statemachine.StateMachine;
import org.ggp.base.util.statemachine.cache.CachedStateMachine;
import org.ggp.base.util.statemachine.exceptions.GoalDefinitionException;
import org.ggp.base.util.statemachine.exceptions.MoveDefinitionException;
import org.ggp.base.util.statemachine.exceptions.TransitionDefinitionException;
import org.ggp.base.util.statemachine.implementation.prover.ProverStateMachine;


public class Mikes_angels extends StateMachineGamer {

	//Limit for depth of our recursive search.
	private static final int LIMIT = 12;
	private static final int DEPTH_CHARGES = 5;
	Player p;
	Boolean ranOutOfTime;
	double mobHeuristic;
	double goalHeuristic;

	@Override
	public StateMachine getInitialStateMachine() {
		// TODO Auto-generated method stub
		return new CachedStateMachine(new ProverStateMachine());
	}

	@Override
	public void stateMachineMetaGame(long timeout)
			throws TransitionDefinitionException, MoveDefinitionException,
			GoalDefinitionException {

//		int numGames = 0;
//		goalHeuristic = 0.5;
//		mobHeuristic = 0.5;
//
//		List<Double> mobVector = new ArrayList<Double>();
//		List<Double> goalVector = new ArrayList<Double>();
//		List<Double> terminalVector = new ArrayList<Double>();
//
//		while(true){
//			List<List<Double>> result = runRandomGame(timeout);
//			if(result == null){
//				System.out.println("Out of time");
//				break;
//			}
//			mobVector.addAll(result.get(0));
//			goalVector.addAll(result.get(1));
//			terminalVector.addAll(result.get(2));
//
//			numGames++;
//		}
//
//		//If we dont play more than one game, we only have one data point, so its meaningless
//		if(numGames > 1){
//			//System.out.println(mobVector);
//			//System.out.println(goalVector);
//			//System.out.println(terminalVector);
//
//			//X is val of mobility or goalutility, y is the final state
//			double B1 = covar(mobVector, mean(mobVector), terminalVector, mean(terminalVector)) / var(mobVector, mean(mobVector));
//			double B0 = mean(terminalVector) - B1*mean(mobVector);
//
//			double C1 = covar(goalVector, mean(goalVector), terminalVector, mean(terminalVector)) / var(goalVector, mean(goalVector));
//			double C0 = mean(terminalVector) - C1*mean(goalVector);
//
//			//The first val is for mobility
//			//The 2nd val is for goal utility
//			System.out.println("Mobility: " + B0 + " Goal utility: " + C0);
//
//			//Now set the heuristic val based on the analysis
//			if(B0 <= 0 && C0 > 0){
//				goalHeuristic = .9;
//				mobHeuristic = .1;
//			} else if (C0 <= 0 && B0 > 0){
//				mobHeuristic = .9;
//				goalHeuristic = .1;
//			} else if (!(C0 <= 0 && B0 <= 0)){
//				goalHeuristic += C0/(C0 + B0);
//				mobHeuristic += B0/(C0 + B0);
//			}
//
//			System.out.println("Mobility heuristic: " + mobHeuristic + " Goal heuristic: " + goalHeuristic);
//
//		}
	}

	/*
	 **********************METAGAME CODE**************************
	 *
	 */

	public List<List<Double>> runRandomGame(long timeout) throws MoveDefinitionException, GoalDefinitionException, TransitionDefinitionException{
		StateMachine machine = getStateMachine();
		List< List< Double > > mobilityRecord = new ArrayList<List<Double>>(machine.getRoles().size());
		List< List< Double > > goalRecord = new ArrayList<List<Double>>(machine.getRoles().size());
		List< Double > terminalValue = new ArrayList<Double>(machine.getRoles().size());

		/*
		 * We need to initialize the lists bc Java is dumb
		 */
		for(int i = 0; i < machine.getRoles().size(); i++){
			mobilityRecord.add(new ArrayList<Double>());
			goalRecord.add(new ArrayList<Double>());
			terminalValue.add((double)0);
		}

		MachineState curr = machine.findInits();
		while(true){
			if(machine.findTerminalp(curr)){
				//System.out.println("At terminal state");
				for(int i = 0; i < machine.getRoles().size(); i++){
					terminalValue.set(i, (double)machine.findReward(machine.getRoles().get(i), curr));
				}
				break;
			}
			//System.out.println("Playing another move");
			List<Move> currJointMove = machine.getRandomJointMove(curr);
			for(int i = 0; i < machine.getRoles().size(); i++){
				if(metaOOT(timeout)){
					return null;
				}

				Role r = machine.getRoles().get(i);

				//Record data
				double mob = mobility(r, curr);
				mobilityRecord.get(i).add(mob);

				int goal = goalUtility(r, curr);
				goalRecord.get(i).add((double)goal);
			}
			curr = machine.getNextState(curr, currJointMove);
		}

		//@David: now we have all the mobility and goal data points for each player stored,
		//in addition to final state - can you fill in the math part here?
		List< Double > netMobility = new ArrayList<Double>(machine.getRoles().size());
		List< Double > netGoal = new ArrayList<Double>(machine.getRoles().size());
		//System.out.println(mobilityRecord);
		//System.out.println(goalRecord);
		//System.out.println(terminalValue);

		for(int j = 0; j < machine.getRoles().size(); j++){
			netMobility.add(mean(mobilityRecord.get(j)));
			netGoal.add(mean(goalRecord.get(j)));
		}

		//System.out.println(netMobility);
		//System.out.println(netGoal);
		List<List<Double>> result = new ArrayList<List<Double>>();
		result.add(netMobility);
		result.add(netGoal);
		result.add(terminalValue);
		return result;
	}

	double mean(List<Double> vals){
		double sum = 0;
		for(int i = 0; i < vals.size(); i++){
			sum+=vals.get(i);
		}
		//System.out.println("Mean: " + sum/vals.size());
		return sum/vals.size();
	}

	double var(List<Double> vals, double mean){
		double sum = 0;
		for(int i = 0; i < vals.size(); i++){
			sum += (vals.get(i) - mean) * (vals.get(i) - mean);
		}
		//System.out.println("Sum in var: " + sum);
		//System.out.println("Var: " + sum/vals.size());
		return sum;
	}

	double covar(List<Double> vals, double mean, List<Double> y_vals, double y_mean){
		double sum = 0;
		//System.out.println("Mean: " + mean + " Y_Mean: " + y_mean);
		for(int i = 0; i < vals.size(); i++){
			sum += (vals.get(i) - mean)*(y_vals.get(i) - y_mean);
		}
		//System.out.println("Covar: " + sum);
		return sum;
	}

	/*
	 **********************END OF METAGAME CODE*******************
	 */

	@Override
	public Move stateMachineSelectMove(long timeout)
			throws TransitionDefinitionException, MoveDefinitionException,
			GoalDefinitionException {

		StateMachine machine = getStateMachine();
		MachineState state = getCurrentState();
		Role role = getRole();

		//This is bestMove without iterative deepening (so a fixed limit combined with variable heuristic)
		//return bestMove(role,state,machine);
//		if (machine.findLegals(role, state).size() == 1) return machine.findLegals(role, state).get(0);
//		ranOutOfTime = false;
//		return bestNetMoveID(role, state, machine, timeout);

		List<node> empty = new ArrayList<node>();
		node start = makeNode(false, state, 0, 0, null, empty, null, true);
		return bestMoveMCTS(timeout, start, role);

	}



	/*
	 *************THESE ARE FUNCTIONS WE ADDED*****************
	 */

	/*
	 * Stuff for MCTS (n-player)
	 *
	 */

	class node{
		int visited;
		List<node> children;
		double utility;
		node parent;
		MachineState state;
		Move move;
		Boolean allChildrenSeen;
		Boolean role; //Tells us if we're at a "min node" or a "max node" aka whether its our move.
	}

	private node makeNode(Boolean childrenSeen, MachineState state, int visited, double utility, node parent, List<node> children, Move move, Boolean ourRole){
		node result = new node();
		result.allChildrenSeen = childrenSeen;
		result.state = state;
		result.children= children;
		result.visited = visited;
		result.utility = utility;
		result.parent = parent;
		result.move = move;
		result.role = ourRole;
		return result;
	}

	private Move bestMoveMCTS(long timeout, node currNode, Role role) throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException{
		node workingNode;
		while(!outOfTime(timeout)){
			workingNode = select(currNode);
			expand(workingNode, role);
			double simVal = monteCarlo(role, workingNode.state, timeout);
			backProp(workingNode, simVal);
		}
		//Now we return a node
		double maxScore = 0;
		node bestNode = null;
		System.out.println(currNode.children.size());
		for(int i = 0; i < currNode.children.size(); i++){
			if(currNode.children.get(i).utility > maxScore){
				maxScore = currNode.children.get(i).utility;
				bestNode = currNode.children.get(i);
			}
		}
		return bestNode.move;
		//Now the best node is stored in bestNode. We need to find the move that gets us there.
	}

	//Select method for MCTS
	private node select(node currNode) throws MoveDefinitionException{
		if(currNode.visited == 0) return currNode;
		for(int i = 0; i < currNode.children.size(); i++){
			if(currNode.children.get(i).visited == 0) return currNode.children.get(i);
		}
		double score = 0;
		node result = currNode;
		for(int i = 0; i < currNode.children.size(); i++){
			double newscore = selectfn(currNode.children.get(i));
			if(newscore > score){
				score = newscore;
				result = currNode.children.get(i);
			}
		}
		return result;
	}

	//selectfn method for MCTS
	double selectfn(node curr){
		int coefficient = (curr.role) ? 1 : -1;
		return coefficient * curr.utility/curr.visited + Math.sqrt(2*Math.log(curr.parent.visited)/curr.visited);
	}

	//Expand function for MCTS
	private void expand(node currNode, Role role) throws TransitionDefinitionException, MoveDefinitionException, GoalDefinitionException{
		StateMachine machine = getStateMachine();
		if(machine.findTerminalp(currNode.state)){
			currNode.utility = machine.findReward(role, currNode.state);
			currNode.allChildrenSeen = true;
		} else if (currNode.role){
			List<Move> moves = machine.findLegals(role, currNode.state);
			for(int i = 0; i < moves.size(); i++){
				List<Move> nextMoves = new ArrayList<Move>();
				nextMoves.add(moves.get(i));
				MachineState newState = machine.getNextState(currNode.state, nextMoves);
				List<node> emptyChildren = new ArrayList<node>();
				node newNode = makeNode(false, newState, 0, 0, currNode, emptyChildren, moves.get(i), !currNode.role);
				currNode.children.add(newNode);
			}
		} else {
			List<List<Move>> moves = machine.getLegalJointMoves(currNode.state, role, currNode.move);
			for(int i = 0; i < moves.size(); i++){
				MachineState newState = machine.getNextState(currNode.state, moves.get(i));
				List<node> emptyChildren = new ArrayList<node>();
				node newNode = makeNode(false, newState, 0, 0, currNode, emptyChildren, null, !currNode.role);
				currNode.children.add(newNode);
			}
		}
	}

	//backpropogate function for MCTS
	private void backProp(node curr, double score){
		List<node> children = curr.children;
		Boolean allVisited = false;;
		for(int i = 0; i < children.size(); i++){
			allVisited = children.get(i).allChildrenSeen;
			if(!allVisited) break;
		}
		if(allVisited) curr.allChildrenSeen = true;

		//So here, basically if we recieve a score of 100 from a child, we KNOW
		//that we can win from here, so EMPHASIZE this node by adding more score/visits
		if(allVisited){
			curr.visited+=10;
			curr.utility += 10*score;
		} else {
			curr.visited++;
			curr.utility += score;
		}

		if(curr.parent != null)  backProp(curr.parent, score);
	}

	/*
	 * returns the best net move in iterative deepening (we assume to be the deepest returned move)
	 */

	private Move bestNetMoveID(Role role,MachineState state, StateMachine machine, long timeout) throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException{
		Stack<Move> selectedMoves = new Stack<Move>();
		int currDepth = 0; //We default one to search to at least a depth of one, right?
		while (currDepth <= LIMIT && !ranOutOfTime){ //We want 4s (4000ms) remaining
			if(!ranOutOfTime) selectedMoves.push(bestMoveID(role, state, machine, currDepth, timeout)); //The stack stores the current "best" move
			currDepth++;
			//timeDiff = System.currentTimeMillis() - temp;
			//System.out.println(timeout - System.currentTimeMillis());
		}
		return selectedMoves.pop();
	}

	/*
	 * Returns the best move, searching every move and seeing which one maximizes the result from minscore.
	 * Identical to best move from minimax
	 *
	 * USES ITERATIVE DEEPENING
	 */
	private Move bestMoveID(Role role, MachineState state, StateMachine machine, int depth, long timeout)
			throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException {
		List<Move> legals = machine.findLegals(role,state);
		Move action = legals.get(0);
		double score = 0;
		for (int i = 0; i < legals.size(); i++){
			double result = minscoreID(role, legals.get(i), state, 0, 0, depth, timeout); //Initially the mobility is 0, so we always recurse
			if (result > score) {
				score = result;
				action = legals.get(i);
			}
		}
		return action;
	}

	private Boolean outOfTime(long timeout){
		return ((timeout - System.currentTimeMillis()) < 2000);
	}

	private Boolean metaOOT(long timeout){
		return ((timeout - System.currentTimeMillis()) < 2000);
	}

	/*
	 * maxScore - essentially the max score fn from minimax; I removed alpha and beta to follow the code from the
	 * notes very strictly, although we could certainly try with those as well (shouldnt matter as much since
	 * we will never search the whole tree anyway).
	 *
	 * Also, note that a lot of the method names changed to directly match those in the notes (e.g. findReward)
	 *
	 * USES ITERATIVE DEEPENING
	 */
	private double maxScoreID(Role role, MachineState state, int level, double prevMobility, int depthLimit, long timeout)
			throws GoalDefinitionException, MoveDefinitionException, TransitionDefinitionException {
		StateMachine machine = getStateMachine();
		if (machine.findTerminalp(state)) return machine.findReward(role, state);
		//if ((level >= depthLimit)) return weightedHeuristic(role,state);
		if ((level >= depthLimit)) return monteCarlo(role,state,timeout);
		if(outOfTime(timeout)){
			ranOutOfTime = true;
			return weightedHeuristic(role, state);
		}
		//if (!varDepthMobility(role, state, prevMobility)) return goalUtility(role, state);
		//We need to get the curr mobility to pass to minscore (which basically saves it for the next maxscore)
		double currMobility = mobility(role, state);
		List<Move> actions = machine.findLegals(role, state);
		double score = 0;
		for (int i = 0; i<actions.size(); i++){
			double result = minscoreID(role, actions.get(i), state, level, currMobility, depthLimit, timeout);
			if (result == 100) return 100;
			if (result > score) score = result;
		}
		return score;
	}

	private double monteCarlo(Role role, MachineState state, long timeout) throws GoalDefinitionException, TransitionDefinitionException, MoveDefinitionException {
		StateMachine machine = getStateMachine();
		int total = 0;
		for (int i = 0; i < DEPTH_CHARGES; i++) {
			if (outOfTime(timeout)) break;
			int[] distance = new int[1];
			distance[0] = 0;
			total += machine.findReward(role, machine.performDepthCharge(state, distance));
		}
		if (outOfTime(timeout)) return 0;
		return (double) total / DEPTH_CHARGES;
	}

	/*
	 * minScore - essentially the min score fn from minimax; I removed alpha and beta to follow the code from the
	 * notes very strictly, although we could certainly try with those as well.
	 *
	 * USES ITERATIVE DEEPENING
	 */
	private double minscoreID(Role role, Move action, MachineState state, int level, double currMobility, int depthLimit, long timeout)
			throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException {
		StateMachine machine = getStateMachine();
		List<List<Move>> legals = machine.getLegalJointMoves(state, role, action);
		double score = 100;
		for (int i = 0; i<legals.size(); i++) {
			List<Move> move = legals.get(i);
			MachineState newState = machine.getNextState(state, move);
			double result = maxScoreID(role, newState, level + 1, currMobility, depthLimit, timeout);
			if(ranOutOfTime) return score;
			if (result == 0) return 0;
			if(result < score) score = result;
		}
		return score;
	}

	/*
	 * Returns the best move, searching every move and seeing which one maximizes the result from minscore.
	 * Identical to best move from minimax
	 */
//	private Move bestMove(Role role, MachineState state, StateMachine machine)
//			throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException {
//		List<Move> legals = machine.findLegals(role,state);
//		Move action = legals.get(0);
//		double score = 0;
//		for (int i = 0; i < legals.size(); i++){
//			double result = minscore(role, legals.get(i), state, 0, 0); //Initially the mobility is 0, so we always recurse
//			if (result > score) {
//				score = result;
//				action = legals.get(i);
//			}
//		}
//		return action;
//	}

	/*
	 * maxScore - essentially the max score fn from minimax; I removed alpha and beta to follow the code from the
	 * notes very strictly, although we could certainly try with those as well (shouldnt matter as much since
	 * we will never search the whole tree anyway).
	 *
	 * Also, note that a lot of the method names changed to directly match those in the notes (e.g. findReward)
	 */
//	private double maxScore(Role role, MachineState state, int level, double prevMobility)
//			throws GoalDefinitionException, MoveDefinitionException, TransitionDefinitionException {
//		StateMachine machine = getStateMachine();
//		if (machine.findTerminalp(state)) return machine.findReward(role, state);
//		if (level >= LIMIT || !varDepthMobility(role, state, prevMobility)) return weightedHeuristic(role,state);
//		//if (!varDepthMobility(role, state, prevMobility)) return goalUtility(role, state);
//		//We need to get the curr mobility to pass to minscore (which basically saves it for the next maxscore)
//		double currMobility = mobility(role, state);
//		List<Move> actions = machine.findLegals(role, state);
//		double score = 0;
//		for (int i = 0; i<actions.size(); i++){
//			double result = minscore(role, actions.get(i), state, level, currMobility);
//			if (result == 100) return 100;
//			if (result > score) score = result;
//		}
//		return score;
//	}

	/*
	 * minScore - essentially the min score fn from minimax; I removed alpha and beta to follow the code from the
	 * notes very strictly, although we could certainly try with those as well.
	 */
//	private double minscore(Role role, Move action, MachineState state, int level, double currMobility)
//			throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException {
//		StateMachine machine = getStateMachine();
//		//if (machine.isTerminal(state)) return machine.getGoal(state, role);
//		List<List<Move>> legals = machine.getLegalJointMoves(state, role, action);
//		double score = 100;
//		for (int i = 0; i<legals.size(); i++) {
//			List<Move> move = legals.get(i);
//			MachineState newState = machine.getNextState(state, move);
//			double result = maxScore(role, newState, level + 1, currMobility);
//			if (result == 0) return 0;
//			if(result < score) score = result;
//		}
//		return score;
//	}

	/*
	 * Variable depth heuristic: NOT IN USE & UNTESTED
	 *
	 * I tried to capture the gradient of mobility
	 *
	 * Returns true if we should continue searching, false if we should not continue searching.
	 * We continue if the mobility in the next state changes by a factor of >1.5x current mobility
	 */

	private Boolean varDepthMobility(Role role, MachineState state, double prevMobility) throws MoveDefinitionException {
		double currMobility = mobility(role, state);
		if(currMobility == 0|| prevMobility == 0) return true;
		double ratio = (double)prevMobility / (double)currMobility;
		return (ratio >= 1.5 || 1/ratio >= 1.5); //Returns true if there is a big change.
	}

	/*
	 * mobility - Returns a percentage equal to the ratio of possible actions in the current state to total
	 * available actions
	 *
	 * Maybe we should return a double?
	 */
	private double mobility(Role role, MachineState state) throws MoveDefinitionException {
		StateMachine machine = getStateMachine();
		List<Move> actions = machine.findLegals(role, state);
		List<Move> feasible = machine.findActions(role);
		return 100*((double)actions.size()/(double)feasible.size());
	}

	/*
	 * focus - Returns a percentage equal to 100 - ratio of possible actions in the current state to total
	 * available actions
	 *
	 * Maybe we should return a double?
	 */
	private int focus(Role role, MachineState state) throws MoveDefinitionException {
		StateMachine machine = getStateMachine();
		List<Move> actions = machine.findLegals(role, state);
		List<Move> feasible = machine.findActions(role);
		return 100 - (100*(actions.size()/feasible.size()));
	}

	/*
	 * goal proximity - Returns the utility of the current state (should be increasing as we get
	 * close to goal assuming the game is monotonic)
	 *
	 * Maybe we should return a double?
	 *
	 * This one almost nailed hunter. We had one pawn left
	 */
	private int goalUtility(Role role, MachineState state) throws MoveDefinitionException, GoalDefinitionException {
		StateMachine machine = getStateMachine();
		return machine.getGoal(state, role);
	}

	/*
	 * Weighted heuristic fn which combines focus, mobility, and goal utility. We need weights
	 */
	private double weightedHeuristic(Role role, MachineState state) throws MoveDefinitionException, GoalDefinitionException{
		return goalHeuristic * goalUtility(role, state) + mobHeuristic * mobility(role, state);
	}

	/*
	 *************END OF FUNCTIONS WE ADDED********************
	 */

	@Override
	public void stateMachineStop() {
		// TODO Auto-generated method stub

	}

	@Override
	public void stateMachineAbort() {
		// TODO Auto-generated method stub

	}

	@Override
	public void preview(Game g, long timeout) throws GamePreviewException {
		// TODO Auto-generated method stub

	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return "Mike's Angels";
	}

}
