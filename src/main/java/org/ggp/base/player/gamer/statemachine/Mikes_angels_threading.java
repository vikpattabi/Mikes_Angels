package org.ggp.base.player.gamer.statemachine;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.ggp.base.apps.player.Player;
import org.ggp.base.player.gamer.exception.GamePreviewException;
import org.ggp.base.util.game.Game;
import org.ggp.base.util.statemachine.MachineState;
import org.ggp.base.util.statemachine.Move;
import org.ggp.base.util.statemachine.Role;
import org.ggp.base.util.statemachine.StateMachine;
import org.ggp.base.util.statemachine.exceptions.GoalDefinitionException;
import org.ggp.base.util.statemachine.exceptions.MoveDefinitionException;
import org.ggp.base.util.statemachine.exceptions.TransitionDefinitionException;
import org.ggp.base.util.statemachine.implementation.propnet.SamplePropNetStateMachine;


public class Mikes_angels_threading extends StateMachineGamer {

	//Limit for depth of our recursive search.
	private static final int DEPTH_CHARGES = 6;
	private int numCharges;
	Player p;
	private ArrayList<StateMachine> threadMachines;

	private boolean greedyMode;
	private ExecutorService pool;

	@Override
	public StateMachine getInitialStateMachine() {
		// TODO Auto-generated method stub
		SamplePropNetStateMachine machine = new SamplePropNetStateMachine();
		machine.initialize(getMatch().getGame().getRules());
		return machine;
	}

	public void performDepthChargeTimed(MachineState state, long ref, StateMachine machine) throws TransitionDefinitionException, MoveDefinitionException {
		while(!machine.isTerminal(state)) {
			state = machine.getNextStateDestructively(state, machine.getRandomJointMove(state));
			if (System.currentTimeMillis() - ref > 1500) {
				greedyMode = true;
				System.out.println("Let's get greedy.");
			}
		}
		greedyMode = false;
		System.out.println("No greed today.");
	}

	@Override
	public void stateMachineMetaGame(long timeout)
			throws TransitionDefinitionException, MoveDefinitionException,
			GoalDefinitionException {

		System.out.println("hi meta");
		StateMachine machine = getStateMachine();
		long ref = System.currentTimeMillis();
		pool = Executors.newFixedThreadPool(DEPTH_CHARGES);
		performDepthChargeTimed(getCurrentState(),ref,machine);
		if(greedyMode){
			System.out.println("Not initializing threadsafe machines. Done with metagame.");
			return;
		}
		threadMachines = new ArrayList<StateMachine>();
		for(int i = 0; i < DEPTH_CHARGES; i++){
			StateMachine curr = new SamplePropNetStateMachine();
			curr.initialize(getMatch().getGame().getRules());
			threadMachines.add(curr);
		}
		System.out.println("Done with metagame.");
	}

	/*
	 **********************METAGAME CODE**************************
	 *
	 */

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
			GoalDefinitionException, InterruptedException, ExecutionException {

		MachineState state = getCurrentState();
		Role role = getRole();
		StateMachine machine = getStateMachine();

		if (machine.findLegals(role, state).size() == 1) return machine.findLegals(role, state).get(0);

		if(greedyMode) {
			Random rdn = new Random();
			int choice = rdn.nextInt(machine.getLegalMoves(state, role).size());
			return machine.getLegalMoves(state, role).get(choice);
		}

		List<node> emptyChildren = new ArrayList<node>();
		node treeRoot = makeNode(state, 0, 0, null, emptyChildren, null);

		if(machine.getRoles().size() == 1){
			return bestMoveMCTS(treeRoot, timeout, role);
		}else {
			return bestMoveMCTS_Multiplayer(treeRoot, timeout, role);
		}
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
	}

	private node makeNode(MachineState state, int visited, double utility, node parent, List<node> children, Move move){
		node result = new node();
		result.state = state;
		result.children= children;
		result.visited = visited;
		result.utility = utility;
		result.parent = parent;
		result.move = move;
		return result;
	}

	private Move bestMoveMCTS(node treeRoot, long timeout, Role role) throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException, InterruptedException, ExecutionException{
		StateMachine machine = getStateMachine();
		node workingNode;

		while(!outOfTime(timeout)){
			workingNode = select(treeRoot);
			expandOnePlayer(workingNode, role);
			backProp(workingNode, monteCarlo(role, workingNode.state, timeout));
		}
		System.out.println("Out of time for MCTS loop.");

		//Now we return a node
		double maxScore = -1.0;
		node bestNode = null;
		System.out.println("Finished current round of MCTS.");
		System.out.println(treeRoot.children.size());
		System.out.println(machine.getLegalMoves(treeRoot.state, role).size());
		System.out.println("");
		for(int i = 0; i < treeRoot.children.size(); i++){
			if(treeRoot.children.get(i).utility > maxScore){
				maxScore = treeRoot.children.get(i).utility;
				bestNode = treeRoot.children.get(i);
			}
		}
		treeRoot = bestNode;
		return bestNode.move;
	}

	private Move bestMoveMCTS_Multiplayer(node treeRoot, long timeout, Role role) throws MoveDefinitionException, TransitionDefinitionException, GoalDefinitionException, InterruptedException, ExecutionException{
		numCharges = 0;
		StateMachine machine = getStateMachine();
		node workingNode;

		while(!outOfTime(timeout)){
			workingNode = selectMultiplayer(treeRoot);
			if(machine.isTerminal(workingNode.state)) continue;
			expandMultiplayer(workingNode, role);
			backPropMultiplayer(workingNode, monteCarlo(role, workingNode.state, timeout));
		}

		//Now we return a node
		double maxScore = -1.0;
		node bestNode = treeRoot.children.get(0);
		System.out.println("Finished current round of MCTS.");
		System.out.println(treeRoot.children.size());
		System.out.println(machine.getLegalMoves(treeRoot.state, role).size());
		System.out.println("");
		for(int i = 0; i < treeRoot.children.size(); i++){
			if(treeRoot.children.get(i).utility > maxScore){
				maxScore = treeRoot.children.get(i).utility;
				bestNode = treeRoot.children.get(i);
			}
		}


		System.out.println("Depth charges: " + numCharges);
		return bestNode.move;
	}

	//Select method for MCTS
	private node select(node currNode) throws MoveDefinitionException{
		StateMachine machine = getStateMachine();
		if(currNode.visited == 0 || machine.isTerminal(currNode.state)) return currNode;

		if(currNode.children.size() == 0) return currNode;

		for(int i = 0; i < currNode.children.size(); i++){
			if(currNode.children.get(i).visited == 0 && !machine.isTerminal(currNode.children.get(i).state)) return currNode.children.get(i);
		}

		double score = -1.0;
		node result = currNode.children.get(0);
		for(int i = 0; i < currNode.children.size(); i++){
			double newscore = selectfn(currNode.children.get(i));
			if(newscore > score && !machine.isTerminal(currNode.children.get(i).state)){
				score = newscore;
				result = currNode.children.get(i);
			}
		}
		return select(result);
	}

	//Select method for MCTS
	private node selectMultiplayer(node currNode) throws MoveDefinitionException{
		StateMachine machine = getStateMachine();
		if(currNode.visited == 0 || machine.isTerminal(currNode.state)) return currNode;

		for(int i = 0; i < currNode.children.size(); i++){
			List<node> grandchildren = currNode.children.get(i).children;
			for(int j = 0; j < grandchildren.size(); j++){
				if(grandchildren.get(j).visited == 0 || machine.isTerminal(grandchildren.get(j).state)) return grandchildren.get(j);
			}
		}

		double score = -1.0;
		node result = currNode.children.get(0).children.get(0);
		for(int i = 0; i < currNode.children.size(); i++){
			List<node> grandchildren = currNode.children.get(i).children;
			for(int j = 0; j < grandchildren.size(); j++){
				double newscore = selectfnmp(grandchildren.get(j));
				//System.out.println("newscore: " + newscore);
				if(newscore > score){
					score = newscore;
					result = grandchildren.get(j);
				}
			}
		}

		return selectMultiplayer(result);
	}

	//selectfn method for MCTS
	double selectfn(node curr){
		 return curr.utility/curr.visited+Math.sqrt(2*Math.log(curr.parent.visited)/curr.visited);
	}

	double selectfnmp(node curr){
		 return curr.utility/curr.visited+Math.sqrt(2*Math.log(curr.parent.parent.visited)/curr.visited);
	}

	//Expand function for MCTS
	private void expandOnePlayer(node currNode, Role role) throws TransitionDefinitionException, MoveDefinitionException, GoalDefinitionException{
		StateMachine machine = getStateMachine();
		if(machine.isTerminal(currNode.state)) return;
		List<Move> moves = machine.findLegals(role, currNode.state);
		for(int i = 0; i < moves.size(); i++){
			List<Move> nextMoves = new ArrayList<Move>();
			nextMoves.add(moves.get(i));
			MachineState newState = machine.getNextState(currNode.state, nextMoves);
			List<node> emptyChildren = new ArrayList<node>();
			node newNode = makeNode(newState, 0, 0, currNode, emptyChildren, moves.get(i));
			currNode.children.add(newNode);
		}
	}

	//Expand function for MCTS
	private void expandMultiplayer(node currNode, Role role) throws TransitionDefinitionException, MoveDefinitionException, GoalDefinitionException{
		StateMachine machine = getStateMachine();
		if(machine.isTerminal(currNode.state)) return;
		List<Move> moves = machine.findLegals(role, currNode.state);

		for(int i = 0; i < moves.size(); i++){
			List<node> emptyChildren = new ArrayList<node>();
			node newNode = makeNode(null, 0, 0, currNode, emptyChildren, moves.get(i));
			currNode.children.add(newNode);

			List<List<Move>> jointMoves = machine.getLegalJointMoves(currNode.state, role, moves.get(i));
			for(int j = 0; j < jointMoves.size(); j++){
				MachineState next = machine.getNextState(currNode.state, jointMoves.get(j));
				List<node> emptyChildrenNew = new ArrayList<node>();
				node newStateNode = makeNode(next, 0, 0, newNode, emptyChildrenNew, moves.get(i));
				newNode.children.add(newStateNode);
			}
		}
	}

	//backpropogate function for MCTS
	private void backProp(node curr, double score){
		curr.visited++;
		curr.utility += score;
		if(curr.parent != null)  backProp(curr.parent, score);
	}

	private void backPropMultiplayer(node curr, double score){
		curr.visited++;
		curr.utility += score;
		if(curr.parent != null){
			backPropMultiplayer(curr.parent, score);
		}
	}

	private Boolean outOfTime(long timeout){
		return ((timeout - System.currentTimeMillis()) < 3000);
	}

	private double monteCarlo(Role role, MachineState state, long timeout) throws GoalDefinitionException, TransitionDefinitionException, MoveDefinitionException, InterruptedException, ExecutionException {
		StateMachine machine = getStateMachine();
		if(machine.isTerminal(state)) return machine.getGoal(state, role);

		ArrayList<Future<Integer>> futures = new ArrayList<Future<Integer>>();
		int total = 0;
		for (int i = 0; i < DEPTH_CHARGES; i++) {
			if (outOfTime(timeout)) break;
			DepthCharge charge = new DepthCharge(state, threadMachines.get(i));
			@SuppressWarnings("unchecked")
			Future<Integer> result = pool.submit(charge);
			futures.add(result);
			numCharges++;
		}
		if (outOfTime(timeout)) return goalUtility(role, state);
		for(int i = 0; i < futures.size(); i++){
			total += futures.get(i).get();
		}
		return (double) total / futures.size();
	}

	/*
	 * --------------------CODE FOR THREADING-----------------------
	 */

	public class DepthCharge implements Callable{
		private StateMachine thread_machine;
		private MachineState ourState;


		public DepthCharge(MachineState state, StateMachine in){
			ourState = state;
			thread_machine = in;
		}

		@Override
		public Integer call() throws GoalDefinitionException, TransitionDefinitionException, MoveDefinitionException{
			int[] distance = new int[1];
			distance[0] = 0;
			int result = thread_machine.findReward(getRole(), thread_machine.performDepthCharge(ourState, distance));
			return result;
		}
	}


	/*
	 * -------------------------------------------------------------
	 */

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
		return "Mike's Angels w/ threading";
	}

}
