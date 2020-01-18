/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package quasylab.fabs.companion;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class App {
	
    public static void main(String[] args) throws FileNotFoundException, InterruptedException {

    	boolean[] flags = getSimulationFlag(args);
    	int[] scales = scale(args);
    	System.out.println(Arrays.toString(scales));
    	if (flags[0]) {
    		simulateRBModel( scales );
    	}
    	if (flags[1]) {
    		simulateGossipUnicastModel( scales );
    	}
    	if (flags[2]) {
    		simulateGossipBroadcastModel( scales );
    	}
    	System.exit(0);
    }
    
    private static void simulateGossipBroadcastModel(int[] scales) throws FileNotFoundException, InterruptedException {
		System.out.println("SIMULATION OF GOSSIP-BROADCAST MODEL STARTED!");    	
		GossipBroadcast model = new GossipBroadcast(GossipBroadcast.K,GossipBroadcast.DIFF_RATE,GossipBroadcast.PASS_RATE);
		for( int i=0 ; i<scales.length ; i++ ) {
			System.out.println("SIMULATING SCALE "+scales[i]);
			model.run(scales[i], 1, GossipBroadcast.DEADLINE, GossipBroadcast.SAMPLINGS, "./data/");
		}
		System.out.println("GOSSIP-BROADCAST MODEL DONE!\n\n");    	
	}

	private static void simulateGossipUnicastModel(int[] scales) throws FileNotFoundException, InterruptedException {
		System.out.println("SIMULATION OF GOSSIP-UNICAST MODEL STARTED!");    	
		GossipUnicast model = new GossipUnicast(GossipUnicast.DIFF_RATE,GossipUnicast.PASS_RATE);
		for( int i=0 ; i<scales.length ; i++ ) {
			System.out.println("SIMULATING SCALE "+scales[i]);
			model.run(scales[i], 1, GossipUnicast.DEADLINE, GossipUnicast.SAMPLINGS, "./data/");
		}
		System.out.println("GOSSIP-BROADCAST MODEL DONE!\n\n");    	
	}

	private static void simulateRBModel(int[] scales) throws FileNotFoundException, InterruptedException {
		System.out.println("SIMULATION OF RED-BLUE MODEL STARTED!");    	
    	RBModel model = new RBModel(10,1,1,0.5,0.5);	
		for( int i=0 ; i<scales.length ; i++ ) {
			System.out.println("SIMULATING SCALE "+scales[i]);
			model.run(scales[i], 1, RBModel.DEADLINE, RBModel.SAMPLINGS, "./data/");			
		}
		System.out.println("RED-BLUE MODEL DONE!\n\n");    	
	}

	public static int[] scale(String[] args) {
		return new int[] { 1, 10, 100, 1000 };
    }

	private static boolean[] getSimulationFlag(String[] args) {
		if ((args.length==0)||("all".equals(args[0]))) {
			return new boolean[] { true, true, true };
		}
		if ("rb".equals(args[0])) {
			return new boolean[] { true, false, false };
		}
		if ("ugossip".equals(args[0])) {
			return new boolean[] { false, true, false };
		}
		if ("bgossip".equals(args[0])) {
			return new boolean[] { false, false, true };
		}
		throw new IllegalArgumentException("Unexpected parameter "+args[0]+"!");
	}
}
