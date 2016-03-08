package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import align2.ListNum;
import align2.Shared;

/**
 * This class is designed for distributed environments.
 * The 'master' reads from the filesystem, creates reads, and broadcasts them.
 * The 'slaves' listen for broadcasts.
 * @author Brian Bushnell
 * @date Oct 7, 2014
 *
 */
public class ConcurrentReadInputStreamD implements ConcurrentReadStreamInterface {
	
	public ConcurrentReadInputStreamD(ConcurrentReadStreamInterface cris_, boolean master_){
		source=cris_;
		master=master_;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		assert(master==(cris_!=null));
		
		if(master){
			paired=source.paired();
			broadcastPaired(paired);
		}else{
			paired=listenPaired();
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("**************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("**************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try {
				list=depot.full.take();
				assert(list!=null);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("**************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	public void returnList(ListNum<Read> ln, boolean poison){
		if(ln!=null){
			ln.list.clear();
		}else{
			System.err.println("Warning, null list returned: ");  //System.err.println("Warning from class "+getClass().getName()+", null list returned: ");
			new Exception().printStackTrace();
		}
		if(poison){
			if(verbose){System.err.println("A: Adding empty list to full.");}
			depot.full.add(ln==null ? new ArrayList<Read>(0) : ln.list);
		}else{
			if(ln!=null){depot.empty.add(ln.list);}
//			depot.empty.add(ln==null ? new ArrayList<Read>(0) : ln.list);
		}
	}
	
	public void returnList(boolean poison){
		if(poison){
			depot.full.add(new ArrayList<Read>(0));
		}else{
			depot.empty.add(new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
		}
	}
	
	@Override
	public void run() {
		synchronized(running){
			assert(!running[0]) : "This cris was started by multiple threads.";
			running[0]=true;
		}
		if(verbose){System.err.println("cris started.");}
		threads=new Thread[] {Thread.currentThread()};
		
		if(master){
			readLists_master();
		}else{
			readLists_slave();
		}

		addPoison();
		
		//End thread

		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("Ending");
			if(verbose){System.err.println("B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
		if(verbose){System.err.println("cris thread syncing before shutdown.");}
		
		synchronized(running){//TODO Note: for some reason syncing on 'this' instead of 'running' causes a hang.  Something else must be syncing improperly on this.
			assert(running[0]);
			running[0]=false;
		}
		if(verbose){System.err.println("cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());}
	}
	
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
//					System.err.println("Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("D: Adding list("+list.size()+") to full.");}
				depot.full.add(list);
			}
		}
		if(verbose){System.err.println("Added poison.");}
	}
	
	private final void readLists_master(){

		if(verbose){System.err.println("Entered readLists_master().");}
		for(ListNum<Read> ln=source.nextList(); !shutdown && ln.list!=null; ln=source.nextList()){
			final ArrayList<Read> reads=ln.list;
			final int count=(reads==null ? 0 : reads.size());
			
			if(verbose){System.err.println("Master fetched "+count+" reads.");}
			
			try {
				depot.full.put(reads);
				if(verbose){System.err.println("Master added reads to depot.");}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			broadcast(reads);
			if(verbose){System.err.println("Master broadcasted.");}
			source.returnList(ln, count<1);
			if(verbose){System.err.println("Master returned a list.");}
			if(count<1){break;}
		}
		if(verbose){System.err.println("Finished readLists_master().");}
	}
	
	private final void readLists_slave(){
		
		if(verbose){System.err.println("Entered readLists_slave().");}
		for(ArrayList<Read> reads=listen(); !shutdown && reads!=null; reads=listen()){
			final int count=(reads==null ? 0 : reads.size());
			
			if(verbose){System.err.println("Slave fetched "+count+" reads.");}

			try {
				depot.full.put(reads);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(count<1){break;}
		}
		if(verbose){System.err.println("Finished readLists_slave().");}
	}
	
	/*--------------------------------------------------------------*/
	
	private void broadcast(ArrayList<Read> ln){
		if(verbose){System.err.println("Broadcasting reads.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	private void broadcastPaired(boolean b){
		if(verbose){System.err.println("Broadcasting pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	private ArrayList<Read> listen(){
		if(verbose){System.err.println("Listening for reads.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	private boolean listenPaired(){
		if(verbose){System.err.println("Listening for pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public void shutdown(){
		if(verbose){System.out.println("Called shutdown.");}
		
		shutdown=true;
		if(!shutdown){
			
			if(master){
				source.shutdown();
			}
			for(Thread t : threads){
				if(t!=null && t.isAlive()){
					t.interrupt();
				}
			}
		}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		basesIn=0;
		readsIn=0;
		listnum=0; //Added Oct 9, 2014
		if(master){
			source.restart();
		}
	}
	
	@Override
	public synchronized void close(){
		shutdown();
		
		if(master){
			source.close();
		}else{
			
		}
		
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			
			while(threads[0].isAlive()){
//				System.out.println("B");
				ArrayList<Read> list=null;
				for(int i=0; i<1000 && list==null && threads[0].isAlive(); i++){
					try {
						list=depot.full.poll(200, TimeUnit.MILLISECONDS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						System.err.println("Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
				
//				System.out.println("isAlive? "+threads[0].isAlive());
			}
			
		}
		
		if(threads!=null){
			for(int i=1; i<threads.length; i++){
				while(threads[i]!=null && threads[i].isAlive()){
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/

	@Override
	public boolean paired() {
		return paired;
	}
	
	@Override
	public void setSampleRate(float rate, long seed){
		if(master){source.setSampleRate(rate, seed);}
	}
	
	public long basesIn(){return basesIn;}
	public long readsIn(){return readsIn;}
	
	@Override
	public boolean errorState(){
		if(master){return errorState|source.errorState();}
		return errorState;
	}
	
	public Object[] producers(){
		if(master){return source.producers();}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadStreamInterface source;
	private final boolean master;
	
	private boolean errorState=false;
	
	private boolean[] running=new boolean[] {false};
	
	private boolean shutdown=false;

	private ConcurrentDepot<Read> depot;

	private Thread[] threads;
	
	private long basesIn=0;
	private long readsIn=0;
	
	private long listnum=0;
	
	/** This should be set in the first broadcast */
	private final boolean paired;

	private final int BUF_LEN=Shared.READ_BUFFER_LENGTH;
	private final int NUM_BUFFS=Shared.READ_BUFFER_NUM_BUFFERS;
	
	public static boolean verbose=false;
	
}
