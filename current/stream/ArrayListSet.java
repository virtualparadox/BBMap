package stream;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;

/**
 * @author Brian Bushnell
 * @date Apr 2, 2015
 *
 */
public class ArrayListSet {

	public ArrayListSet(boolean ordered_){
		ordered=ordered_;
		assert(!ordered) : "Ordered output does not currently work with this program.";
	}
	
	public ArrayListSet(boolean ordered_, Iterable<String> names_){
		ordered=ordered_;
		for(String s : names_){Pack p=new Pack(s);}
	}
	
	public void add(Read r, Iterable<String> names_){
		for(String s : names_){add(r, s);}
	}
	
	public void add(Read r, String name){
		final Pack p=getPack(name, true);
		p.add(r);
	}
	
	public void add(Read r, int id){
		final Pack p=getPack(id, true);
		p.add(r);
	}
	
	public ArrayList<Read> getAndClear(String name){
		final Pack p=getPack(name, false);
		return p==null ? null : p.getAndClear();
	}
	
	public ArrayList<Read> getAndClear(int id){
		final Pack p=getPack(id, false);
		return p==null ? null : p.getAndClear();
	}
	
	public Collection<String> getNames(){
		return nameList;
	}
	
	private Pack getPack(String name, boolean add){
		Pack p=stringMap.get(name);
		if(p==null && add){p=new Pack(name);}
		return p;
	}
	
	private Pack getPack(int id, boolean add){
		Pack p=packList.size()>id ? packList.get(id) : null;
		if(p==null && add){p=new Pack(id);}
		return p;
	}
	
	private class Pack {
		
		Pack(String s){
			assert(s==null || !stringMap.containsKey(s));
			name=s;
			id=packList.size();
			nameList.add(s);
			packList.add(this);
			if(s!=null){stringMap.put(s, this);}
		}
		
		Pack(int x){
			name=null;
			id=x;
			while(packList.size()<=x){packList.add(null);}
			assert(packList.get(x)==null);
			packList.set(x, this);
		}
		
		public void add(Read r){
			if(list==null){list=new ArrayList<Read>();}
			list.add(r);
		}
		
		public ArrayList<Read> getAndClear(){
			ArrayList<Read> temp=list;
			list=null;
			return temp;
		}
		
		final String name;
		final int id;
		private ArrayList<Read> list;
	}
	
	private final boolean ordered;
	private final ArrayList<String> nameList=new ArrayList<String>();
	private final ArrayList<Pack> packList=new ArrayList<Pack>();
	private final LinkedHashMap<String, Pack> stringMap=new LinkedHashMap<String, Pack>();
	
}
