package MPM_CudaSim;

import java.io.*;
import java.util.*;

//this class will manage file io
public class FileIOManager{
	//owning map manager
	protected MPM_Abs_CUDASim expMgr;
	//name of owning class of the instance of this object, for display
	protected String owner;
	
	public FileIOManager(MPM_Abs_CUDASim _expMgr, String _owner) {expMgr = _expMgr; owner=_owner;}	
	
	private String buildClrStr(ConsoleCLR bk, ConsoleCLR clr, String str) {return bk.toString() + clr.toString() + str + ConsoleCLR.RESET.toString();	}
	private String _processMsgCode(String src, MsgCodes useCode) {
		if (!MPM_Abs_CUDASim.supportsANSITerm) {return src;}
		switch(useCode) {//add background + letter color for messages
			//info messages
			case info1 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.WHITE, src);}		//basic informational printout
			case info2 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.CYAN, src);}
			case info3 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.YELLOW, src);}		//informational output from som EXE
			case info4 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.GREEN, src);}
			case info5 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.CYAN_BOLD, src);}	//beginning or ending of processing chuck/function
			//warning messages                                                 , 
			case warning1 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLACK_BOLD, src);}
			case warning2 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_BOLD, src);}	//warning info re: ui does not exist
			case warning3 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLACK_UNDERLINED, src);}
			case warning4 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_UNDERLINED, src);}	//info message about unexpected behavior
			case warning5 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_BRIGHT, src);}
			//error messages                                                   , 
			case error1 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_UNDERLINED, src);}//try/catch error
			case error2 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_BOLD, src);}		//code-based error
			case error3 : {		return  buildClrStr(ConsoleCLR.RED_BACKGROUND_BRIGHT, ConsoleCLR.BLACK_BOLD, src);}	//file load error
			case error4 : {		return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND_BRIGHT, ConsoleCLR.RED_BRIGHT, src);}	//error message thrown by som executable
			case error5 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_BOLD_BRIGHT, src);}
		}
		return src;
	}//_processMsgCode	
	private void dispMessage(String srcClass, String srcMethod, String msgText, MsgCodes useCode) {_dispMessage_base(expMgr.getTimeStrFromProcStart() +"|" + srcClass,srcMethod,msgText, useCode,true);	}
	private void dispMessage(String srcClass, String srcMethod, String msgText, MsgCodes useCode, boolean onlyConsole) {_dispMessage_base(expMgr.getTimeStrFromProcStart() +"|" + srcClass,srcMethod,msgText, useCode,onlyConsole);	}	
	private void _dispMessage_base(String srcClass, String srcMethod, String msgText, MsgCodes useCode, boolean onlyConsole) {		
		String msg = _processMsgCode(srcClass + "::" + srcMethod + " : " + msgText, useCode);
		if((onlyConsole) || (expMgr.pa == null)) {		System.out.println(msg);	} else {		expMgr.pa.outStr2Scr(msg);	}
	}//dispMessage
	
	//write data to file
	public void saveStrings(String fname, String[] data) {
		PrintWriter pw = null;
		try {
		     File file = new File(fname);
		     FileWriter fw = new FileWriter(file, false);
		     pw = new PrintWriter(fw);
		     for (int i=0;i<data.length;++i) { pw.println(data[i]);}
		     
		} catch (IOException e) {	e.printStackTrace();}
		finally {			if (pw != null) {pw.close();}}
	}//saveStrings

	public void saveStrings(String fname, ArrayList<String> data) {
		PrintWriter pw = null;
		try {
		     File file = new File(fname);
		     FileWriter fw = new FileWriter(file, false);
		     pw = new PrintWriter(fw);
		     for (int i=0;i<data.size();++i) { pw.println(data.get(i));}
		     
		} catch (IOException e) {	e.printStackTrace();}
		finally {			if (pw != null) {pw.close();}}
	}//saveStrings
	
	public String[] loadFileIntoStringAra(String fileName, String dispYesStr, String dispNoStr) {try {return _loadFileIntoStringAra(fileName, dispYesStr, dispNoStr);} catch (Exception e) {e.printStackTrace(); } return new String[0];}
	//stream read the csv file and build the data objects
	private String[] _loadFileIntoStringAra(String fileName, String dispYesStr, String dispNoStr) throws IOException {		
		FileInputStream inputStream = null;
		Scanner sc = null;
		List<String> lines = new ArrayList<String>();
		String[] res = null;
	    //int line = 1, badEntries = 0;
		try {
		    inputStream = new FileInputStream(fileName);
		    sc = new Scanner(inputStream);
		    while (sc.hasNextLine()) {lines.add(sc.nextLine()); }
		    //Scanner suppresses exceptions
		    if (sc.ioException() != null) { throw sc.ioException(); }
		    dispMessage("fileIOManager:"+owner, "_loadFileIntoStringAra",dispYesStr+"\tLength : " +  lines.size(), MsgCodes.info3);
		    res = lines.toArray(new String[0]);		    
		} catch (Exception e) {	
			e.printStackTrace();
			expMgr.dispMessage("fileIOManager:"+owner, "_loadFileIntoStringAra","!!"+dispNoStr, MsgCodes.error3);
			res= new String[0];
		} 
		finally {
		    if (inputStream != null) {inputStream.close();		    }
		    if (sc != null) { sc.close();		    }
		}
		return res;
	}//loadFileContents	
	
	//load into multiple arrays for multi-threaded processing
	public String[][] loadFileIntoStringAra_MT(String fileName, String dispYesStr, String dispNoStr, int numHdrLines, int numThds) {
		try {return _loadFileIntoStringAra_MT(fileName, dispYesStr, dispNoStr, numHdrLines, numThds);} 
		catch (Exception e) {e.printStackTrace(); } 
		return new String[0][];
	}
	//load files into multiple arrays for multi-threaded processing
	private String[][] _loadFileIntoStringAra_MT(String fileName, String dispYesStr, String dispNoStr, int numHdrLines, int numThds) throws IOException {		
		FileInputStream inputStream = null;
		Scanner sc = null;
		List<String>[] lines = new ArrayList[numThds];
		for (int i=0;i<numThds;++i) {lines[i]=new ArrayList<String>();	}
		String[][] res = new String[numThds+1][];
		String[] hdrRes = new String[numHdrLines];
		int idx = 0, count = 0;
		try {
		    inputStream = new FileInputStream(fileName);
		    sc = new Scanner(inputStream);
		    for(int i=0;i<numHdrLines;++i) {    	hdrRes[i]=sc.nextLine();   }		    
		    while (sc.hasNextLine()) {
		    	lines[idx].add(sc.nextLine()); 
		    	idx = (idx + 1)%numThds;
		    	++count;
		    }
		    //Scanner suppresses exceptions
		    if (sc.ioException() != null) { throw sc.ioException(); }
		    expMgr.dispMessage("fileIOManager:"+owner, "_loadFileIntoStringAra_MT",dispYesStr+"\tLength : " +  count + " distributed into "+lines.length+" arrays.", MsgCodes.info1);
		    for (int i=0;i<lines.length;++i) {res[i] = lines[i].toArray(new String[0]);	 }
		    res[res.length-1]=hdrRes;
		} catch (Exception e) {	
			e.printStackTrace();
			expMgr.dispMessage("fileIOManager:"+owner, "_loadFileIntoStringAra_MT","!!"+dispNoStr, MsgCodes.error2);
			res= new String[0][];
		} 
		finally {
		    if (inputStream != null) {inputStream.close();		    }
		    if (sc != null) { sc.close();		    }
		}
		return res;
	}//_loadFileIntoStringAra_MT

}//class fileIOManager



enum ConsoleCLR {
    RESET("\033[0m"),
    // Regular Colors
    BLACK("\033[0;30m"),    // BLACK
    RED("\033[0;31m"),      // RED
    GREEN("\033[0;32m"),    // GREEN
    YELLOW("\033[0;33m"),   // YELLOW
    BLUE("\033[0;34m"),     // BLUE
    MAGENTA("\033[0;35m"),  // MAGENTA
    CYAN("\033[0;36m"),     // CYAN
    WHITE("\033[0;37m"),    // WHITE
    // Bold
    BLACK_BOLD("\033[1;30m"),   // BLACK
    RED_BOLD("\033[1;31m"),     // RED
    GREEN_BOLD("\033[1;32m"),   // GREEN
    YELLOW_BOLD("\033[1;33m"),  // YELLOW
    BLUE_BOLD("\033[1;34m"),    // BLUE
    MAGENTA_BOLD("\033[1;35m"), // MAGENTA
    CYAN_BOLD("\033[1;36m"),    // CYAN
    WHITE_BOLD("\033[1;37m"),   // WHITE
    // Underline
    BLACK_UNDERLINED("\033[4;30m"),     // BLACK
    RED_UNDERLINED("\033[4;31m"),       // RED
    GREEN_UNDERLINED("\033[4;32m"),     // GREEN
    YELLOW_UNDERLINED("\033[4;33m"),    // YELLOW
    BLUE_UNDERLINED("\033[4;34m"),      // BLUE
    MAGENTA_UNDERLINED("\033[4;35m"),   // MAGENTA
    CYAN_UNDERLINED("\033[4;36m"),      // CYAN
    WHITE_UNDERLINED("\033[4;37m"),     // WHITE
    // Background
    BLACK_BACKGROUND("\033[40m"),   // BLACK
    RED_BACKGROUND("\033[41m"),     // RED
    GREEN_BACKGROUND("\033[42m"),   // GREEN
    YELLOW_BACKGROUND("\033[43m"),  // YELLOW
    BLUE_BACKGROUND("\033[44m"),    // BLUE
    MAGENTA_BACKGROUND("\033[45m"), // MAGENTA
    CYAN_BACKGROUND("\033[46m"),    // CYAN
    WHITE_BACKGROUND("\033[47m"),   // WHITE
    // High Intensity
    BLACK_BRIGHT("\033[0;90m"),     // BLACK
    RED_BRIGHT("\033[0;91m"),       // RED
    GREEN_BRIGHT("\033[0;92m"),     // GREEN
    YELLOW_BRIGHT("\033[0;93m"),    // YELLOW
    BLUE_BRIGHT("\033[0;94m"),      // BLUE
    MAGENTA_BRIGHT("\033[0;95m"),   // MAGENTA
    CYAN_BRIGHT("\033[0;96m"),      // CYAN
    WHITE_BRIGHT("\033[0;97m"),     // WHITE
    // Bold High Intensity
    BLACK_BOLD_BRIGHT("\033[1;90m"),    // BLACK
    RED_BOLD_BRIGHT("\033[1;91m"),      // RED
    GREEN_BOLD_BRIGHT("\033[1;92m"),    // GREEN
    YELLOW_BOLD_BRIGHT("\033[1;93m"),   // YELLOW
    BLUE_BOLD_BRIGHT("\033[1;94m"),     // BLUE
    MAGENTA_BOLD_BRIGHT("\033[1;95m"),  // MAGENTA
    CYAN_BOLD_BRIGHT("\033[1;96m"),     // CYAN
    WHITE_BOLD_BRIGHT("\033[1;97m"),    // WHITE
    // High Intensity backgrounds
    BLACK_BACKGROUND_BRIGHT("\033[0;100m"),     // BLACK
    RED_BACKGROUND_BRIGHT("\033[0;101m"),       // RED
    GREEN_BACKGROUND_BRIGHT("\033[0;102m"),     // GREEN
    YELLOW_BACKGROUND_BRIGHT("\033[0;103m"),    // YELLOW
    BLUE_BACKGROUND_BRIGHT("\033[0;104m"),      // BLUE
    MAGENTA_BACKGROUND_BRIGHT("\033[0;105m"),   // MAGENTA
    CYAN_BACKGROUND_BRIGHT("\033[0;106m"),      // CYAN
    WHITE_BACKGROUND_BRIGHT("\033[0;107m");     // WHITE
	
    private final String value;
	private static Map<String, ConsoleCLR> map = new HashMap<String, ConsoleCLR>(); 
    private ConsoleCLR(String _val){value = _val;}     
	static { for (ConsoleCLR enumV : ConsoleCLR.values()) { map.put(enumV.value, enumV);}}
	public String getVal(){return value;}
	public static ConsoleCLR getVal(String idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	@Override
    public String toString() { return value; }
}//console printout colors enum

//enum delineates each kind of message to be displayed - various information and error codes 
enum MsgCodes{
	info1(0),info2(1),info3(2), info4(3), info5(4),
	warning1(5),warning2(6),warning3(7),warning4(8),warning5(9),
	error1(10),error2(11),error3(12),error4(13),error5(14);
	
	private int value; 
	private static Map<Integer, MsgCodes> map = new HashMap<Integer, MsgCodes>(); 
	static { for (MsgCodes enumV : MsgCodes.values()) { map.put(enumV.value, enumV);}}
	private MsgCodes(int _val){value = _val;} 
	public int getVal(){return value;}
	public static MsgCodes getVal(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	@Override
    public String toString() { return ""+value; }	
}
