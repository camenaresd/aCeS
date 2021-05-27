
import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Devin
 */
public class aces1cl {

    /**
     * @param args the command line arguments
     */
    @SuppressWarnings("empty-statement")
    public static void main(String[] args)
    {

     Scanner in = new Scanner(System.in);
     String fileName = "";
     boolean done = false;

     while(!done)
      {
       System.out.print("Enter the  filename:");
       fileName = in.nextLine();

       if(fileName.equals("") || fileName == null)
        {
         done = true;
         System.out.println("I'm done!");
         break;
        }
       
        /**
         * A timestamp, used for generating unique file IDs
         */			 
        final long timeUnique = System.currentTimeMillis();
        String dirName = Long.toString(timeUnique);
        String jobIDtxt = dirName;
        String filePath = "";  
        
        final long timeStart = System.currentTimeMillis();
        
        int biSimpleWint = 0;
        String biSimpleStatus = "";
        boolean biSimpleTF = false;
        
        boolean debugModeOn = false;        

       /**
        * This section will take parameters as arguments for the rest of the program. 
        */
       if(fileName.contains(", "));
       {
           String[] parameters = fileName.split(", ");
           fileName = parameters[0];
           
           for (int i = 1; i < parameters.length; i++)
           {
               String[] subpar = parameters[i].split(":");
               
               if(subpar[0].contains("id"))
               {
                   jobIDtxt = subpar[1];
               }
               else if(subpar[0].contains("output"))
               {
                   dirName = subpar[1];
               }
               else if(subpar[0].contains("phylogeny"))
               {
                    try{biSimpleWint = Integer.parseInt(subpar[1]);}
                    catch(Exception e){biSimpleWint = 0; biSimpleTF = false;};

                    if(biSimpleWint > 100)
                    {
                        biSimpleWint = 100;
                        biSimpleTF = true;
                    }
                    else if (biSimpleWint < 0)
                    {
                        biSimpleWint = 0;
                        biSimpleTF = false;
                    }                   
               }
               else if(subpar[0].contains("debug"))
               {
                   if(subpar[1].equals("t"))
                   {
                       debugModeOn = true;
                   }
               }               
           }
         }        

        File dir = new File(dirName);
        dir.mkdir();         
       
             //Initalize final files
             /**
              * Recipient file for runtime and processing information
              */
             File runTime = new File(dirName + "\\Runtime_Info_" + jobIDtxt + ".txt");                                          

             File debugLog = new File(dirName + "\\Debug_log_" + jobIDtxt + ".txt");

             File molA = new File(dirName + "\\Molecule_A_" + jobIDtxt + ".txt");                                        
             File molB = new File(dirName + "\\Molecule_B_" + jobIDtxt + ".txt");                                        
             File molC = new File(dirName + "\\Molecule_C_" + jobIDtxt + ".txt");     
             File molD = new File(dirName + "\\Molecule_D_" + jobIDtxt + ".txt");     
             File molE = new File(dirName + "\\Molecule_E_" + jobIDtxt + ".txt");                  
              
             moleculeCL[] molList = new moleculeCL[5];
             
             /*
             The number of organisms or species for which there will be sequence
             */
             int orgNumber = 0;
             
             int setMol = 0;
             
             int[] catNumList = new int[5];
             
             boolean analyzeLine = false;
             
             boolean molDefine = false;
             
             boolean miDefine = false;
             
             int[] miSettings = new int[3];
             
             int lineCounter = 0;
             
             int miCounter = 0;
             
             ArrayList<mutualCL> miPairs = new ArrayList<mutualCL>();        
        
       try (BufferedReader br = new BufferedReader(new FileReader(fileName))){
                  String line;
                  while ((line = br.readLine()) != null) {
                
                lineCounter++;  
                  
                /*
                  Allows reader to skip ahead if line contains @@ symbols. Nothing to commit to memory
                */  
                if (line.contains("@@START"))
                 {
                  analyzeLine = true;
                  line = br.readLine();
                 }
                else if (line.contains("@@END"))
                {
                  analyzeLine = false;
                  molDefine = false;
                  miDefine = false;
                  line = br.readLine();                    
                }
                
                if (analyzeLine)
                {
                    /*
                    Start of processing line that has useful information. Recognize what element of the line it is.
                    */
                    if (line.contains("##"))
                     {
                      /*
                      In the input file, the first block of instructions after @@START will contain the number of organisms in the MSA
                      followed by information on how to create the underlying molecules (how much conservation, what residues, etc).
                      */
                      line = line.replaceAll("##", "");
                      orgNumber = Integer.parseInt(line);
                     }               
                    else if (line.contains(">"))
                     {
                      /*
                         This block determines that this is the header for a new molecule and gathers information accordingly.
                         This creates a new molecule object
                      */
                      molDefine = true;
                      miDefine = false;
                      line = line.replaceAll(">", "");                         
                      String molLet = line.substring(0,1);
                      setMol = molNumber(molLet);
                      line = line.replaceFirst(molLet, "");
                      String[] sublines = line.split(": ");
                      int molLen = Integer.parseInt(sublines[0]);
                      molList[setMol] = new moleculeCL(molLen, sublines[1], setMol);
                      catNumList[setMol] = molList[setMol].getCatSize();
                      molList[setMol].biSimpleSet(biSimpleTF, biSimpleWint);
                      

                      if (debugModeOn)
                      {
                          molList[setMol].debugSet(debugLog);
                          writeFileN("New Molecule Created: " + molLet + "\t Lenght: " + molLen + "\t Categories: " + sublines[1], debugLog);
                      }

                     }
                    else if(molDefine)
                    {
                      /*
                        This block contains instructions on how to pull information for each line, regardless of what molecule it is going to.
                        First step is to split the line into component parts, sL.

                        Residue information indicated in this fashion: 
                        0 - Residue Number
                        1- Identity (group)
                        3- Identity (residue)
                        4- Conservation
                       # characters (or any non-integer) for a particular value indicate wildcard.          
                      */
                      String[] sL = line.split(", ");
                      
                      int minCons1 = (int)Math.ceil(100.0 / catNumList[setMol]);
                      
                      /*
                      Process basic information that should be present in each line (at least three)
                      */
                      int resNo = checkInteger(sL[0], molList[setMol].getSize(), 1) - 1;         
                      int resIdent1 = checkInteger(sL[1], catNumList[setMol], 0);                      
                      int resIdent2 = checkInteger(sL[2], molList[setMol].getCatSize(resIdent1), 0);
                      int resCons1 = checkInteger(sL[3], 101, minCons1);

                      int minCons2 = (int)Math.ceil(100 / molList[setMol].getCatSize(resIdent1));
                      int resCons2 = checkInteger(sL[4], 101, minCons2);
                      
                      // This function will set the information for each residue
                      molList[setMol].setRes(resNo, resIdent1, resIdent2, resCons1, resCons2); 
                      
                    }
                    else if(line.contains("MI#"))
                    {
                        molDefine = false;
                        miDefine = true;
                        
                        line = line.replaceAll("MI#", "");
                        String[] sublinesInt = line.split(", ");
                        miSettings = strArrIntArr(sublinesInt);
                        
                        /*
                        Block of code to populate the molecules
                        */


                        if (debugModeOn)
                        {
                          writeFileN("Fill in unspecified positions:", debugLog);
                        }

                        
                        for (int i = 0; i < molList.length; i++)
                         {
                          molList[i].setOrg(orgNumber);
                          molList[i].fillProb();
                         }                                     
                        
                    }
                    else if(miDefine = true)
                    {
                       /*
                        0 - Molecule
                        1 - Residue
                        2 - Interacting molecule
                        3 - Interacting residue
                        4 - Mutual Information (100%) target      
                      */
                      String[] sL = line.split(", ");
                      int[] sLv = strArrIntArr(sL);   
                      
                      // This creates a new mutual information distribution given a particular set of constraints
                      mutualCL miTemp = new mutualCL(miSettings, molList[sLv[0]], sLv[1]-1, molList[sLv[2]], sLv[3]-1, sLv[4], orgNumber);


                      if (debugModeOn)
                      {
                          miTemp.debugSet(debugLog);
                      }

                      
                      miCounter++;
                      
                      if(sLv[0] == sLv[2])
                      {
                         molList[sLv[0]].setMIpair(Math.max(sLv[1]-1, sLv[3]-1), miCounter);                         
                      }
                      else
                      {
                      molList[sLv[0]].setMIpair(sLv[1]-1, miCounter);
                      molList[sLv[2]].setMIpair(sLv[3]-1, miCounter);                          
                      }
                      
                      miPairs.add(miTemp);

                    }
                   }
              
       
                      
                      
                  }   
                 } 
                 catch (FileNotFoundException ex) {
                  // Logger.getLogger(aces1cl.class.getName()).log(Level.SEVERE, null, ex);
                  System.out.println("Invalid filename selected or file not found. Try again!");
                 } 
                 catch (IOException ex) {
                  Logger.getLogger(aces1cl.class.getName()).log(Level.SEVERE, null, ex);
                 }
       
 /*
             Code to run after all the information has been collected; Run through each organism, make sequence and print
             
             The code to follow uses several steps of logical to select a particular residue identity
             It starts by looping through all organisms
             */
             
             for (int i = 0; i < orgNumber; i++)
              {
                 int[][] orgOmeInt = new int[5][];
                 boolean[][] orgOmeFill = new boolean[5][];
                 
                 // Initialize the orgOme array to allow for call aheads later
                   for (int j = 0; j < molList.length; j++)
                   {
                      int thisMolSize = molList[j].getSize();
                      orgOmeInt[j] = new int[thisMolSize];
                      orgOmeFill[j] = new boolean[thisMolSize];                      
                   }    
                  
                   
                   // Cycle through all 5 molecules
                   for (int j = 0; j < molList.length; j++)
                   {   
                      /*
                       This gives the size of a particular molecule's length
                       */ 
                      int thisMolSize = molList[j].getSize();  
                       
                      // Cycle through all residues in the molecule
                       for (int k = 0; k < thisMolSize; k++)
                       {    
                          int theChoice;
                          
                          // Check to see if there is mutual information at all for this residue
                          if(molList[j].checkForMI(k))
                          {
                             int[] resMIpairs = molList[j].getMIpairPos(k);
                             
                             // This gives the current distributions of residues across the MSA for that residue (?)
                             int[] obsK = molList[j].getObs(k);
                            
                             int[] probSumMI = new int[obsK.length];
                             
                             boolean useMIarr = false;

                             int runningTotal = 0;
                             
                             for (int m = 0; m < resMIpairs.length; m++)
                             {
                              int[] obsKmi = new int[obsK.length];    
                                 
                              // Cycle through all MI Pairs for this residue
                                     int accessNum = resMIpairs[m] - 1;

                                     int[] tNa = miPairs.get(accessNum).getNames();
                                     int[] tNu = miPairs.get(accessNum).getNums();
                                     
                                     // Code in case one is filled but the other isn't
                                     if(logicXOR(orgOmeFill[tNa[0]][tNu[0]], orgOmeFill[tNa[1]][tNu[1]]))
                                      {
                                        useMIarr = true;
                                        int[][] obsXY = miPairs.get(accessNum).getMI();
                                        
                                        int previousChoice = 0;
                                        int orientation = 0;
                                         if(orgOmeFill[tNa[0]][tNu[0]])
                                         {
                                          previousChoice = orgOmeInt[tNa[0]][tNu[0]];
                                          obsKmi = obsXY[previousChoice];
                                          orientation = 0;
                                         }
                                         else
                                         {
                                          orientation = 1;
                                          previousChoice = orgOmeInt[tNa[1]][tNu[1]];                                         
                                          obsKmi = getColumn(obsXY, previousChoice);
                                         }
                                         
                                        for (int r = 0; r < obsKmi.length; r++)
                                         {
                                          runningTotal += obsKmi[r];
                                          probSumMI[r] += runningTotal;
                                         }                                         
                                         

                                         if(debugModeOn)
                                         {
                                          writeFileN("MI needed! Molecule: " + tNa[0] + "\t Residue: " + tNu[0] + "\t Pair Mol: " + tNa[1] + "\t Pair Res: " + tNu[1] + "\t Orientation: " + orientation + "\t Previous Choice: " + previousChoice, debugLog);
                                          writeFileN("\t obsKMI: " + printArray(obsKmi), debugLog);
                                         }

                                      }                          

                                // End code block for loop m, MI pairs
                                }                            
                             if(useMIarr)
                             {
                                for (int r = 0; r < probSumMI.length; r++)
                                 {
                                    if(runningTotal > 0)
                                    {
                                    probSumMI[r] = (100*probSumMI[r])/runningTotal;   
                                    }                                          
                                 }                                 
                                 
                              theChoice = molList[j].chooseResidue(k, probSumMI);                                   
                             }
                             else
                             {
                              theChoice = molList[j].chooseResidue(k);  
                             }
                            }
                           else
                            {
                             theChoice = molList[j].chooseResidue(k);                              
                            }                    
                          orgOmeInt[j][k] = theChoice;
                          orgOmeFill[j][k] = true;
                        // Residue code block ends
                      }

                    // Molecule code block ends
                  }
               
               // Convert residue identity to characters

               for (int w = 0; w < orgOmeInt.length; w++)
               {
                   String sequence = "";
                   for (int z = 0; z < orgOmeInt[w].length; z++)
                   {
                       int groupChar = molList[w].getInGroupChoice(z, orgOmeInt[w][z]);
                       sequence += molList[w].convertResidue(orgOmeInt[w][z], groupChar);
                   }
                   
                   switch (w)
                   {
                    case 0: writeFileN(">Organism_" + i + "Molecule_A", molA);
                            writeFileN(sequence + System.lineSeparator(), molA);
                     break;
                    case 1: writeFileN(">Organism_" + i + "Molecule_B", molB);
                            writeFileN(sequence + System.lineSeparator(), molB);                  
                     break;
                    case 2: writeFileN(">Organism_" + i + "Molecule_C", molC);
                            writeFileN(sequence + System.lineSeparator(), molC);                    
                     break;
                    case 3: writeFileN(">Organism_" + i + "Molecule_D", molD);
                            writeFileN(sequence + System.lineSeparator(), molD);                    
                     break;
                    case 4: writeFileN(">Organism_" + i + "Molecule_E", molE);
                            writeFileN(sequence + System.lineSeparator(), molE);                    
                     break;                     
                   }
                   
               }            
              // Organism code block ends                  
              }

		final long timeEnd = System.currentTimeMillis();
		long runMinutesL = (timeEnd - timeStart) / 60000;
		int runMinutesI = (int) runMinutesL;             
             
                writeFileN("Job ID#: " + jobIDtxt, runTime);
                writeFileN("Input file " + filePath, runTime);                
		writeFileN(Long.toString(timeEnd - timeStart) + " milliseconds of processing time", runTime);
                writeFileN("- or about " + Integer.toString(runMinutesI) + " minutes", runTime);
                writeFileN("Weight for bifurcating tree was " + biSimpleWint + " " + biSimpleStatus, runTime);
                
                
//                writeFileN(, runTime);

                System.out.println("Completed!");  
                   
       
       }
    }

    /**
     * Used for converting a molecule name assignment (A through E) to array positions
     * @param line
     * @return 
     */
    public static int molNumber(String letter)
    {
        int pos = 0;
        switch(letter)
        {
            case "A": pos = 0;
                    break;
            case "B": pos = 1;
                    break;                            
            case "C": pos = 2;
                    break;        
            case "D": pos = 3;
                    break;        
            case "E": pos = 4;
                    break;        
        }
        return pos; 
    }
    
    public static String molNumber(int num)
    {
        String res = "";
        switch(num)
        {
            case 0: res = "A";
                    break;
            case 1: res = "B";
                    break;                            
            case 2: res = "C";
                    break;        
            case 3: res = "D";
                    break;        
            case 4: res = "E";
                    break;        
        }
        return res; 
    }    
    
    /**
     * Checks to see if the input string is an integer - uses it if it is good and within range, otherwise it will choose randomly according to boundary set.
     * @param input
     * @param mol1
     * @return 
     */
    public static int checkInteger(String input, int upperLimit, int lowerLimit)
     {
        int value = 0;
        double rando = Math.random();
        double newValue = (upperLimit - lowerLimit) * rando; 
        
        try
        {
         value = Integer.parseInt(input);  
          if(value < upperLimit && value > lowerLimit)
           {             
            return value;
           }
          else
           {
            value = (int)newValue;
            value += lowerLimit;         
            return value;             
           }
        }
        catch(Exception e)
        {
         value = (int)newValue;
         value += lowerLimit;         
         return value;
        }

    }    
    
    public static int[] strArrIntArr(String[] array)
    {
        int[] res = new int[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            String tS = array[i].replaceAll(" ", "");
            try{
            res[i] = Integer.parseInt(tS);                   
            }
            catch (NumberFormatException e) {
            res[i] = molNumber(tS);
            };

        }
        
        return res;
    }
    
    public static String printArray(int[][] array)
    {
        String res = "";
        for (int i = 0; i < array.length; i++)
        {
           res += printArray(array[i]);
           res += System.lineSeparator();
        }
        res += "***";        
        return res;
    }
    
    public static void printArray(String[] array)
    {
            for (int j = 0; j < array.length; j++)
            {
                System.out.print(array[j] + "\t");
            }
                System.out.println("***");        
    }

    public static String printArray(int[] array)
    {
            String res = "";
            int sum1 = 0;
            for (int j = 0; j < array.length; j++)
            {
                res += array[j];
                res += "\t";
                sum1 += array[j];
            }
                res += "* ";
                res += sum1;
            return res;
    }
    
    public static boolean logicXOR(boolean x, boolean y)
    {
        return ((x || y) && !(x && y));
    }
    
    public static int[] getColumn(int[][] arr, int column)
    {
        ArrayList<Integer> temp = new ArrayList<>();
        
        for (int i = 0; i < arr.length; i++)
        {
            temp.add(arr[i][column]);
        }
        
        int[] res = new int[temp.size()];
        
        for (int i = 0; i < temp.size(); i++)
        {
            res[i] = temp.get(i);
        }
        
        return res;
    }    

    	/**
         * This method was written by CC and modified by DJC, based upon code from StackOverflow user Kip.
         * It takes a string and writes it to a specified file, appending said file.
         * The purpose of this method is to save information and results as they are generated, reducing memory load.
         * 
         * @param content The string to be written
         * @param file The destination file for the string output
         */
	private static void writeFileN(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    content += "" + System.lineSeparator();
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(aces1cl.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            }        

        private static Desktop desktop = Desktop.getDesktop();            
        
	/**
         * This method was written by CC
         * 
         * @param file The desired file to be opened.
         */
        private static void openFile(File file) {
            try {
		 desktop.open(file);
		} catch (IOException ex) {
		Logger.getLogger(
                 aces1cl.class.getName()).log(
                 Level.SEVERE, null, ex
                 );
		 }
            }
    
    
}
