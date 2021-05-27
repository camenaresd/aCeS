import java.awt.Desktop;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;
import javafx.scene.control.TextField;
import javafx.scene.control.Label;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.scene.control.RadioButton;
import javafx.scene.control.ToggleGroup;
/**
 *
 * @author camenaresd
 */
public class aces2 extends Application {
    
        private Desktop desktop = Desktop.getDesktop();    
    
    	/**
         * This method was written by CC and modified by DJC, based upon code from StackOverflow user Kip.
         * It takes a string and writes it to a specified file, appending said file.
         * The purpose of this method is to save information and results as they are generated, reducing memory load.
         * 
         * @param content The string to be written
         * @param file The destination file for the string output
         */
	private void writeFileN(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    content += "" + System.lineSeparator();
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(aces2.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            }        
		
	/**
         * This method was written by CC
         * 
         * @param file The desired file to be opened.
         */
        private void openFile(File file) {
            try {
		 desktop.open(file);
		} catch (IOException ex) {
		Logger.getLogger(aces2.class.getName()).log(
                 Level.SEVERE, null, ex
                 );
		 }
            }
        
        // Configures the File Chooser, from Oracle Documentation Example 26-5
	private static void configureFileChooser(
            final FileChooser fileChooser) {      
		fileChooser.setTitle("View Files");
		fileChooser.setInitialDirectory(
			new File(System.getProperty("user.dir"))
			);                 
		fileChooser.getExtensionFilters().addAll(
		new FileChooser.ExtensionFilter("Plain Text", "*.txt"),
		new FileChooser.ExtensionFilter("FASTA", "*.fasta"),
		new FileChooser.ExtensionFilter("Rich Text Format", "*.rtf"),
		new FileChooser.ExtensionFilter("All Files", "*.*")
		);
            }        

     /**
     * A timestamp, used for generating unique file IDs
     */			 
    final long timeUnique = System.currentTimeMillis();
    public String dirName = Long.toString(timeUnique);
    public String filePath = "";

    // NEED TO CREATE VARIABLE OR OBJECT FOR MUTATION MATRIX
    
    @Override
    public void start(Stage primaryStage) {
        primaryStage.setTitle("Generate Simulated MSA");

	// Create Grid pane, FileChooser, and Button
        final GridPane inputGridPane = new GridPane();		
        final FileChooser fileChooser = new FileChooser();	        
        final Button processingButton = new Button("Open MSA Parameters File");	
        final Button idButton = new Button("Generate New Job ID");	
 
 		// Text Fields and Labels. Some of this framework contributed by Christopher Camenares
		Label lbl1 = new Label("Job ID#:");
		lbl1.setMinHeight(50);
		lbl1.setMinWidth(250);                 
                
		Label lbl2 = new Label("Awaiting File Selection");
		lbl2.setMinHeight(50);
		lbl2.setMinWidth(100);               
		
		TextField jobID = new TextField();
		jobID.setText(dirName);
		jobID.setMinHeight(50);
		jobID.setMinWidth(200);                  

		TextField inputZ = new TextField();
		inputZ.setText("Enter Zeta Factor");
		inputZ.setMinHeight(50);
		inputZ.setMinWidth(200);                                  

		TextField inputW = new TextField();
		inputW.setText("Enter Omega Factor");
		inputW.setMinHeight(50);
		inputW.setMinWidth(200);   

		TextField inputE = new TextField();
		inputE.setText("Enter Epsilon Factor");
		inputE.setMinHeight(50);
		inputE.setMinWidth(200);  

		TextField inputM = new TextField();
		inputM.setText("Enter Mu Factor");
		inputM.setMinHeight(50);
		inputM.setMinWidth(200);

		TextField inputG = new TextField();
		inputG.setText("Enter Gamma Factor");
		inputG.setMinHeight(50);
		inputG.setMinWidth(200);
                
		CheckBox chck1;
		chck1 = new CheckBox("Debugging Mode");            
                
                final ToggleGroup outputChoice = new ToggleGroup();
                
		RadioButton rb1 = new RadioButton();
                rb1.setText("Folder named job ID");
                rb1.setSelected(true);
                rb1.setToggleGroup(outputChoice);
                rb1.setUserData("id");
                
		RadioButton rb2 = new RadioButton();
                rb2.setText("Common output folder");
                rb2.setToggleGroup(outputChoice);
                rb2.setUserData("common");
                
       		inputGridPane.add(lbl1, 0, 0);
       		inputGridPane.add(jobID, 0, 1);
       		inputGridPane.add(idButton, 0, 2);
                inputGridPane.add(rb1, 0, 3);
                inputGridPane.add(rb2, 0, 4);                
                inputGridPane.add(chck1, 0, 6);
                inputGridPane.add(inputZ, 0, 7);
                inputGridPane.add(inputW, 0, 8);
                inputGridPane.add(inputE, 0, 9);                     
                inputGridPane.add(inputM, 0, 10);
                inputGridPane.add(inputG, 0, 11);
                inputGridPane.add(lbl2, 0, 12);
       		inputGridPane.add(processingButton, 0, 13);
        
        final Pane rootGroup = new VBox(12);
        rootGroup.getChildren().addAll(inputGridPane);
        rootGroup.setPadding(new Insets(12, 12, 12, 12));
 
        /**
         * Defines button action: generates new job ID
         */
        idButton.setOnAction(
			new EventHandler<ActionEvent>() {
                        @Override
                        public void handle(final ActionEvent e) {
						final long timeButton = System.currentTimeMillis();                   
						jobID.setText(Long.toString(timeButton));
                    }
                }
            );
       
        
       processingButton.setOnAction(new EventHandler<ActionEvent>() {
         @Override
          public void handle(final ActionEvent e) {
           File file1 = fileChooser.showOpenDialog(primaryStage);
                
// Begin processing of file!!

            if (file1 != null) {
             filePath = file1.getPath();
             
             final long timeStart = System.currentTimeMillis();
             
             boolean debugModeOn = chck1.isSelected();
             
             //Defaults
             /*
             Zeta: the amount of mutants to create each round
             */
             int factorZ = 1000;
             
             /*
             Omega: the amount of branching for each round
             */
             int factorW = 1;
             
             /*
             Epsilon: the cutoff for each round
             */
             int factorE = 90;

             /*
             Gamma: the global weight for co-evolution versus conservation
             */
             int factorG = 1;
             
             /*
             Mu: the basic mutation rate. To be replaced, eventually!
             */
             double factorM = 1;
             
             String inputStatus = "";

                // Get the value from input boxes
                try
                {   
                  // the String to int conversion happens here
                  factorZ = Integer.parseInt(inputZ.getText().trim());
                }
                catch (NumberFormatException nfe)
                {
                    inputStatus += "(Input for Zeta not valid, default used) ";
                };

                try
                {   
                  // the String to int conversion happens here
                  factorW = Integer.parseInt(inputW.getText().trim());
                }
                catch (NumberFormatException nfe)
                {
                    inputStatus += "(Input for Omega not valid, default used) ";
                };

                try
                {   
                  // the String to int conversion happens here
                  factorE = Integer.parseInt(inputE.getText().trim());
                  
                  if(factorE > 100)
                  {
                      factorE = 100;
                  }
                  else if(factorE < 0)
                  {
                      factorE = 0;
                  }
                }
                catch (NumberFormatException nfe)
                {
                    inputStatus += "(Input for Epsilon not valid, default used) ";
                };                  

                try
                {   
                  // the String to int conversion happens here
                  factorG = Integer.parseInt(inputG.getText().trim());
                  
                  if(factorG > 100)
                  {
                      factorG = 100;
                  }
                  else if(factorG < 0)
                  {
                      factorG = 0;
                  }
                }
                catch (NumberFormatException nfe)
                {
                    inputStatus += "(Input for Gamma not valid, default used) ";
                };                
                
                try
                {   
                  // the String to int conversion happens here
                  factorM = Double.parseDouble(inputM.getText().trim());

                  if(factorM > 100.0)
                  {
                      factorM = 100.0;
                  }
                  else if(factorM < 0.0)
                  {
                      factorM = 0.0;
                  }                  
                }
                catch (NumberFormatException nfe)
                {
                    inputStatus += "(Input for Mu not valid, default used) ";
                };
                
             
             // Get new name for directory, initialize directory
             String jobIDtxt = jobID.getText();
             String dirName = "";
             if(outputChoice.getSelectedToggle().getUserData().toString().equals("common")){
                 dirName = "output";
             }
             else
             {
                 dirName = jobIDtxt;
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
              
             molecule2[] molList = new molecule2[5];
             
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
             
             ArrayList<mutual2> miPairs = new ArrayList<mutual2>();
             
             try (BufferedReader br = new BufferedReader(new FileReader(file1))) {
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
                         This block determines that this is the header for a new molecule2 and gathers information accordingly.
                         This creates a new molecule2 object
                      */
                      molDefine = true;
                      miDefine = false;
                      line = line.replaceAll(">", "");                         
                      String molLet = line.substring(0,1);
                      setMol = molNumber(molLet);
                      line = line.replaceFirst(molLet, "");
                      String[] sublines = line.split(": ");
                      int molLen = Integer.parseInt(sublines[0]);
                      molList[setMol] = new molecule2(molLen, sublines[1], setMol);
                      catNumList[setMol] = molList[setMol].getCatSize();
                      
                      if (debugModeOn)
                      {
                          molList[setMol].debugSet(debugLog);
                          writeFileN("New Molecule Created: " + molLet + "\t Lenght: " + molLen + "\t Categories: " + sublines[1], debugLog);
                      }
                     }
                    else if(molDefine)
                    {
                      /*
                        This block contains instructions on how to pull information for each line, regardless of what molecule2 it is going to.
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
                     // I NEED TO FIX THIS PART TOO!
                       /*
                        0 - Molecule
                        1 - Residue
                        2 - Interacting molecule2
                        3 - Interacting residue
                        4 - Mutual Information (100%) target      
                      */
                      String[] sL = line.split(", ");
                      int[] sLv = strArrIntArr(sL);   
                      
                      // This creates a new mutual2 information distribution given a particular set of constraints
                      mutual2 miTemp = new mutual2(miSettings, molList[sLv[0]], sLv[1]-1, molList[sLv[2]], sLv[3]-1, sLv[4], orgNumber);

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
              Logger.getLogger(aces2.class.getName()).log(Level.SEVERE, null, ex);
             } 
             catch (IOException ex) {
              Logger.getLogger(aces2.class.getName()).log(Level.SEVERE, null, ex);
             }                              

             /*
             Code to run after all the information has been collected; Run through each organism, make sequence and print
             
             The code to follow uses several steps of logical to select a particular residue identity
             It starts by looping through all organisms
             */
             
             int[][] orgOmePrev = new int[5][];             
             
                   // Cycle through all 5 molecules                    
                   for (int j = 0; j < molList.length; j++)
                   {   
                      int thisMolSize = molList[j].getSize();
                      orgOmePrev[j] = new int[thisMolSize];   
 
                      // Cycle through all residues in the molecule2
                       for (int k = 0; k < thisMolSize; k++)
                       {    
                          int theChoice = 0;

                        //Use most common residue for top choice
                        theChoice = molList[j].chooseTopResidue(k);
                        orgOmePrev[j][k] = theChoice;
                        // Residue code block ends
                      } 
                   }
             
             
             for (int i = 0; i < orgNumber; i++)
              {
                  int[][][] orgZ = new int[factorZ][][];
                  int[] orgZscore1 = new int[factorZ];
                  int[] orgZscore2 = new int[factorZ];

                  int maxScore = 0;             
                  
                   // THIS PART IS CHANGED FOR ACES2
                   // INITIALIZE TEMPORARY ARRAYS
                   
                   // Cycle through iterations of making the molecule2
                for (int q = 0; q < factorZ; q++){
                    
                 int[][] orgOmeInt = new int[5][];
                 
                 // Initialize the orgOme array to allow for call aheads later
                   for (int j = 0; j < molList.length; j++)
                   {
                      int thisMolSize = molList[j].getSize();
                      orgOmeInt[j] = new int[thisMolSize];                   
                   }    
                                      
                   // Cycle through all 5 molecules                    
                   for (int j = 0; j < molList.length; j++)
                   {   
                      /*
                       This gives the size of a particular molecule2's length
                       */ 
                      int thisMolSize = molList[j].getSize();  
                       
                      // Cycle through all residues in the molecule2
                       for (int k = 0; k < thisMolSize; k++)
                       {    
                         
                          orgOmeInt[j][k] = orgOmePrev[j][k];
                          //MUTATION BLOCK
                          
                          int catSize = molList[j].getCatSize();
                          double random100d = (Math.random() * 100.0); 

                          // Make a decision if a mutation should be made. This is done by comparing the mutation rate versus the current observed residues
                          // In the future, factorM will be an aggregate of all mutation rates for a particular residue if not otherwise the default
                          if(Double.compare(factorM, random100d) >= 0 || molList[j].getObs(k)[orgOmeInt[j][k]] < 1)
                          {
                             // It has been decided that a mutation is needed. Now, select one that is not the same, and also allowed given previous constraints.
                             // Which choice is made is completely random - this can be addressed with some sort of mutation table? Otherwise, it is Jukes-Cantor
                             int[] randomChoice = randomWalk1D(catSize);      
                            
                             for (int p = 0; p < randomChoice.length; p++)
                             {
                                 int newChoice = randomChoice[p];
                                 
                                 if(newChoice != orgOmePrev[j][k] || molList[j].getObs(k)[newChoice] > 0)
                                 {
                                     orgOmeInt[j][k] = newChoice;
                                     break;
                                 }
                             }
                          }
                          
                          // Residue code block ends
                      }

                    // Molecule code block ends
                  }
                   //All organisms for that mutant are created - now it's time to score it!
                   int score = 0;

                   
                   for(int j = 0; j < orgOmeInt.length; j++)
                   {
                       for (int k = 0; k < orgOmeInt[j].length; k++)
                       {
                          int thisResidue = orgOmeInt[j][k];
                       
                          if(molList[j].checkForMI(k))
                          {
                         
                             int[] resMIpairs = molList[j].getMIpairPos(k);
                             // This gives the current distributions of residues across the MSA for that residue (?)
                             
                             for (int m = 0; m < resMIpairs.length; m++)
                             {
                                 
                              // Cycle through all MI Pairs for this residue
                                     int accessNum = resMIpairs[m] - 1;

                                     int[] tNa = miPairs.get(accessNum).getNames();
                                     int[] tNu = miPairs.get(accessNum).getNums();
                                     int[][] obsMIz = miPairs.get(accessNum).getMI();
                                     
                                     int mol1Choice = orgOmeInt[tNa[0]][tNu[0]];
                                     int mol2Choice = orgOmeInt[tNa[1]][tNu[1]];
                                     
                                     /* For some reason, the following was not helpful.
                                     if(k == tNa[1])
                                     {
                                     mol1Choice = orgOmeInt[tNa[1]][tNu[1]];
                                     mol2Choice = orgOmeInt[tNa[0]][tNu[0]];                                          
                                     }
                                     */
                                     
                                     score += (obsMIz[mol1Choice][mol2Choice] * factorG);
      
                                     score += molList[j].getObs(k, thisResidue);
                                     
                                // End code block for loop m, MI pairs
                              }
                            }
                          else
                          {
                              score += molList[j].getObs(k, thisResidue);
                          }
                         }
                     }
                   
                     maxScore = Math.max(maxScore, score);

if(debugModeOn){writeFileN("For organism : " + i + ", mutant species " + q + " created with score of " + score, debugLog);}                   
                     orgZscore1[q] = score;
                     orgZscore2[q] = score;                     
                     orgZ[q] = orgOmeInt;                     
                }
               // PROCESS TEMPORARY ARRAYS - ELMINATE BASED ON CONSTRAINTS, SELECT W CANDIDATES
                int[] coords1 = new int[orgZscore1.length];
                int[] coords2 = new int[orgZscore1.length];
                
                Arrays.sort(orgZscore1);
                
                for (int z1 = 0; z1 < orgZscore1.length; z1++)
                {
                    for (int z2 = 0; z2 < orgZscore2.length; z2++)
                    {
                        // A potential issue here - same score across two sequences. Very rare?
                        if(orgZscore1[z1] == orgZscore1[z2])
                        {
                            coords2[z1] = z2;
                            coords1[z2] = z1;
                        }
                    }
                }

                for (int w2 = 0; w2 < factorW; w2++)
                {
                    int topCrop = orgZ.length - (orgZ.length * factorE/100);
                    int randomSelection = (int)Math.round(Math.random() * Math.max(topCrop, factorW));
if(debugModeOn){writeFileN("RandomSelection: " + randomSelection + ", score: " + orgZscore2[randomSelection] + ", organism: " + i, debugLog);}
                    int[][] orgOmeInt = orgZ[randomSelection];

                    // Convert residue identity to characters

                    for (int w = 0; w < orgOmeInt.length; w++)
                    {
                        String sequence = "";
                        for (int z = 0; z < orgOmeInt[w].length; z++)
                        {
                            int groupChar = molList[w].getInGroupChoice(z, orgOmeInt[w][z]);
                            
                          //subtract selected sequences from obs
                          if(molList[w].checkForMI(z))
                          {
                             int[] resMIpairs = molList[w].getMIpairPos(z);
                             // This gives the current distributions of residues across the MSA for that residue (?)
                             
                             for (int m = 0; m < resMIpairs.length; m++)
                             {
                                 
                              // Cycle through all MI Pairs for this residue
                                     int accessNum = resMIpairs[m] - 1;

                                     int[] tNa = miPairs.get(accessNum).getNames();
                                     int[] tNu = miPairs.get(accessNum).getNums();
                                     
                                     int mol1Choice = orgOmeInt[tNa[0]][tNu[0]];
                                     int mol2Choice = orgOmeInt[tNa[1]][tNu[1]];                                     

                                     miPairs.get(accessNum).useMI(mol1Choice, mol2Choice);
                                     
                                // End code block for loop m, MI pairs
                              }
                            }
                          else
                           {
                            molList[w].useObs(z, orgOmeInt[w][z]);
                           }
                          
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
                   orgOmePrev = orgOmeInt;
                   
                   i += w2;
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
                writeFileN("Zeta: " + factorZ + ", Omega: " + factorW + ", Epsilon: " + factorE + ", Mu: " + factorM, runTime);
                writeFileN(inputStatus, runTime);
                
//                writeFileN(, runTime);

                openFile(dir);    
            }
          }
        });
        
        Scene scene = new Scene(rootGroup, 350, 550);
        
        primaryStage.setTitle("aCeS2 - beta version");
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        launch(args);
    }

    /**
     * Used for converting a molecule2 name assignment (A through E) to array positions
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
    public int checkInteger(String input, int upperLimit, int lowerLimit)
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
    
    public int[] strArrIntArr(String[] array)
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
    
    public String printArray(int[][] array)
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
    
    public void printArray(String[] array)
    {
            for (int j = 0; j < array.length; j++)
            {
                System.out.print(array[j] + "\t");
            }
                System.out.println("***");        
    }

    public String printArray(int[] array)
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
    
    public boolean logicXOR(boolean x, boolean y)
    {
        return ((x || y) && !(x && y));
    }
    
    public int[] getColumn(int[][] arr, int column)
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
    
    public int[][] obsProb (int[][] arr){
        int sum = 0;
        int[][] res = new int[arr.length][arr[0].length];

        for (int i = 0; i < arr.length; i++)
        {
            for (int j = 0; j < arr[i].length; j++)
            {
                sum += arr[i][j];
            }
        }
        
        if(sum <= 0)
        {
            sum = 1;
        }
        
        for (int i = 0; i < arr.length; i++)
        {
            for (int j = 0; j < arr[i].length; j++)
            {
                arr[i][j] *= (100/sum);
                if(arr[i][j] > 100)
                {
                    arr[i][j] = 100;
                }
            }
        }
        return res;
    }
    
public static int[] randomWalk1D(int size)
    {
        int[] coords = new int[size];
        
        double[] rando1 = new double[size];
        double[] rando2 = new double[size];
        
        for (int i = 0; i < rando1.length; i++){
            double randoNum = Math.random();
            rando1[i] = randoNum;
            rando2[i] = randoNum;
        }
        
        // Sort one of the random array clones
        Arrays.sort(rando2);

        for (int i = 0; i < coords.length; i++){
            for (int j = 0; j < rando1.length; j++){
                if (rando1[j] == rando2[i])
                {
                    coords[i] = j;
                }
            }
        }
        
        return coords; 
        
    }    

}
