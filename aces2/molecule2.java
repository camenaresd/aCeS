
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author camenaresd
 */
public class molecule2{
    
    /*
    Indicates how many sequences to generate
    */
    int species = 0;

    /*
    A 2D Array that contains the different categories and residue identities
    */    
    public char[][] molCat;
    
    /*
    An array of probabilities for conservation of residues. 
    */
    public int[][] probA;

    /*
    An array of probabilities for conservation of residues, for identities. 
    */
    public int[][] probB;

    /*
    An array of probabilities for conservation of residues, for identities. 
    */
    public int[][] groupCharTotals;
    
    /*
    An array of probabilities for conservation of residues. Summated probabilities for proper choosing
    */    
    public int[][] probSumA;

    /*
    An array of occurences for residues. 
    */
    public int[][] obsA;     
    
    /*
    An array describing which residues have been filled with user information
    */
    public boolean[] filledA;

    /*
    An array to store the prescence of mutual information interactions
    */
    public boolean[] hasMI;
    
    public File debugMol;
    
    public boolean debugMode; 
    
    /*
    An array to store position information for a mutual information object, per residue of the molecule2.
    First list is for the residue in molecule2 A
    Second list is for the array position for the MI object
    */
    public ArrayList<Integer> miPairPos1 = new ArrayList<>();
    public ArrayList<Integer> miPairPos2 = new ArrayList<>();
    
    
    /*
    An array of selection of residue choice within a specific category. 
    For mutual information purposes, this is unimportant but helps to define specifics about the sequence for the users benefit
    */
    public int[][] residueChoice;
    
    /*
    The lenght of the molecule2
    */
    public int molSize;
    
    /*
    The name or ID of the molecule2
    */
    public int molName;    
    
    /*
    Initialize the molecule2 object, the lenght of the sequence and the categories possible
    */
    public molecule2(int len, String categories, int givenName){
        
        String[] molCat1 = categories.split(", ");        
	molCat = new char[molCat1.length][];  
	for (int i = 0; i < molCat1.length; i++)
	 {
	   molCat[i] = molCat1[i].toCharArray();
	 }

        probA = new int[len][molCat1.length];
        probB = new int[len][molCat1.length];
        probSumA = new int[len][molCat1.length];
        obsA = new int[len][molCat1.length];        
        residueChoice = new int[len][molCat1.length];
        filledA = new boolean[len];
        hasMI = new boolean[len];
        molSize = len;
        molName = givenName;
    }
    
    public void debugSet(File filename)
    {
        debugMol = filename;
        debugMode = true;
    }

    
    public int getName()
    {
        return molName;
    }
    
    /*
    Initially set the number of organisms (number of sequences) that will be generated
    */
    public void setOrg(int num)
    {
        species = num;

        for (int i = 0; i < obsA.length; i++)
         {
            for (int j = 0; j < obsA[i].length; j++)
            {
                obsA[i][j] *= (int)Math.round(species/100);
            }
         }
        
    }
    
    /**
     * A function to store important information about mutual information interactions
     * @param residue The residue of this molecule2 that experiences MI
     * @param molecule The interacting molecule2
     * @param arrayPosition  The array position for the MI object
     */
    public void setMIpair(int residue, int arrayPosition)
    {
        hasMI[residue] = true;
        miPairPos1.add(residue);
        miPairPos2.add(arrayPosition);
    }

    public boolean checkForMI(int residue)
    {
        return hasMI[residue];
    }
    
    /*
    For every position in molecule 1 for a particular pair, deliver the corresponding position in molecule 2
    */
    public int[] getMIpairPos(int residue)
    {
        ArrayList<Integer> subArr = new ArrayList<>();
        
        for (int i = 0; i < miPairPos1.size(); i++)
        {
            if(miPairPos1.get(i) == residue)
            {
                subArr.add(miPairPos2.get(i));
            }
        }
        
        int[] res = new int[subArr.size()];
        
        for (int i = 0; i < subArr.size(); i++)
        {
            res[i] = subArr.get(i);
        }
        return res;
    }

    /*
    A method to call for the category array for a molecule2
    */
    public int getCatSize()
    {
        return molCat.length;
    }
    
    /*
    A method to call for a residue array within a category for a molecule2
    */
    public int getCatSize(int subArr)
    {
        return molCat[subArr].length;
    }    
    
    /*
    A method to call for the category array for a molecule2
    */
    public int getSize()
    {
        return molSize;
    }    

    public int getResCon(int res)
    {
        int max = 0;
        
        for (int i = 0; i < probA[res].length; i++)
        {
            if(probA[res][i] > max)
            {
                max = probA[res][i];
            }
        };
        
        return max;     
    }
    
    /*
    A method to call for the array of instances for a molecule2
    */
    public int[][] getObs()
    {
        return obsA;
    }
    
    /*
    A method to call for the array of instances for a molecule2
    */
    public int[] getObs(int res)
    {
        return obsA[res];
    }

    /*
    A method to call for the instances for a molecule2
    */
    public int getObs(int res, int choice)
    {
        return obsA[res][choice];
    }

        /*
    A method to call for the instances for a molecule2
    */
    public void useObs(int res, int choice)
    {
        if(debugMode)
        {
            writeFileN("Mol: " + molName + ", res: " + res + ".", debugMol);
            String temp = "";
            for(int i = 0; i < obsA[res].length; i++)
            {
                temp += obsA[res][i] + "";
                if(i == choice){
                    temp += "!" + "\t";
                }
                else
                {
                    temp += "\t";
                }
            }
            writeFileN(temp, debugMol);

            writeFileN("(ProbA)", debugMol);
            temp = "";
            for(int i = 0; i < probA[res].length; i++)
            {
                temp += probA[res][i] + "";
                if(i == choice){
                    temp += "!" + "\t";
                }
                else
                {
                    temp += "\t";
                }
            }
            writeFileN(temp, debugMol);            
            
        }
        
        obsA[res][choice]--;
        
    }
    
    public int getProb(int res, int choice)
    {
        return probA[res][choice];
    }       
    
    /**
     * Given information pulled from main class, set this information in the appropriate array
     * @param num
     * @param cat1
     * @param cat2
     * @param cons 
     */
    public void setRes(int num, int cat1, int cat2, int cons1, int cons2)
    {   
        probA[num][cat1] = cons1;
        residueChoice[num][cat1] = cat2;
        probB[num][cat1] = cons2;
        
        int remainCons = 100 - cons1;        

        if(debugMode)
        {
            writeFileN("Molecule: " + molName + "\t Residue Number: " + num + "\t Category #1: " + cat1 + "\t conservation: " + cons1 + "\t Catgeory #2: " + cat2 + "\t conservation: " + cons2, debugMol);
        }
        
        /*
        Fill in remaining conservation information randomly.
        */
        int[] others = distributeEvents(molCat.length - 1, remainCons, cons1 - 1);
        
        int[] coord1 = randomWalk1D(molCat.length);
        int k = 0;
        
        for (int i = 0; i < coord1.length; i++)
        {
            int j = coord1[i];
            if(j != cat1)
            {
             probA[num][j] = others[k];
             k++;  
            }
        }
        
        filledA[num] = true;
    }

    
    public void fillProb()
    {   
        double fractionalSize = (double)(species);
        fractionalSize /= 100.0;
        
        for (int i = 0; i < probA.length; i++)
        {
            if(!filledA[i])
            {
                      int resIdent1 = (int)(Math.random() * molCat.length);                      
                      int resIdent2 = (int)(Math.random() * molCat[resIdent1].length);
                      int resCons1 = (int)(Math.random() * 100);
                      int resCons2 = (int)(Math.random() * 100);                      
                      setRes(i, resIdent1, resIdent2, resCons1, resCons2);
            }
        }
        
        for (int i = 0; i < probSumA.length; i++)
        {
            int runningTotal = 0;
         for (int j = 0; j < probSumA[i].length; j++)
         {
            runningTotal += probA[i][j];
            probSumA[i][j] = runningTotal;
         }
        }
        
        for (int i = 0; i < probA.length; i++)
        {
          int sumRow = 0;
            for (int j = 0; j < probA[i].length; j++)
            {
                obsA[i][j] = (int)(probA[i][j] * fractionalSize);
                sumRow += obsA[i][j];
            }
          
          if(sumRow < species)
          {
              int randomPos = (int)Math.round(Math.random() * (obsA[i].length-1));
              obsA[i][randomPos] += (species - sumRow);
          }
            
        }
        
    }
    
 /**
     * Method for choosing a residue when no Mutual Information is present - probabilities are global, not specified.
     * @param residueNo
     * @return 
     */
    public int chooseResidue(int residueNo)
    {           
                // Choose residue randomly, using summed probabilities
                int charChoice = chooseResidue(residueNo, probSumA[residueNo]);         
                return charChoice;
    }
    
    public int chooseTopResidue(int residueNo)
    {           
                //Choose residue based on highest observed count left - not on probability!
                    int topChoice = 0;

                    for (int i = 0; i < obsA[residueNo].length; i++)
                    { 
                        if(obsA[residueNo][i] > obsA[residueNo][topChoice])
                        {
                            topChoice = i;
                        };
                    }
                    
                return topChoice;
    }    
    
    /**
     * Method for choosing a residue when mutual information is present
     * @param residueNo
     * @param partnerRes
     * @param MIpairing
     * @return 
     */
    public int chooseResidue(int residueNo, int[] probSumT)
    {
                
                int charChoice = randomIdentity(probSumT);
                
                // Check to make sure that maximum number of observed counts has not been exceeded
                // If residue choosen by the probability the 2nd time is no good, select the them in a random walk order until it works
                int[] randoWalk = randomWalk1D(obsA[residueNo].length);
                
                if(obsA[residueNo][charChoice] <= 0)
                {
                    for (int i = 0; i < randoWalk.length; i++)
                    {
                        if(obsA[residueNo][randoWalk[i]] > 0)
                        {
                            charChoice = randoWalk[i];
                            break;
                        }
                    }
                }
                
                if(debugMode)
                {
                    writeFileN("\t OBS-A: " + printArr(obsA[residueNo]), debugMol);
                    writeFileN("\t PbS-T: " + printArr(probSumT), debugMol);
                    writeFileN("\t choice: " + charChoice, debugMol);
                }
                
                //Subtract from Observed Count
                obsA[residueNo][charChoice]--;
                
                return charChoice;
    }
    
    public char convertResidue(int charChoice1, int charChoice2)
    {
        return molCat[charChoice1][charChoice2];
    }
    
    public int getInGroupChoice(int res, int cat)
    {
        double rando = Math.random() * 100;

        int groupSize = molCat[cat].length;      
        if(probB[res][cat] >= rando)
        {
        return residueChoice[res][cat];            
        }
        else if(groupSize < 2)
        {
        return residueChoice[res][cat];               
        }
        else
        {   
            int[] remainingChoices = new int[groupSize -1];
            int j = 0;
            for (int i = 0; i < remainingChoices.length; i++)
            {
                if (j == residueChoice[res][cat])
                {                    
                    j++;
                }
                remainingChoices[i] = j;
                j++;
            }          
            double rando2 = Math.random() * remainingChoices.length;
            int otherChoice = (int)rando2;
            return remainingChoices[otherChoice];
                        
/*
int otherChoice = 0;
int iterations = 0;

do{
double randoChoice = Math.random() * groupSize;
otherChoice = (int)randoChoice;
iterations ++;         
System.out.print("\t randoChoice: " + randoChoice + "\t otherChoice" + otherChoice + "\t iterations" + iterations);            
}
while((otherChoice != residueChoice[res][cat]) && (iterations < groupSize * 3));

*/
        }
    }
        
    
    /*
    Choose a character identify based upon the probability distribution
    */
    public int randomIdentity(int[] probSums)
    {  

        int random100 = (int)Math.round(Math.random() * 100); 
     
            for (int k = 0; k < probSums.length; k++)
                {
                    if(random100 <= probSums[k])
                    {
                       return k;        
                    }
                }
        return 0;
         
    }    

    /*
    Debugging Block
    */
                            public void printArr(String[] arr){
                                for(int i = 0; i < arr.length; i++)
                                {
                                    System.out.print(arr[i] + "\t");
                                };
                                System.out.println("*");
                            }

                            public void printArr(char[] arr){
                                for(int i = 0; i < arr.length; i++)
                                {
                                    System.out.print(arr[i]);
                                    System.out.print("\t"); 
                                };
                                System.out.println("*");
                            }

                            public String printArr(int[] arr){
                                String res = "";
                                int sum1 = 0;
                                for(int i = 0; i < arr.length; i++)
                                {
                                    sum1 += arr[i];
                                    res += arr[i] + "\t";
                                };
                                res += "* ";
                                res += sum1;
                                return res;
                            }    

                            public void printArr(char[][] arr){
                                for(int i = 0; i < arr.length; i++)
                                 {                                     
                                    for(int j = 0; j < arr[i].length; j++)
                                     {
                                    System.out.print(arr[i][j]);
                                    System.out.print("\t");               
                                     }
                                    System.out.println("");            
                                 };
                                System.out.println("*");
                            }

                            public String printArr(int[][] arr){
                                String res = "";
                                int[] sum1 = new int[arr[0].length]; 
                                for(int i = 0; i < arr.length; i++)
                                 {
                                    int sum2 = 0;                                     
                                    res += printArr(arr[i]);                                   
                                    res += System.lineSeparator();
                                 };
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
    
    public int[][] cloneArray(int[][] array)
    {
        int[][] res = array;
        return res;
    }
    

    public static int[] distributeEvents(int size, int total, int absMax)
    {   
        int evenDis = (int)Math.round(total/size);
        
        if(absMax < evenDis)
        {
            absMax = evenDis;
        }
        
        int[] resArr = new int[size];
        int boxLeft = resArr.length - 1;
        int min = 0;
        
        for (int i = 0; i < resArr.length; i++)
        {
            int maximum = (int)Math.min(absMax, total);
            int excess = boxLeft * maximum;
            
            if (excess > total)
            {
              min = 0;  
            }
            else
            {
              min = total - excess;
            }
            
            int temp = (int)(Math.random() * maximum + min);
            
            if (i + 1 == resArr.length)
            {
                temp = total;
            }
            
            resArr[i] += temp;
            total -= temp;
        }
        
        return resArr;
    }

	private void writeFileN(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    content += "" + System.lineSeparator();
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(aces2.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            } 
    
}
