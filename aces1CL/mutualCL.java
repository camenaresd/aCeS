
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
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
public class mutualCL{
    /*
    The array of possibilities to be returned
    */
    public int[][] obsXY;
    
    /*
    Number assignment for first molecule 0-4. i = 0 is first molecule, i = 1 is 2nd
    */
    public int[] molNum = new int[2];
    
    public int[] resNum = new int[2];

    public File debugMI;
    
    public boolean debugMode;
    
    public double finalMI;
    
    public double targetMI;
    
    public int iterationsProcessed;
    
    /**
     * Method for creation of a mutual information object
     * @param sA
     * @param molA
     * @param resA
     * @param molB
     * @param resB
     * @param pctMI
     * @param sizeXY 
     */
    public mutualCL(int[] sA, moleculeCL molA, int resA, moleculeCL molB, int resB, int pctMI, int sizeXY)
    {
        
        molNum[0] = molA.getName();
        molNum[1] = molB.getName();
        
        resNum[0] += resA;
        resNum[1] += resB;
        
        int iterations = sA[0];
        double thresholdMIfit = (double)(sA[1]);
        thresholdMIfit /= 100.0;        
        double multiplierStDev = (double)(sA[2]);
        
        int pX1prob = molA.getCatSize();
        int pY1prob = molB.getCatSize();
        int varyX1 = molA.getResCon(resA);
        int varyY1 = molB.getResCon(resB);
        
        int[] oX1 = molA.getObs(resA);
        int[] oY1 = molB.getObs(resB);
 
        double[] pX1 = obsFreq(oX1);
        double[] pY1 = obsFreq(oY1);        
        
        double percentageMI = (double)(pctMI);
        percentageMI /= 100.0;
        
        /**
         * Create a 2D array reflecting pX and pY values - no mutual information
         */
        double[][] pXY1 = new double[pX1.length][pY1.length];
        
        for (int i = 0; i < pX1.length; i++){
            for (int j = 0; j < pY1.length; j++){
                pXY1[i][j] = pX1[i] * pY1[j];
            }
        }

        double noMI = calculateMI(pXY1);  

        int[][] oXYmax = simulateMI(oX1, oY1, sizeXY, 20000, 2000, iterations, 1); 
        double maxMI = calculateMI(oXYmax);
        
        // System.out.println("resA: " + resA + ", resB: " + resB + ", MaxMI: " + maxMI);

        // printArray(oXYmax);        
        
        // System.out.println("stDev" + "\t" + "target" + "\t" + "final" + "\t" + "devSum" + "\t" + "ups" + "\t" + "delta T");        
        
        targetMI = percentageMI * maxMI;        
                
        obsXY = simulateMI(oX1, oY1, sizeXY, multiplierStDev*percentageMI, targetMI, iterations, thresholdMIfit);  

        int devDis2 = arrayDeviation(obsXY, oX1, oY1);
        finalMI = calculateMI(obsXY);
        double diffMI = (Math.abs(targetMI - finalMI))*100/targetMI;
        
    }
    
    public int[][] getMI()
    {
        return obsXY;
    }
    
    public int[] getNames()
    {
        return molNum;
    }    

    public int[] getNums()
    {
        return resNum;
    }    

    public void debugSet(File filename)
    {
        debugMI = filename;
        debugMode = true;
        writeFileN("MI pair: " + molNum[0] + "-" + resNum[0] + " & " + molNum[1] + "-" + resNum[1] + "\t target: " + Math.round(targetMI*1000.0)/1000.0 + "\t final: " + Math.round(finalMI*1000.0)/1000.0 + "\t iterations: " + iterationsProcessed, debugMI);
// Remove just for now..     writeFileN(printArray(obsXY), debugMI);
    }
    
    /*
    THE BELOW CODE IS COPIED FROM A SANDBOX - Last updated 190126
    */
    
    /**
     * Converts an observed count into a frequency
     * @param arr
     * @return 
     */
    public static double[] obsFreq(int[] arr){

        double[] pZ = new double[arr.length];
        int sum = 0;
        
        for (int i = 0; i < arr.length; i++){
            sum += arr[i];
        }
            double sumD = (double)sum; 
        
        for (int i = 0; i < arr.length; i++){
            double arrD = (double)arr[i];

            pZ[i] = arrD/sumD;
        }
        
                
        return pZ;
    }
    
    /**
     * Converts an observed count into a frequency
     * @param arr
     * @return 
     */
    public static double[][] obsFreq(int[][] arr){

        double[][] pZ = new double[arr.length][arr[0].length];
        double sum = 0.0;
        
        for (int i = 0; i < arr.length; i++)
            {
             for (int j = 0; j < arr.length; j++)
             {
              sum += arr[i][j];              
             }

        }
        
        for (int i = 0; i < arr.length; i++){
            for (int j = 0; j < arr[i].length; j++)
                {
                 double arrD = (double)arr[i][j];
                 pZ[i][j] = arrD/sum;
                }

        }
        
                
        return pZ;
    }    

    /**
     * A function to create a duplicate of a 1 dimensional double array;
     * @param array
     * @return An exact replica array;
     */
    public static int[] cloneArray(int[] array){
        
        int[] array2 = new int[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            array2[i] = array[i] + 0;
        }
        
        return array2;
    }
    
        /**
     * A function to create a duplicate of a 1 dimensional double array;
     * @param array
     * @return An exact replica array;
     */
    public static int[][] cloneArray(int[][] array){
        
        int[][] array2 = new int[array.length][array[0].length];
        
        for (int i = 0; i < array.length; i++)
        {
            for (int j = 0; j < array[i].length; j++)
             {
             array2[i][j] = array[i][j] + 0;
             }
        }
        
        return array2;
    }
    
    
    /**
     * A function to create a duplicate of a 1 dimensional double array;
     * @param array
     * @return An exact replica array;
     */
    public static double[] cloneArray(double[] array){
        
        double[] array2 = new double[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            array2[i] = array[i] + 0.0;
        }
        
        return array2;
    }
    
    public static double[][] cloneArray(double[][] array){
        
        double[][] array2 = new double[array.length][array[0].length];
        
        for (int i = 0; i < array.length; i++)
        {
            for (int j = 0; j < array[i].length; j++){
             array2[i][j] = array[i][j] + 0.0;                
            }
        }
        
        return array2;
    }     
    
    public static int[] initArray(int[] array, int dummy){
        for(int i = 0; i < array.length; i++){
            array[i] = dummy;
        }
        return array;
    }
    
    /**
     * Create a distribution matrix of occurrences given a target amount of covariance
     * @param pX1 The probability distribution for identities in the first or X residue
     * @param pY1 The probability distribution for identities in the second or Y residue
     * @param stdevMI A setpoint for MI calculations
     * @param targetMI The target amount of MI to obtain
     * @param iterations The number of times this calculation will be run
     * @param alpha The ratio between target accuracy and distribution accuracy that is desired
     * @return 
     */
    public int[][] simulateMI(int[] oX1, int[] oY1, int setSize, double stdevMI, double targetMI, int iterations, double fit1){

        int[][] oXY = new int[oX1.length][oY1.length];


            // System.out.println("New MI Simulation");
            // System.out.println(printArray(oX1));
            // System.out.println(printArray(oY1));            

        
        int upgradeNumber = 0; 
        iterationsProcessed = 0;

        /*
          initialize the no MI array
          This will be used as the average value from which the adjustments will be made (if targetMI is greater than 0)
        */
        int[][] zeroMI = new int[oX1.length][oY1.length];
        
        for (int m = 0; m < oX1.length; m++)
         {
           for (int n = 0; n < oY1.length; n++)
             {
               double pX1 = ((double)oX1[m])/((double)setSize);
               double pY1 = ((double)oY1[n])/((double)setSize);
               
               double pXY1 = pX1*pY1;
               int oXY1 = (int)(pXY1*setSize);
               
               // System.out.println("pX1: " + pX1 + ", pY1: " + pY1 + ", pXY1: " + pXY1 + ", setSize: " + setSize + ", oXY1: " + oXY1);               
               zeroMI[m][n] = oXY1;
              }
          };
        
        // System.out.println("Zero MI ini:");
        // System.out.println(printArray(zeroMI));
          
        fillGaps(zeroMI, setSize, oX1, oY1, true);  

        // System.out.println("Zero MI, fillGaps:");
        // System.out.println(printArray(zeroMI));        
        
        oXY = cloneArray(zeroMI);
        
        // System.out.println("Zero MI oXY, fillGaps:");
        // System.out.println(printArray(oXY));
        
        /**
         * Iterate the creation of a mutual information grid in order to select from best version
         */
        for (int r = 0; r < iterations; r++){

                int[][] oXYt = new int[oX1.length][oY1.length];
            
                /**
                 * Create a randomized walk-through the 2D grid
                 */
                int[][] walk1 = randomWalk2D(oX1.length, oY1.length);

                /**
                 * Create duplicate arrays of probability values to determine the maximum allowed in each cell
                 */
                int[] oX2 = cloneArray(oX1);
                int[] oY2 = cloneArray(oY1);

                /*
                initialize the maximum array
                */
                int[][] maxXY = new int[oX1.length][oY1.length];

                for (int m = 0; m < oX1.length; m++)
                 {
                   for (int n = 0; n < oY1.length; n++)
                     {
                       maxXY[m][n] = Math.min(oX2[m], oY2[n]);
                      }
                  };

                /**
                 * Using the randomized walk-through, fill out values in the 2D array so that some mutual information is simulated
                 */            
                for (int k = 0; k < walk1.length; k++)
                {

                        int i = walk1[k][0];
                        int j = walk1[k][1];

                        int[] sumMaxX = sumRows(maxXY);
                        int[] sumMaxY = sumColumns(maxXY); 

                        /*
                        initialize the minimum array
                        */
                        int[][]minXY = new int[oX1.length][oY1.length];

                        for (int m = 0; m < oX1.length; m++)
                            {
                             for (int n = 0; n < oY1.length; n++)
                                {
                                 int diffX = oX1[m] - (sumMaxX[m] - maxXY[m][n]);
                                 int diffY = oY1[n] - (sumMaxY[n] - maxXY[m][n]);                   
                                 minXY[m][n] = Math.max(diffX, diffY);
                                 if(minXY[m][n] < 0)
                                 {
                                     minXY[m][n] = 0;
                                 }
                                }
                        };
                       int spaceRemaining = setSize - sumArr(oXYt);            

                       int currentMax = Math.min(oX2[i], oY2[j]);

                       if (currentMax > spaceRemaining)
                        {
                          currentMax = spaceRemaining + 0;  
                        }

                       // First putative adjustment figure
                       int temp1 = randomAdjustment(currentMax, minXY[i][j], oXY[i][j], stdevMI);
                       
                       // Prevent negative numbers from being used
                       temp1 = Math.max(temp1, 0);
                       
                       oXYt[i][j] += temp1;
                       maxXY[i][j] = temp1 + 0;

                       oX2[i] -= temp1;
                       oY2[j] -= temp1;

                 }      
                
                fillGaps(oXYt, setSize, oX1, oY1, false);                    
                
                double tempMI = calculateMI(oXYt);         
                double previousMI = calculateMI(oXY);

                double devMI = Math.abs(previousMI - targetMI);
                double devMIt = Math.abs(tempMI - targetMI);

                /* Array Deviation check - possibily unnecessary
                int devDIS = arrayDeviation(oXY, oX1, oY1);
                int devDISt = arrayDeviation(oXYt, oX1, oY1);
                boolean checkDIS = (devDISt < devDIS) || (devDISt == 0);         
                */
                
                iterationsProcessed++;

                /* Iteration replacement condition - choice for the below is specified by user thresholds and constants above */
                    if(devMIt < targetMI*fit1)
                     {
                        upgradeNumber++;
                        oXY = cloneArray(oXYt);
                        break;
                     }
                    else if (devMIt < devMI)
                     {
                        upgradeNumber++;
                        oXY = cloneArray(oXYt);
                     }
         }

       // System.out.println("New oXY array for MI");        
       // System.out.println(printArray(oXY));

       double finalMI = calculateMI(oXY);                
       
       return oXY;
            
    }
    
    /**
     * This method will take an array in which events have been distributed and make sure that the columns and rows have the correct totals
     * @param array
     * @param setSize
     * @param x1
     * @param y1
     * @return 
     */
    public static int[][] fillGaps(int[][] array, int setSize, int[] x1, int[] y1, boolean debug){

        int xL = x1.length;
        int yL = y1.length;
        int[][] coords1 = randomWalk2D(xL, yL);        
        
        for(int k = 0; k < xL*yL; k++){
            int i = coords1[k][0];
            int j = coords1[k][1];
            
            int[] xSum = sumRows(array);
            int[] ySum = sumColumns(array);
            int newSize = sumArr(array);
            
            int maxDiff = -1*Math.max((xSum[i]-x1[i]), (ySum[j]-y1[j]));
            
            int adjustmentValue = Math.min(maxDiff, (setSize-newSize));
                        
            array[i][j] += adjustmentValue;
            
            if(debug)
            {
            // System.out.println("k: " + k + ", i: " + i + ", j: " + j + ", xSum: " + printArray(xSum) + ", ySum: " + printArray(ySum) + ", setSize: " + setSize + ", newSize:" + newSize + ", maxDiff: " + maxDiff + ", adjustValue: " + adjustmentValue + ", arrayValue: " + array[i][j]);               
            }
            
            if((newSize+adjustmentValue) >= setSize)
            {
                k = xL*yL;
            }
            
        }
        
        return array;
    }    


    public static double calculateMI(int[][] oArray){
        
        double[][] array = obsFreq(oArray); 
        
        double res = calculateMI(array);
        
        return res;
        
    }
    
    /**
     * A method for calculating the mutual information of an array, given the starting probabilities of each row and column
     * @param array The array of pXY values
     * @param pX The array of pX values
     * @param pY The array of pY values
     */
    public static double calculateMI(double[][] array){
        
        double[] pX = sumRows(array);
        double[] pY = sumColumns(array);
        
        double sumMI = 0.0;
        
        double xyMI = 0.0;
        
        for (int i = 0; i < array.length; i++){
         for (int j = 0; j < array[i].length; j++){
             double pXpY = pX[i] * pY[j];
             
             double pXY = array[i][j] + 0.0;
             
             xyMI = pXY * Math.log(pXY/pXpY);                 
             
             if (Double.isNaN(xyMI)){
                 xyMI = 0.0;
             }
             
             sumMI += xyMI;
            }
        }
        
        return sumMI;
        
    }

    /**
     * Collapses a 2D array into columns
     * @param array
     * @return 
     */    
    public static int[] sumRows(int[][] array){
        int[] res = new int[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            int report = 0;
            for (int j = 0; j < array[i].length; j++)
            {
                report += array[i][j];
            }
            res[i] = report;
        }
        
        return res;
    }    
    
    /**
     * Collapses a 2D array into columns
     * @param array
     * @return 
     */    
    public static double[] sumRows(double[][] array){
        double[] res = new double[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            double report = 0.0;
            for (int j = 0; j < array[i].length; j++)
            {
                report += array[i][j];
            }
            res[i] = report;
        }
        
        return res;
    }

    /**
     * Collapses a 2D array into rows
     * @param array
     * @return 
     */
    public static int[] sumColumns(int[][] array){
        int[] res = new int[array[0].length];
        
        for (int i = 0; i < array.length; i++)
        {
            for (int j = 0; j < array[i].length; j++)
            {
                int report = 0;
                 for (int k = 0; k < array.length; k++)
                 {
                     report += array[k][j];
                 }
                res[j] = report; 
            }
        }
        
        return res;
    }   

    /**
     * Collapses a 2D array into rows
     * @param array
     * @return 
     */
    public static double[] sumColumns(double[][] array){
        double[] res = new double[array[0].length];
        
        for (int i = 0; i < array.length; i++)
        {
            for (int j = 0; j < array[i].length; j++)
            {
                double report = 0.0;
                 for (int k = 0; k < array.length; k++)
                 {
                     report += array[k][j];
                 }
                res[j] = report; 
            }
        }
        
        return res;
    }    
    
    /**
     * Collapses a 2D array into rows
     * @param array
     * @return 
     */
    public static int sumArr(int[][] array){
       int res = 0;
       
       for (int i = 0; i < array.length; i++)
       {
           for (int j = 0; j < array[i].length; j++)
           {
               res += array[i][j];
           }
       }
        
        return res;
    }      

    /**
     * Collapses a 2D array into rows
     * @param array
     * @return 
     */
    public static int sumArr(int[] array){
       int res = 0;
       
       for (int i = 0; i < array.length; i++)
       {
               res += array[i];
           
       }
        
        return res;
    }      
    
    /**
     * Gives the deviation between two arrays. Arrays must be the same size!
     * @param array1
     * @param array2
     * @return A number that totals the deviation in counts between two arrays
     */
    public static double arrayDeviation(double[] array1, double[]array2){
        double res = 0.0;
        
        for (int i = 0; i < array1.length; i++){
            res += Math.abs(array1[i] - array2[i]);
        }
        
        return res;
    }

    /**
     * Returns the absolute deviation between the values of a 2D array and the 1D row and column totals that are targeted
     * @param array2D
     * @param arrayX
     * @param arrayY
     * @return 
     */
    public static int arrayDeviation(int[][] array2D, int[] arrayX, int[] arrayY){
        int res = 0;
        
        int[] sumX = sumRows(array2D);
        int[] sumY = sumColumns(array2D);
        
        for (int i = 0; i < sumX.length; i++)
         {
          res += Math.abs(sumX[i] - arrayX[i]);  
         }

        for (int i = 0; i < sumY.length; i++)
         {
          res += Math.abs(sumY[i] - arrayY[i]);  
         }
        
        return res;
    }    
    
    /**
     * Randomly generate a number within certain bounds and with a particular distribution so as to make adjustments to another distribution of events
     * @param upperLimit The maximum value that is allowed to be returned
     * @param lowerLimit The minimum value that is allowed to be returned
     * @param pXY The mean for the distribution
     * @param stDev The standard deviation for the distribution (how tight is it around the mean)
     * @return 
     */
    public static int randomAdjustment(int upperLimit, int lowerLimit, int mean, double stDev)
    {

        int adjustment = 0;
        
        if (lowerLimit > upperLimit)
         {
            adjustment = upperLimit;
            return adjustment;
         }
        
        Random r = new Random();
        double randoGauss = r.nextGaussian();
        double adjustmentD = (randoGauss * stDev) + mean;
        
        adjustment = (int)adjustmentD;
        
        if (adjustment > upperLimit)
        {
            adjustment = upperLimit;
        }
        
        if (adjustment < lowerLimit)
        {
            adjustment = lowerLimit;
        }
                
        return adjustment;
    }
    
    /**
     * Randomly generate a number within certain bounds and with a particular distribution so as to make adjustments to another distribution of events
     * @param upperLimit The maximum value that is allowed to be returned
     * @param lowerLimit The minimum value that is allowed to be returned
     * @param pXY The mean for the distribution
     * @param stDev The standard deviation for the distribution (how tight is it around the mean)
     * @return 
     */
    public static double randomAdjustment(double upperLimit, double lowerLimit, double mean, double stDev)
    {
        double adjustment = 0.0;
        
        if (lowerLimit > upperLimit)
        {
            adjustment = upperLimit;
            return adjustment;
        }
        
        Random r = new Random();
        double randoGauss = r.nextGaussian();
        adjustment = (randoGauss * stDev) + mean;
        
        if (adjustment > upperLimit)
        {
            adjustment = upperLimit;
        }
        else if (adjustment < lowerLimit)
        {
            adjustment = lowerLimit;
        };
                
        return adjustment;
    }
 
    public static int[] randomWalk1D (int size){
        int[] coords1 = new int[size];
        
        /**
         * Initialize two identical arrays, apply the same random number to each cell
         */
        double[] rando1 = new double[size];
        double[] rando2 = new double[size];
        
        for (int i = 0; i < size; i++){
            double randoNum = Math.random();
            rando1[i] = randoNum;
            rando2[i] = randoNum;
        }
        
        // Sort one of the random array clones
        Arrays.sort(rando2);
        
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                if(rando2[i] == rando1[j]){
                    coords1[i] = j;
                }
            }
        }
       
        return coords1;
    }
    
    /**
     * This method will create a random walk through a 2D array
     * @param xSize The size of the 1st dimension of the array
     * @param ySize The size of the 2nd dimension of the array
     * @return An array which lists coordinates (x, y) in a randomized order
     */
    public static int[][] randomWalk2D (int xSize, int ySize){
                int[][] coords1 = new int[xSize * ySize][2];
        
        /**
         * A counter for a position among the 2D array
         */
        int z = 0;
        
        for (int i = 0; i < xSize; i++){
            for (int j = 0; j < ySize; j++){
                
                coords1[z][0] = i;
                coords1[z][1] = j;
                
                z++;
            }
        }
        
        /**
         * Initialize two identical arrays, apply the same random number to each cell
         */
        double[] rando1 = new double[xSize * ySize];
        double[] rando2 = new double[xSize * ySize];
        
        for (int i = 0; i < rando1.length; i++){
            double randoNum = Math.random();
            rando1[i] = randoNum;
            rando2[i] = randoNum;
        }
        
        // Sort one of the random array clones
        Arrays.sort(rando2);
        
        // Create a new array to translate random number position
        int[] rando3 = new int[xSize * ySize];
        
        for (int i = 0; i < rando3.length; i++){
            for (int j = 0; j < rando1.length; j++){
                if (rando1[j] == rando2[i])
                {
                    rando3[i] = j;
                }
            }
        }

        int[][] coords2 = new int[xSize * ySize][2];        
        int[][] coords3 = new int[xSize][ySize];
        
        z = 0;
        
        /**
         * Create a 2D version of the coordinate array, just to report the walk procedure back out (if necessary)
         */
        for (int i = 0; i < xSize; i++){
            for (int j = 0; j < ySize; j++){
                
                int iPos = coords1[rando3[z]][0];
                coords2[z][0] = iPos ;
                int jPos = coords1[rando3[z]][1];
                coords2[z][1] = jPos;
                
                coords3[iPos][jPos] = z;
                
                z++;
            }
        } 
        
        return coords2;
    }
    
    public static double doubleCheck (double input, double upperBound, double lowerBound)
    {
                double res = 0.0;
                
                res = input;
                
                if(input > upperBound)
                { 
                    res = upperBound;
                }
                else if(input < lowerBound)
                {
                    res = lowerBound;
                };
                
                return res;
    }

    public static String printArray(int[] aR)
    {
        String res = "";
        int runningTot = 0;
        for (int i = 0; i < aR.length; i++)
        {
            runningTot += aR[i];
            res += aR[i];
            res += "\t";
        }
        res += "* " + runningTot;
        return res;
    }

    public String printArray(int[][] aR)
    {
        String res = "";
        for (int i = 0; i < aR.length; i++)
        {
            res += printArray(aR[i]) + System.lineSeparator();
        }
        return res;
    }
    
	private void writeFileN(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    content += "" + System.lineSeparator();
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(aces1cl.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            }
    /*
        Build a method for testing the enginge (how well does it match the MI)
    public void testEngine(int trials, int cons1, int cons2, int testSize, double testMSD, double testIts, double testThres)
    {
        for (int i = 0; i < trials; i++)
        {
        
        int[][] testMI = simulateMI(tX1, tY1, testSize, targetMI, testIts, testThres); 
        
        }
    }
   */
}
