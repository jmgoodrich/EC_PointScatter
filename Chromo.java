/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class Chromo
{
/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	public String chromo;
	public double rawFitness;
	public double sclFitness;
	public double proFitness;

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	private static double randnum;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public Chromo(){

	switch (Parameters.chromosomeType){

		case 1:     //  Set gene values to a randum sequence of 1's and 0's
			char geneBit;
			chromo = "";
			for (int i=0; i<Parameters.numGenes; i++){
				for (int j=0; j<Parameters.geneSize; j++){
					randnum = Search.r.nextDouble();
					if (randnum > 0.5) geneBit = '0';
					else geneBit = '1';
					this.chromo = chromo + geneBit;
				}
			}
			break;

		case 2: 	//genetate binary number and then convert to gray
			char currBit;
			String binChromo = "";
			//Generate binary string
			for (int i=0; i<Parameters.numGenes; i++){
				for (int j=0; j<Parameters.geneSize; j++){
					randnum = Search.r.nextDouble();
					if (randnum > 0.5) currBit = '0';
					else currBit = '1';
					binChromo = binChromo + currBit;
				}
			}
			//Transform to gray string
			//System.out.println("Binary representation: " + binChromo + System.lineSeparator());

			this.chromo = binarytoGray(binChromo);
			//System.out.println("Gray representation: " + this.chromo + System.lineSeparator());

			break;


		case 3: 	//  Set gene values to random numbers for radius and angle.
			chromo = "";
			for (int i=0; i<Parameters.numGenes; i++){
				double radius = Search.r.nextDouble();
				double angle = Search.r.nextDouble() * 2 * Math.PI;
				this.chromo = chromo + String.format("%.4f", radius) + "," + String.format("%.4f", angle) + "/";
			}
			break;

		default:
			System.out.println("ERROR - No chromosome type selected");
		}	

		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized

	}


/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

	//  Get Alpha Represenation of a Gene **************************************

	public String getGeneAlpha(int geneID){
		int start = geneID * Parameters.geneSize;
		int end = (geneID+1) * Parameters.geneSize;
		String geneAlpha = this.chromo.substring(start, end);
		return (geneAlpha);
	}

	//  Get Integer Value of a Gene (Positive or Negative, 2's Compliment) ****

	public int getIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneSign;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=1; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		geneSign = geneAlpha.charAt(0);
		if (geneSign == '1') geneValue = geneValue - (int)Math.pow(2.0, Parameters.geneSize-1);
		return (geneValue);
	}

	//  Get Integer Value of a Gene (Positive only) ****************************

	public int getPosIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=0; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		return (geneValue);
	}

	//  Mutate a Chromosome Based on Mutation Type *****************************

	public void doMutation(){

		String mutChromo = "";
		char x;

		switch (Parameters.mutationType){

		case 1:     //  Bit flip mutation
			switch(Parameters.chromosomeType){
				case 1: //Binary Mutation
					for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
						x = this.chromo.charAt(j);
						randnum = Search.r.nextDouble();
						if (randnum < Parameters.mutationRate){
							if (x == '1') x = '0';
							else x = '1';
						}
						mutChromo = mutChromo + x;
					}
					this.chromo = mutChromo;
					break;
				case 2: 	//Gray scale mutation
					String binString = graytoBinary(this.chromo);
		
					for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
						x = binString.charAt(j);
						randnum = Search.r.nextDouble();
						if (randnum < Parameters.mutationRate){
							if (x == '1') x = '0';
							else x = '1';
						}
						mutChromo = mutChromo + x;
					}
					this.chromo = binarytoGray(mutChromo);	
					break;			
				default:
					System.out.println("ERROR - Mutation type and Representation dont match");
			}
			break;
		case 3: 	// Value mutation
			double[][] chromosome = new double[Parameters.numGenes][2];
			for (int j=0; j<(Parameters.geneSize); j++){
				String[] c = this.chromo.split("/");
				for (int i=0; i<Parameters.numGenes; i++){
					String[] gene = c[i].split(",");
					chromosome[i][0] = Double.parseDouble(gene[0]);
					chromosome[i][1] = Double.parseDouble(gene[1]);	
				}
			}
			/*For value mutation we are currently replacing the value of the chromosome with a new random number which is a big mutation
			 * for binary representation we shift only one bit. perhaps we can scale the change and increase or decrease by a differential?
			 */
			for (int k=0; k<Parameters.numGenes; k++){
				randnum = Search.r.nextDouble();
				if (randnum < Parameters.mutationRate){
					double anotherRand = Search.r.nextDouble();
					double delta = chromosome[k][0]*0.1;
					if(anotherRand > 0.5){
						chromosome[k][0] = chromosome[k][0] + delta;
					}
					else{
						chromosome[k][0] = chromosome[k][0] - delta;
					}
				}
				randnum = Search.r.nextDouble();
				if (randnum < Parameters.mutationRate){
					double anotherRand = Search.r.nextDouble();
					double delta = chromosome[k][1]*0.1;
					if(anotherRand > 0.5){
						chromosome[k][1] = chromosome[k][1] + delta;
					}
					else{
						chromosome[k][1] = chromosome[k][1] - delta;
					}
				}
				mutChromo = mutChromo + String.format("%.4f", chromosome[k][0]) + "," + String.format("%.4f", chromosome[k][1]) + "/";
			}
			break;

		default:
			System.out.println("ERROR - No mutation method selected");
		}
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	//  Select a parent for crossover ******************************************

	public static int selectParent(){

		double rWheel = 0;
		int j = 0;
		int k = 0;
		double best = 0;
		int victor = 0;
		int tournSize = 3;

		switch (Parameters.selectType){

		case 1:     // Proportional Selection
			randnum = Search.r.nextDouble();
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				if (randnum < rWheel) return(j);
			}
			break;

		case 2:     //  Tournament Selection
			
			for (j=0; j<Parameters.popSize; j++){
				for (k=0; k<tournSize; k++){
					randnum = Search.r.nextDouble();
					int contestant = (int) (randnum * Parameters.popSize);
					if (best < Search.member[contestant].proFitness){
						best = Search.member[contestant].proFitness;
						victor = contestant;
					}
				}
				return(victor);
			}

		case 3:     // Random Selection
			randnum = Search.r.nextDouble();
			j = (int) (randnum * Parameters.popSize);
			return(j);

		case 4:		// Rank Selection



		default:
			System.out.println("ERROR - No selection method selected");
		}
	return(-1);
	}

	//  Produce a new child from two parents  **********************************

	public static void mateParents(int pnum1, int pnum2, Chromo parent1, Chromo parent2, Chromo child1, Chromo child2){

		int xoverPoint1;
		int xoverPoint2;

		switch (Parameters.xoverType){

		case 1:     //  Single Point Crossover
			//  Select crossover operator based on encoding type
			switch(Parameters.chromosomeType){
				case 1:/*Binary encoding */
					xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));

					//  Create child chromosome from parental material
					child1.chromo = parent1.chromo.substring(0,xoverPoint1) + parent2.chromo.substring(xoverPoint1);
					child2.chromo = parent2.chromo.substring(0,xoverPoint1) + parent1.chromo.substring(xoverPoint1);
					break;
				case 2: /*Gray scale encoding */
					xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));
					//Convert parent gray representation to binary
					String p1Binchromo = graytoBinary(parent1.chromo);
					String p2Binchromo = graytoBinary(parent2.chromo);

					//  Create child chromosome from parental material and transform back to gray scale
					child1.chromo = binarytoGray(p1Binchromo.substring(0,xoverPoint1) + p2Binchromo.substring(xoverPoint1));
					child2.chromo = binarytoGray(p2Binchromo.substring(0,xoverPoint1) + p1Binchromo.substring(xoverPoint1));
				break;
				case 3: /*Value encoding crossover*/				
					//System.out.print("Printing Parent 1: "+ parent1.chromo + System.lineSeparator());
					//System.out.print("Printing Parent 2: "+ parent2.chromo + System.lineSeparator());

						String[] c1 = parent1.chromo.split("/");
						String[] c2 = parent2.chromo.split("/");
						String c1Chromo = "";
						String c2Chromo = "";
						for (int i=0; i<Parameters.numGenes; i++){
							String[] p1genes = c1[i].split(",");
							double p1radius = Double.parseDouble(p1genes[0]);
							double p1angle = Double.parseDouble(p1genes[1]);
							String[] p2genes = c2[i].split(",");
							double p2radius = Double.parseDouble(p2genes[0]);
							double p2angle = Double.parseDouble(p2genes[1]);
							double random = Search.r.nextDouble();
							double child1Rad = (1-random)*p1radius + (random)*p2radius;
							double child2Rad = (1-random)*p2radius + (random)*p1radius;
							double child1Angle = (1-random)*p1angle + (random)*p2angle;
							double child2Angle = (1-random)*p2angle + (random)*p1angle;
							c1Chromo = c1Chromo + String.format("%.4f", child1Rad) + "," + String.format("%.4f", child1Angle) + "/";
							c2Chromo= c2Chromo + String.format("%.4f", child2Rad) + "," + String.format("%.4f", child2Angle) + "/";
						}

					//System.out.print("Printing Child 1: "+ c1Chromo + System.lineSeparator());
					//System.out.print("Printing Child 2: "+ c2Chromo + System.lineSeparator());
					child1.chromo = c1Chromo;
					child2.chromo = c2Chromo;

					break;
			}
			break;


		case 2:     //  Two Point Crossover

		case 3:     //  Uniform Crossover

		case 4:		// Custom Crossover 

		default:
			System.out.println("ERROR - Bad crossover method selected");
		}

		//  Set fitness values back to zero
		child1.rawFitness = -1;   //  Fitness not yet evaluated
		child1.sclFitness = -1;   //  Fitness not yet scaled
		child1.proFitness = -1;   //  Fitness not yet proportionalized
		child2.rawFitness = -1;   //  Fitness not yet evaluated
		child2.sclFitness = -1;   //  Fitness not yet scaled
		child2.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Produce a new child from a single parent  ******************************

	public static void mateParents(int pnum, Chromo parent, Chromo child){

		//  Create child chromosome from parental material
		child.chromo = parent.chromo;

		//  Set fitness values back to zero
		child.rawFitness = -1;   //  Fitness not yet evaluated
		child.sclFitness = -1;   //  Fitness not yet scaled
		child.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Copy one chromosome to another  ***************************************

	public static void copyB2A (Chromo targetA, Chromo sourceB){

		targetA.chromo = sourceB.chromo;

		targetA.rawFitness = sourceB.rawFitness;
		targetA.sclFitness = sourceB.sclFitness;
		targetA.proFitness = sourceB.proFitness;
		return;
	}
	//HELPER FUNCTIONS TO CONVERT FROM BINARY TO GRAY ENCODING AND VICE VERSA
	//  X-OR bits   					***************************************
	public static char bit_xor(char a, char b){
        return (a == b) ? '0' : '1';
	}
	//  Flip Bits
	public static char bit_flip(char a){
        return (a == '0') ? '1' : '0';
	}
	//  Convert gray code to binary  	***************************************

	public static String graytoBinary(String gray)
    {
        String bin = "";

        bin += gray.charAt(0);
 
        for (int i = 1; i < gray.length(); i++) 
        {
            if (gray.charAt(i) == '0')
				bin += bin.charAt(i - 1);
            else
				bin += bit_flip(bin.charAt(i - 1));
        }
        return bin;
    }

	//  Convert binary to gray code   	***************************************
	public static String binarytoGray(String bin)
    {
        String gray = "";
        gray += bin.charAt(0);

        for (int i = 1; i < bin.length(); i++) 
        {
            gray += bit_xor(bin.charAt(i - 1),
                          bin.charAt(i));
        }
       return gray;
    }



}   // End of Chromo.java ******************************************************
