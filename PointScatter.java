import java.io.*;
import java.util.*;
import java.text.*;

public class PointScatter extends FitnessFunction{
/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public PointScatter(){
		name = "PointScatter Problem";
	}

/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

//  COMPUTE A CHROMOSOME'S RAW FITNESS *************************************

	public void doRawFitness(Chromo X){
        X.rawFitness = 0;
		double minDistance = 100;
		double[][] chromosome = new double[Parameters.numGenes][2];

		// Convert chromosome string to 2D arrray of polar coordinates
		switch (Parameters.chromosomeType){
			case 1: //  Decode binary string to polar values
				for (int i=0; i<Parameters.numGenes; i++){
					String gene = X.chromo.substring(i*Parameters.geneSize, (i+1)*Parameters.geneSize);
					String radius = gene.substring(0, 4);
					String angle = gene.substring(4, 8);
					chromosome[i][0] = Double.valueOf(Integer.parseInt(radius, 2)*0.0625);
					double degrees = Double.valueOf(Integer.parseInt(angle, 2)*22.5);
					double radians = Math.toRadians(degrees);
					chromosome[i][1] = radians;
				}
				break;
			case 2: //  Decode gray string to polar values
				//Transform gray to binary
				String binChromo = Chromo.graytoBinary(X.chromo);
				//Do operations in binary
				for (int i=0; i<Parameters.numGenes; i++){
					String gene = binChromo.substring(i*Parameters.geneSize, (i+1)*Parameters.geneSize);
					String radius = gene.substring(0, 4);
					String angle = gene.substring(4, 8);
					chromosome[i][0] = Double.valueOf(Integer.parseInt(radius, 2)*0.0625);
					double degrees = Double.valueOf(Integer.parseInt(angle, 2)*22.5);
					double radians = Math.toRadians(degrees);
					chromosome[i][1] = radians;
				}
				break;

			case 3: // Decode value string to polar values
				String[] c = X.chromo.split("/");
				for (int i=0; i<Parameters.numGenes; i++){
					String[] gene = c[i].split(",");
					chromosome[i][0] = Double.parseDouble(gene[0]);
					chromosome[i][1] = Double.parseDouble(gene[1]);	
				}
				break;
		}

		// Calculate fitness
		for (int a=0; a<Parameters.numGenes; a++){
			for (int b=0; b<Parameters.numGenes; b++){
				if (a != b && a > b){
					// Calculate the distance between two points, save the max. (Distance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta1 - theta2))
					minDistance = Math.min(minDistance, Math.sqrt(Math.pow(chromosome[a][0], 2) + Math.pow(chromosome[b][0], 2) - (2 * chromosome[a][0] * chromosome[b][0] * Math.cos(chromosome[b][1] - chromosome[a][1]))));
				}
			}
		}
		X.rawFitness = minDistance;
	}

//  PRINT OUT AN INDIVIDUAL GENE TO THE SUMMARY FILE *********************************

	public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{

		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getGeneAlpha(i),11,output);
		}
		output.write("   RawFitness");
		output.write("\n        ");
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getPosIntGeneValue(i),11,output);
		}
		Hwrite.right((int) X.rawFitness,13,output);
		output.write("\n\n");
		return;
	}


}