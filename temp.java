import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Picture;
import java.util.Arrays;
import java.awt.Color;

public class SeamCarver {
   private Picture picture; 
   
   // create a seam carver object based on the given picture
   public SeamCarver(Picture picture) {
       if (picture == null){
            throw new IllegalArgumentException();
       }
       this.picture = picture; 
   }

   // current picture
   public Picture picture() {
   	return picture;
   }

   // width of current picture
   public int width() {
       return picture.width();
   }

   // height of current picture
   public int height() {
   	return picture.height();
   }

   // energy of pixel at column x and row y
   /*public double energy(int x, int y) {
   	if (x < 0 || x >= width() || y < 0 || y >= height()) {
	    throw new IllegalArgumentException();
	}
	if (x == 0 || x == picture.width()-1 || y == 0 || y == picture.height()-1) { 
		return (double) 1000.0;
	}
	double ret = 0.0;
	Color color1 = picture.get(x+1, y); 
	int r1 = color1.getRed(); 
	int g1 = color1.getGreen(); 
	int b1 = color1.getBlue(); 
	Color color2 = picture.get(x-1, y); 
	int r2 = color2.getRed(); 
        int g2 = color2.getGreen(); 
        int b2 = color2.getBlue(); 

	double deltaX = Math.pow(r1-r2, 2) + Math.pow(g1-g2, 2) + Math.pow(b1-b2, 2);

	color1 = picture.get(x, y+1);
        r1 = color1.getRed(); 
        g1 = color1.getGreen(); 
        b1 = color1.getBlue(); 
	color2 = picture.get(x, y-1);
        r2 = color2.getRed();    
        g2 = color2.getGreen();    
        b2 = color2.getBlue();     

        double deltaY = Math.pow(r1-r2, 2) + Math.pow(g1-g2, 2) + Math.pow(b1-b2, 2);

	ret = Math.sqrt(deltaX + deltaY);
	return ret; 
   }*/
   public double energy(int x, int y) {
        if (x < 0 || x >= width() || y < 0 || y >= height()) {
            throw new IllegalArgumentException("Coordinates out of bounds: (" + x + ", " + y + ")");
        }
        if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1) {
            return 1000.0;
        }

        double deltaX2 = deltaSquared(picture.get(x - 1, y), picture.get(x + 1, y));
        double deltaY2 = deltaSquared(picture.get(x, y - 1), picture.get(x, y + 1));

        return Math.sqrt(deltaX2 + deltaY2);
    }

    private double deltaSquared(Color c1, Color c2) {
        int dr = c1.getRed() - c2.getRed();
        int dg = c1.getGreen() - c2.getGreen();
        int db = c1.getBlue() - c2.getBlue();
        return dr * dr + dg * dg + db * db;
    }

   // sequence of indices for horizontal seam
   /*public int[] findHorizontalSeam() { 
   	int[] ret = new int[picture.width()]; 
        for (int i = 0; i < picture.width(); i++) {
            double min = energy(i, 0);
	    int index = 0;
            for (int j = 0; j < picture.height(); j++) {
                if (energy(j, i) < min) {
                    min = energy(j, i);
		    index = j;
                }
            }
            ret[i] = index; 
        }
	return ret;
   }

   // sequence of indices for vertical seam
   public int[] findVerticalSeam() {
   	int[] ret = new int[picture.height()]; 
	for (int i = 0; i < picture.height(); i++) { 
	    double min = energy(0, i);
	    int index = 0;
	    for (int j = 0; j < picture.width(); j++) { 
	        if (energy(j, i) < min) { 
		    min = energy(j, i);
	 	    index = j;
		}
	    }
	    ret[i] = index; 
	}
	return ret;
   }*/

  // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        double[][] energyMatrix = new double[width()][height()];
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                energyMatrix[x][y] = energy(x, y);
            }
        }
        return findSeam(energyMatrix, false);
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        double[][] energyMatrix = new double[height()][width()];
        for (int x = 0; x < height(); x++) {
            for (int y = 0; y < width(); y++) {
                energyMatrix[x][y] = energy(x, y);
            }
        }
        return findSeam(energyMatrix, true);
    }

    private int[] findSeam(double[][] energyMatrix, boolean isVertical) {
        int height = energyMatrix.length;
        int width = energyMatrix[0].length;
        double[][] distTo = new double[height][width];
        int[][] edgeTo = new int[height][width];

        // Initialize the first row of distTo with the energy values
        for (int x = 0; x < width; x++) {
            distTo[0][x] = energyMatrix[0][x];
        }

        // Initialize the rest of distTo with positive infinity
        for (int y = 1; y < height; y++) {
            for (int x = 0; x < width; x++) {
                distTo[y][x] = Double.POSITIVE_INFINITY;
            }
        }

        // Relax edges in topological order
        for (int y = 0; y < height - 1; y++) {
            for (int x = 0; x < width; x++) {
                if (x > 0) {
                    if (distTo[y][x] + energyMatrix[y + 1][x - 1] < distTo[y + 1][x - 1]) {
                        distTo[y + 1][x - 1] = distTo[y][x] + energyMatrix[y + 1][x - 1];
                        edgeTo[y + 1][x - 1] = x;
                    }
                }
                if (distTo[y][x] + energyMatrix[y + 1][x] < distTo[y + 1][x]) {
                    distTo[y + 1][x] = distTo[y][x] + energyMatrix[y + 1][x];
                    edgeTo[y + 1][x] = x;
                }
                if (x < width - 1) {
                    if (distTo[y][x] + energyMatrix[y + 1][x + 1] < distTo[y + 1][x + 1]) {
                        distTo[y + 1][x + 1] = distTo[y][x] + energyMatrix[y + 1][x + 1];
                        edgeTo[y + 1][x + 1] = x;
                    }
                }
            }
        }

        // Find the minimum energy path in the bottom row
        double minEnergy = Double.POSITIVE_INFINITY;
        int minIndex = 0;
        for (int x = 0; x < width; x++) {
            if (distTo[height - 1][x] < minEnergy) {
                minEnergy = distTo[height - 1][x];
                minIndex = x;
            }
        }

        // Reconstruct the seam path from bottom to top
        int[] seam = new int[height];
        int y = height - 1;
        for (int x = minIndex; y >= 0; y--) {
            seam[y] = x;
            x = edgeTo[y][x];
        }

        if (!isVertical) {
            // Transpose the seam for horizontal cases
            int[] transposedSeam = new int[seam.length];
            for (int i = 0; i < seam.length; i++) {
                transposedSeam[i] = seam[seam.length - 1 - i];
            }
            seam = transposedSeam;
        }

        return seam;
    }
  private boolean isValidSeam(int[] seam, boolean vertical){

        if (seam == null){
            return false;
        }

        if ((vertical && seam.length != height()) || (!vertical && seam.length != width())){
            return false;
        }

        for(int i : seam){
            if ((i < 0 ) || (vertical && i >= width()) || (!vertical && i>= height())){
                return false;
            }
        }
        for (int i = 0; i < seam.length - 1; i++){
            if (Math.abs(seam[i] - seam[i + 1]) > 1){
                return false;
            }
        }
        return true;
   }

   // remove horizontal seam from current picture
   public void removeHorizontalSeam(int[] seam) {
	if (!isValidSeam(seam, false)){
            throw new IllegalArgumentException("Illegal seam!");
        }
	
   	Picture seamedPicture = new Picture(width(), height() - 1);	
	for (int col = 0; col < width(); col++){
            int rowBias = 0;
            for (int row = 0; row < height() - 1; row++){
                if (seam[col] == row){
                    rowBias = 1;
                }
                seamedPicture.set(col, row, picture.get(col, row + rowBias));
            }
        }
	this.picture = seamedPicture;
   }

   // remove vertical seam from current picture
   public void removeVerticalSeam(int[] seam) {
	if (!isValidSeam(seam, true)){
            throw new IllegalArgumentException("Illegal seam!");
        }
	Picture seamedPicture = new Picture(width() - 1, height());
        for(int row = 0; row < height(); row++){
            int colBias = 0;
            for(int col = 0; col < width() - 1; col++){
                if (seam[row] == col){
                    colBias = 1;
                }
                seamedPicture.set(col, row, picture.get(col + colBias, row));
            }
        }
	this.picture = seamedPicture; 
   }

   //  unit testing (optional)
   public static void main(String[] args) {
   }
};

