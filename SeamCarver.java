import edu.princeton.cs.algs4.Picture;
import java.awt.Color;
import edu.princeton.cs.algs4.StdOut; 

public class SeamCarver {
    private Picture picture; 
    
    // Create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException();
        }
        this.picture = picture; 
    }

    // Current picture
    public Picture picture() {
        return picture;
    }

    // Width of current picture
    public int width() {
        return picture.width();
    }

    // Height of current picture
    public int height() {
        return picture.height();
    }

    // Energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= width() || y < 0 || y >= height()) {
            throw new IllegalArgumentException("Coordinates out of bounds: (" + x + ", " + y + ")");
        }
        if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1) { 
            return 1000.0;
        }
        double deltaX = 0.0, deltaY = 0.0;

        // Calculate deltaX
        Color color1 = picture.get(x + 1, y);
        Color color2 = picture.get(x - 1, y);
        deltaX += Math.pow(color1.getRed() - color2.getRed(), 2);
        deltaX += Math.pow(color1.getGreen() - color2.getGreen(), 2);
        deltaX += Math.pow(color1.getBlue() - color2.getBlue(), 2);

        // Calculate deltaY
        color1 = picture.get(x, y + 1);
        color2 = picture.get(x, y - 1);
        deltaY += Math.pow(color1.getRed() - color2.getRed(), 2);
        deltaY += Math.pow(color1.getGreen() - color2.getGreen(), 2);
        deltaY += Math.pow(color1.getBlue() - color2.getBlue(), 2);

        return Math.sqrt(deltaX + deltaY);
    }

    public int[] findHorizontalSeam() {
        double[][] energyMatrix = new double[height()][width()];
        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                energyMatrix[y][x] = energy(x, y);
            }
        }
        double[][] transposedEnergyMatrix = transposeMatrix(energyMatrix);
        return findSeam(transposedEnergyMatrix, true);
    }
    // Sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        double[][] energyMatrix = new double[height()][width()];
        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                energyMatrix[y][x] = energy(x, y);
            }
        }
        return findSeam(energyMatrix, true);
    }
    private double[][] transposeMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[][] transposedMatrix = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }
    // General method to find seam
    private int[] findSeam(double[][] energyMatrix, boolean vertical) {
        //int width = energyMatrix[0].length;
        //int height = energyMatrix.length;

	int width = vertical ? energyMatrix[0].length : energyMatrix.length;
        int height = vertical ? energyMatrix.length : energyMatrix[0].length;
        double[][] distTo = new double[height][width];
        int[][] edgeTo = new int[height][width];

        // Initialize distance to infinity
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                distTo[y][x] = Double.POSITIVE_INFINITY;
            }
        }

        // Initialize top row
        for (int x = 0; x < width; x++) {
            distTo[0][x] = energyMatrix[0][x];
        }

        // Relax edges in topological order
        for (int y = 0; y < height - 1; y++) {
            for (int x = 0; x < width; x++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int nx = x + dx;
                    if (nx >= 0 && nx < width) {
                        if (distTo[y + 1][nx] > distTo[y][x] + energyMatrix[y + 1][nx]) {
                            distTo[y + 1][nx] = distTo[y][x] + energyMatrix[y + 1][nx];
                            edgeTo[y + 1][nx] = x;
                        }
                    }
                }
            }
        }

        // Find minimum energy in the bottom row
        double minEnergy = Double.POSITIVE_INFINITY;
        int minIndex = -1;
        for (int x = 0; x < width; x++) {
            if (distTo[height - 1][x] < minEnergy) {
                minEnergy = distTo[height - 1][x];
                minIndex = x;
            }
        }

        // Trace back to find the seam
        int[] seam = new int[height];
        for (int y = height - 1; y >= 0; y--) {
            seam[y] = minIndex;
            minIndex = edgeTo[y][minIndex];
        }

        return seam;
    }

    // Remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (!isValidSeam(seam, false)) {
            throw new IllegalArgumentException("Illegal seam!");
        }
        Picture seamedPicture = new Picture(width(), height() - 1);
        for (int col = 0; col < width(); col++) {
            int rowBias = 0;
            for (int row = 0; row < height() - 1; row++) {
                if (seam[col] == row) {
                    rowBias = 1;
                }
                seamedPicture.set(col, row, picture.get(col, row + rowBias));
            }
        }
        this.picture = seamedPicture;
    }

    // Remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (!isValidSeam(seam, true)) {
            throw new IllegalArgumentException("Illegal seam!");
        }
        Picture seamedPicture = new Picture(width() - 1, height());
        for (int row = 0; row < height(); row++) {
            int colBias = 0;
            for (int col = 0; col < width() - 1; col++) {
                if (seam[row] == col) {
                    colBias = 1;
                }
                seamedPicture.set(col, row, picture.get(col + colBias, row));
            }
        }
        this.picture = seamedPicture;
    }

    // Check if seam is valid
    private boolean isValidSeam(int[] seam, boolean vertical) {
        if (seam == null) {
            return false;
        }
        if (vertical) {
            if (seam.length != height()) {
                return false;
            }
            for (int row = 0; row < height() - 1; row++) {
                if (Math.abs(seam[row] - seam[row + 1]) > 1) {
                    return false;
                }
            }
        } else {
            if (seam.length != width()) {
                return false;
            }
            for (int col = 0; col < width() - 1; col++) {
                if (Math.abs(seam[col] - seam[col + 1]) > 1) {
                    return false;
                }
            }
        }
        return true;
    }
    public static void main(String[] args) {
        if (args.length != 1) {
            StdOut.println("Usage: java SeamCarver [image filename]");
            return;
        }

        Picture picture = new Picture(args[0]);
        SeamCarver seamCarver = new SeamCarver(picture);

        StdOut.println("Original picture:");
        picture.show();

        int[] verticalSeam = seamCarver.findVerticalSeam();
        StdOut.println("Vertical seam: " + java.util.Arrays.toString(verticalSeam));
        seamCarver.removeVerticalSeam(verticalSeam);
        Picture verticalSeamRemoved = seamCarver.picture();
        verticalSeamRemoved.show();
        verticalSeamRemoved.save("verticalSeamRemoved.png");

        int[] horizontalSeam = seamCarver.findHorizontalSeam();
        StdOut.println("Horizontal seam: " + java.util.Arrays.toString(horizontalSeam));
        seamCarver.removeHorizontalSeam(horizontalSeam);
        Picture horizontalSeamRemoved = seamCarver.picture();
        horizontalSeamRemoved.show();
        horizontalSeamRemoved.save("horizontalSeamRemoved.png");

        StdOut.println("Processing complete.");
    }
}

