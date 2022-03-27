import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

public class Main {
	static InputReader ir = new InputReader(
			new File("C:\\Users\\BIBIMONI\\eclipse-workspace\\AI prediction project\\src\\input.txt"));
	static double[] result;

	public static void main(String[] args) throws FileNotFoundException {

		int n = ir.getLength();
		String input = "";
		double[] y = new double[n]; // left input column
		double[] a1 = new double[n]; // right input column
		int[] a2 = new int[n]; // vector 1
		Arrays.fill(a2, 1);
		double[] x = new double[n]; // x is a vector contains a and b for the prediction
		input = ir.readFile().trim();
		splitArr(n, input, y, a1);
		double[][] A = new double[n][2]; // matrix A contains vector a1 and a2
		A = initializeMatrix(a1, a2);
		double[][] AT = transposeMatrix(A);
		double[][] AAT = multiplyMatrix(AT, A);
		double[] ATy = multiplyVector(y, AT);
		double[][] adj = invert(AAT);
		// result[0] is a //result[1] is b
		result = multiplyVector(ATy, adj);
		inputting();
	}

	public static void inputting() {
		Scanner input = new Scanner(System.in);
		System.out.println("nhap cot trai hay cot phai ? | nhap 2 : trai / 1 : phai");
		int in = input.nextInt();
		if (in == 1) {
			solveForL();
		} else if (in == 2) {
			solveForR();
		}

	}

	// find the y value follows the a and b
	public static void solveForL() {
		double a = result[0];
		double b = result[1];
		Scanner input = new Scanner(System.in);
		// input x
		double x = input.nextDouble();
		input.nextLine();
		double result = a * x + b;
		System.out.println(result);
	}

	// find the x value follows the a and b
	public static void solveForR() {
		double a = result[0];
		double b = result[1];
		Scanner input = new Scanner(System.in);
		// input y
		double y = input.nextDouble();
		double result = (y - b) / a;
		System.out.println(result);
	}

	/*
	 * public static int[][] convertArrToInt(double[][] a) { int[][] temp = new
	 * int[a.length][a[0].length]; for(int i = 0; i < a.length; i++) { for(int j =
	 * 0; j < a[0].length; j++) { temp[i][j] = (int) a[i][j]; } } return temp; }
	 */
	// split the output to 2 different array
	public static void splitArr(int length, String input, double[] y, double[] a1) {
		String[] temp = new String[length];
		temp = input.split(" ");
		int m = 0;
		int n = 0;
		// length*2 because we need to access 2 arrays
		for (int i = 0; i < length * 2; i++) {
			if ((i + 1) % 2 == 1) {
				y[m] = Integer.valueOf(temp[i]);
				m++;
			}
			if ((i + 1) % 2 == 0) {
				a1[n] = Integer.valueOf(temp[i]);
				n++;
			}
		}
	}
	//put a1 and a2 into 1 matrix
	public static double[][] initializeMatrix(double[] a1, int[] a2) {
		double[][] result = new double[ir.getLength()][2];
		for (int i = 0; i < ir.getLength(); i++) {
			result[i][0] = a1[i];
			result[i][1] = a2[i];
		}
		return result;
	}

	// multiply a matrix and a vector
	public static double[] multiplyVector(double[] y, double[][] aT) {
		int RY = y.length;
		int CA = aT[0].length;
		int RA = aT.length;
		double[] result = new double[RA];
		if (RY != CA) {
			System.out.println("multiply failed");
			return null;
		} else {
			for (int i = 0; i < RA; i++) {
				result[i] = 0;
				for (int k = 0; k < CA; k++) {
					result[i] += aT[i][k] * y[k];
				}
			}
			return result;
		}
	}

	// multiply 2 matrices
	public static double[][] multiplyMatrix(double[][] aT, double[][] a) {
		int RA = aT.length; // row A
		int CB = a[0].length; // column B
		int RB = a.length; // row B
		int CA = aT[0].length; // collumn A
		int n = ir.getLength();
		double[][] result = new double[RA][CB];
		// if column A not equals to B, return
		if (CA != RB) {
			System.out.println("multiply failed");
			return null;
		} else {
			// goes to each element of the result
			// go to each row of A first because the i element of the result is a row
			for (int i = 0; i < RA; i++) {
				for (int j = 0; j < CB; j++) {
					result[i][j] = 0;
					for (int k = 0; k < CA; k++) {
						// column of A is also row of B
						result[i][j] += aT[i][k] * a[k][j];
					}

				}
			}
			return result;
		}
	}
	//transpose a matrix
	public static double[][] transposeMatrix(double[][] a) {
		int n = ir.getLength();
		double[][] transpose = new double[2][n]; // 2 is a1 a2 but transpose so that 2 will be the row
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < 2; j++) {
				transpose[j][i] = a[i][j];
			}
		}
		return transpose;
	}
	//display matrix
	public static void displayMatrix(double[][] ooAT) {
		for (int i = 0; i < ooAT.length; i++) {
			for (int j = 0; j < ooAT[i].length; j++) {
				System.out.println(String.format("%.10f", ooAT[i][j]));
			}
			System.out.println(" ");
		}
	}
	//display vector
	public static void displayVector(double[] aTy) {
		for (int i = 0; i < aTy.length; i++) {
			System.out.println(String.format("%.10f", aTy[i]));
		}
	}

	/*
	 * public static double determinant(double[][] a) { int R = a.length; int C =
	 * a[0].length; double sum = 0; if(R != C) {
	 * System.out.println("non-square matrix"); } else if(R == C && R == 2) { sum =
	 * a[0][0]*a[1][1] - a[0][1]*a[1][0]; } return 1/sum; } public static double[][]
	 * adjMatrix(double[][] a) { int R = a.length; int C = a[0].length; double[][]
	 * temp = a; double[][] adj = new double[R][C]; //swap position first adj[0][0]
	 * = temp[1][1]; adj[1][1] = temp[0][0]; adj[0][1] = -temp[0][1]; adj[1][0] =
	 * -temp[1][0]; return adj; } public static double[][] invert(double[][] a) {
	 * int R = a.length; int C = a[0].length; double[][] invert = new double[R][C];
	 * double dert = determinant(a); for(int i = 0; i < R; i++) { for(int j = 0; j <
	 * C; j++) { invert[i][j] = dert * a[i][j]; } } return invert; }
	 */
	//invert a matrix 
	public static double[][] invert(double a[][]) {
		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];
		for (int i = 0; i < n; ++i)
			b[i][i] = 1;

		// Transform the matrix into an upper triangle
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
		for (int i = 0; i < n - 1; ++i)
			for (int j = i + 1; j < n; ++j)
				for (int k = 0; k < n; ++k)
					b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];

		// Perform backward substitutions
		for (int i = 0; i < n; ++i) {
			x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
			for (int j = n - 2; j >= 0; --j) {
				x[j][i] = b[index[j]][i];
				for (int k = j + 1; k < n; ++k) {
					x[j][i] -= a[index[j]][k] * x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}
		return x;
	}
	// Method to carry out the partial-pivoting Gaussian
	// elimination.  Here index[] stores pivoting order.
	public static void gaussian(double a[][], int index[]) {
		int n = index.length;
		double c[] = new double[n];
		// Initialize the index
		for (int i = 0; i < n; ++i)
			index[i] = i;
		// Find the rescaling factors, one from each row
		for (int i = 0; i < n; ++i) {
			double c1 = 0;
			for (int j = 0; j < n; ++j) {
				double c0 = Math.abs(a[i][j]);
				if (c0 > c1)
					c1 = c0;
			}
			c[i] = c1;
		}
		// Search the pivoting element from each column
		int k = 0;
		for (int j = 0; j < n - 1; ++j) {
			double pi1 = 0;
			for (int i = j; i < n; ++i) {
				double pi0 = Math.abs(a[index[i]][j]);
				pi0 /= c[index[i]];
				if (pi0 > pi1) {
					pi1 = pi0;
					k = i;
				}
			}
			// Interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i = j + 1; i < n; ++i) {
				double pj = a[index[i]][j] / a[index[j]][j];
				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;
				// Modify other elements accordingly
				for (int l = j + 1; l < n; ++l)
					a[index[i]][l] -= pj * a[index[j]][l];
			}
		}
	}
}
