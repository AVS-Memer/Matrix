import java.util.Arrays;
import java.util.function.*;
public class Matrix {
  public static final Matrix EMPTY = new Matrix();
  double[][] values; public double det; int[] size;
  public Matrix() {
    values = new double[0][0];
    size = new int[]{0,0};
    det = 1;
  }
  public Matrix(int i) {
    if (i < 0) {
      throw new IllegalArgumentException("Invalid Matrix Dimensions");
    }
    size = new int[]{i,i};
    values = new double[i][i];
    for (int n = 0; n < i; n++) Arrays.fill(values[n], 0);
    det = (i==0)?1:0;
  }
  public Matrix(int i, int j) {
    if (i < 0 || j < 0) {
      throw new IllegalArgumentException("Invalid Matrix Dimensions");
    }
    if (i == 0 || j == 0) size = new int[]{0,0};
    else size = new int[]{i,j};
    values = new double[i][j];
    for (int n = 0; n < i; n++) Arrays.fill(values[n], 0);
    det = (i==0||j==0)?1:0;
  }
  public Matrix(int[] a) {
    if (a.length != 2) throw new IllegalArgumentException("Matrix must have 2 Dimensions");
    if (a[0] == 0 || a[1] == 0) size = new int[]{0,0};
    else size = new int[]{a[0],a[1]};
    values = new double[a[0]][a[1]];
    for (int i = 0; i < a[0]; i++) Arrays.fill(values[i], 0);
    det = (a[0]==0||a[1]==0)?1:0;
  }
  public Matrix(double[][] v) {
    if (v.length == 0) {
      values = new double[0][0];
      size = new int[]{0,0};
      det = 1;
    } else {
      values = new double[v.length][v[0].length];
      for (int i = 0; i < v.length; i++) {
        values[i] = Arrays.copyOf(v[i], v[i].length);
      }
      size = new int[]{v.length,v[0].length};
      this.calculateDeterminant();
    }
  }
  public Matrix(Vector v) {
    values = new double[v.size][1];
    for (int i = 0; i < v.size; i++) values[i][0] = v.values[i];
    if (v.size == 0) {
      size = new int[]{0,0};
      det = 1;
    } else {
      size = new int[]{v.size,1};
      this.calculateDeterminant();
    }
  }
  public Matrix(double[][] v, boolean f) {
    values = new double[v.length][v[0].length];
    for (int i = 0; i < v.length; i++) values[i] = Arrays.copyOf(v[i], v[i].length);
    if (v.length == 0) {
      size = new int[]{0,0};
      det = 1;
    } else {
      size = new int[]{v.length,v[0].length};
      if (f) det = Double.NaN;
      else this.calculateDeterminant();
    }
  }
  public static Matrix fill(int a, int b, double x) {
    return new Matrix(a,b).map(m -> x);
  }
  public static Matrix fill(int[] a, double x) {
    if (a.length != 2) throw new IllegalArgumentException("Matrix must have 2 Dimensions");
    return new Matrix(a[0],a[1]).map(m -> x);
  }
  public static Matrix zeros(int a) {
    return new Matrix(a);
  }
  public static Matrix zeros(int a, int b) {
    return new Matrix(a,b);
  }
  public static Matrix zeros(int[] a) {
    return new Matrix(a);
  }
  public static Matrix ones(int a) {
    return new Matrix(a,a).map(m -> 1.0);
  }
  public static Matrix ones(int a, int b) {
    return new Matrix(a,b).map(m -> 1.0);
  }
  public static Matrix ones(int[] a) {
    return Matrix.fill(a,1);
  }
  public static Matrix random(int a) {
    return new Matrix(a,a).map(m -> Math.random());
  }
  public static Matrix random(int a, int b) {
    return new Matrix(a,b).map(m -> Math.random());
  }
  public static Matrix random(int[] a) {
    if (a.length != 2) throw new IllegalArgumentException("Matrix must have 2 Dimensions");
    return new Matrix(a[0],a[1]).map(m -> Math.random());
  }
  public static Matrix identity(int n) {
    double[][] s = new double[n][n];
    for (int i = 0; i < n; i++) {
      s[i][i] = 1;
    }
    Matrix r = new Matrix(s,true);
    r.det = 1;
    return r;
  }
  public static Matrix diagonal(double[] n) {
    double[][] s = new double[n.length][n.length];
    for (int i = 0; i < n.length; i++) {
      Arrays.fill(s[i],0);
      s[i][i] = n[i];
    }
    return new Matrix(s);
  }
  public static Matrix diagonal(Vector n) {
    return diagonal(n.values);
  }
  public static Matrix fromColumnVector(Vector n) {
    return new Matrix(n);
  }
  public static Matrix fromRowVector(Vector n) {
    double[][] s = new double[1][n.size];
    s[0] = Arrays.copyOf(n.values, n.size);
    return new Matrix(s);
  }
  public int[] shape() {
    return size;
  }
  public Matrix reshape(int a) {
    a = Math.abs(a);
    if ((size[0]*size[1])%a!=0) throw new IllegalArgumentException("Matrix reshape must have divisible amount of elements");
    else if (a == size[0]) return this;
    double[][] s = new double[a][size[0]*size[1]/a];
    for (int i = 0; i < size[0]*size[1]; i++) {
      s[a*i/(size[0]*size[1])][i%(size[0]*size[1]/a)] = values[i/size[1]][i%size[1]];
    }
    return new Matrix(s);
  }
  public Matrix reshape(int a, int b) {
    a = Math.abs(a);
    b = Math.abs(b);
    if (a*b != size[0]*size[1]) throw new IllegalArgumentException("Matrix reshape must have same amount of elements");
    else if (a == size[0]) return this;
    double[][] s = new double[a][b];
    for (int i = 0; i < a*b; i++) {
      s[i/b][i%b] = values[i/size[1]][i%size[1]];
    }
    return new Matrix(s);
  }
  public Matrix reshape(int[] a) {
    if (a.length != 2) throw new IllegalArgumentException("Matrix must have 2 Dimensions");
    return this.reshape(a[0],a[1]);
  }
  public Matrix horizontal_augment(double[] m2) {
    if (size[0] != m2.length) {
      throw new IllegalArgumentException("Matrices must have the same number of rows for horizontal augmentation.");
    }
    double[][] r = new double[size[0]][size[1]+1];
    for (int i = 0; i < size[0]; i++) {
      System.arraycopy(values[i], 0, r[i], 0, size[1]);
      r[i][size[1]] = m2[i];
    }
    return new Matrix(r);
  }
  public Matrix horizontal_augment(Matrix m2) {
    if (size[0] != m2.size[0]) {
      throw new IllegalArgumentException("Matrices must have the same number of rows for horizontal augmentation.");
    }
    double[][] r = new double[size[0]][size[1]+m2.size[1]];
    for (int i = 0; i < size[0]; i++) {
      System.arraycopy(values[i], 0, r[i], 0, size[1]);
      System.arraycopy(m2.values[i], 0, r[i], size[1], m2.size[1]);
    }
    return new Matrix(r);
  }
  public Matrix vertical_augment(double[] m2) {
    if (size[1] != m2.length) {
      throw new IllegalArgumentException("Matrices must have the same number of columns for vertical augmentation.");
    }
    double[][] r = new double[size[0]+1][size[1]];
    for (int i = 0; i < size[0]; i++) {
      System.arraycopy(values[i], 0, r[i], 0, size[1]);
    }
    System.arraycopy(m2, 0, r[size[0]], 0, size[1]);
    return new Matrix(r);
  }
  public Matrix vertical_augment(Matrix m2) {
    if (size[1] != m2.size[1]) {
      throw new IllegalArgumentException("Matrices must have the same number of columns for vertical augmentation.");
    }
    double[][] r = new double[size[0]+m2.size[0]][size[1]];
    for (int i = 0; i < size[0]; i++) {
      System.arraycopy(values[i], 0, r[i], 0, size[1]);
    }
    for (int i = 0; i < m2.size[0]; i++) {
      System.arraycopy(m2.values[i], 0, r[i+size[0]], 0, size[1]);
    }
    return new Matrix(r);
  }
  public double sum() {
    double sum = 0;
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        sum += values[i][j];
      }
    }
    return sum;
  }
  public static Matrix block(Matrix[][] blocks) {
    // Construct each row of blocks horizontally
    Matrix[] horizontalRows = new Matrix[blocks.length];
    for (int i = 0; i < blocks.length; i++) {
      Matrix row = blocks[i][0];
      for (int j = 1; j < blocks[i].length; j++) {
        row = row.horizontal_augment(blocks[i][j]);
      }
      horizontalRows[i] = row;
    }
    // Now stack rows vertically
    Matrix full = horizontalRows[0];
    for (int i = 1; i < horizontalRows.length; i++) {
      full = full.vertical_augment(horizontalRows[i]);
    }
    return full;
  }
  public Matrix kronecker(Matrix m2) {
    double[][] s = new double[size[0]*m2.size[0]][size[1]*m2.size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        double[][] t = m2.multiply(values[i][j]).values;
        for (int k = 0; k < m2.size[0]; k++) {
          System.arraycopy(t[k],0,s[m2.size[0]*i+k],size[1]*j,m2.size[1]);
        }
      }
    }
    return new Matrix(s);
  }
  public Matrix submatrix(int rowStart, int rowEnd, int colStart, int colEnd) {
    if (rowStart < 0 || rowEnd > size[0] || colStart < 0 || colEnd > size[1] || rowStart > rowEnd || colStart > colEnd) {
      throw new IllegalArgumentException("Invalid submatrix bounds");
    }
    double[][] sub = new double[rowEnd-rowStart][colEnd-colStart];
    for (int i = rowStart; i < rowEnd; i++) {
      System.arraycopy(values[i], colStart, sub[i-rowStart], 0, colEnd-colStart);
    }
    return new Matrix(sub);
  }
  public Matrix submatrix(int rowEnd, int colStart, int colEnd) {
    return submatrix(0, rowEnd, colStart, colEnd);
  }
  public Matrix submatrix(int rowEnd, int colEnd) {
    return submatrix(0, rowEnd, 0, colEnd);
  }
  public Matrix submatrix(int rowEnd) {
    return submatrix(0, rowEnd, 0, size[1]-1);
  }
  public double get(int r, int c) {
    if (r < 0 || r >= size[0]) throw new IllegalArgumentException("Invalid row index");
    if (c < 0 || c >= size[1]) throw new IllegalArgumentException("Invalid column index");
    return values[r][c];
  }
  public double[] getRow(int r) {
    if (r >= size[0]) throw new IllegalArgumentException("Invalid row index");
    return values[r];
  }
  public double[] getRow(int r, int start) {
    return this.getRow(r,start,size[1]);
  }
  public double[] getRow(int r, int start, int end) {
    if (r >= size[0]) throw new IllegalArgumentException("Invalid row index");
    else if (start < 0 || end > size[1] || start > end) throw new IllegalArgumentException("Invalid start or end index");
    return Arrays.copyOfRange(values[r],start,end);
  }
  public Vector getRowVector(int r) {
    return new Vector(this.getRow(r));
  }
  public Vector getRowVector(int r, int end) {
    return new Vector(this.getRow(r,0,end));
  }
  public Vector getRowVector(int r, int start, int end) {
    return new Vector(this.getRow(r,start,end));
  }
  public double[] getColumn(int c) {
    return this.getColumn(c,0,size[0]);
  }
  public double[] getColumn(int c, int start) {
    return this.getColumn(c,start,size[0]);
  }
  public double[] getColumn(int c, int start, int end) {
    if (c >= size[1]) throw new IllegalArgumentException("Invalid column index");
    else if (start < 0 || end > size[0] || start > end) throw new IllegalArgumentException("Invalid start or end index");
    double[] r = new double[end-start];
    for (int i = start; i < end; i++) r[i-start] = values[i][c];
    return r;
  }
  public Vector getColumnVector(int c) {
    return new Vector(this.getColumn(c,0,size[0]));
  }
  public Vector getColumnVector(int c, int start) {
    return new Vector(this.getColumn(c,start,size[0]));
  }
  public Vector getColumnVector(int c, int start, int end) {
    return new Vector(this.getColumn(c,start,end));
  }
  public double[] getDiagonal() {
    return this.getDiagonal(0,size[0]);
  }
  public double[] getDiagonal(int start) {
    return this.getDiagonal(start,size[0]);
  }
  public double[] getDiagonal(int start, int end) {
    if (!this.isSquare()) throw new IllegalArgumentException("Nonsquare matrix does not have diagonal");
    else if (start < 0 || end > size[0] || start > end) throw new IllegalArgumentException("Invalid start or end index");
    double[] r = new double[end-start];
    for (int i = start; i < end; i++) r[i-start] = values[i][i];
    return r;
  }
  public Vector getDiagonalVector() {
    return new Vector(this.getDiagonal(0,size[0]));
  }
  public Vector getDiagonalVector(int start) {
    return new Vector(this.getDiagonal(start,size[0]));
  }
  public Vector getDiagonalVector(int start, int end) {
    return new Vector(this.getDiagonal(start,end));
  }
  public void set(int r, int c, double x) {
    if (r >= size[0]) throw new IllegalArgumentException("Invalid row index");
    if (c >= size[1]) throw new IllegalArgumentException("Invalid column index");
    values[r][c] = x;
    this.calculateDeterminant();
  }
  public void setRow(int r, double[] v) {
    if (r >= size[0]) throw new IllegalArgumentException("Invalid row index");
    if (v.length != size[1]) throw new IllegalArgumentException("Invalid row length");
    values[r] = v;
    this.calculateDeterminant();
  }
  public void setRow(int r, Vector v) {
    this.setRow(r,v.values);
  }
  public void setColumn(int c, double[] v) {
    if (c >= size[1]) throw new IllegalArgumentException("Invalid row index");
    if (v.length != size[0]) throw new IllegalArgumentException("Invalid column length");
    for (int i = 0; i < size[0]; i++) values[i][c] = v[i];
    this.calculateDeterminant();
  }
  public void setColumn(int c, Vector v) {
    this.setColumn(c,v.values);
  }
  public void swapRows(int r1, int r2) {
    if (r1 >= size[0] || r2 >= size[0]) throw new IllegalArgumentException("Invalid row indexes for swapping");
    if (r1 != r2) {
      double[] temp = Arrays.copyOf(values[r1],size[1]);
      values[r1] = Arrays.copyOf(values[r2],size[1]);
      values[r2] = Arrays.copyOf(temp,size[1]);
    }
    det *= -1;
  }
  public void swapColumns(int r1, int r2) {
    if (r1 >= size[1] || r2 >= size[1]) throw new IllegalArgumentException("Invalid column indexes for swapping");
    double temp;
    if (r1 != r2) {
      for (int i = 0; i < size[0]; i++) {
        temp = values[r1][i];
        values[r1][i] = values[r2][i];
        values[r2][i] = temp;
      }
    }
    det *= -1;
  }
  public Vector[] toVectorArray() {
    return this.toVectorArray(false);
  }
  public Vector[] toVectorArray(boolean row) {
    if (row) {
      Vector[] r = new Vector[size[0]];
      for (int i = 0; i < size[0]; i++) r[i] = new Vector(values[i]);
    } else {
      Vector[] r = new Vector[size[1]];
      for (int i = 0; i < size[1]; i++) r[i] = this.getColumnVector(i);
    }
    return r;
  }
  public Matrix scaleRow(int row, double number) {
    if (row >= size[0]) throw new IllegalArgumentException("Invalid row index for scaling");
    double[][] r = Arrays.copyOf(values,size[0]);
    for (int i = 0; i < size[1]; i++) r[row][i] *= number;
    return new Matrix(r);
  }
  public Matrix scaleColumn(int col, double number) {
    if (col >= size[0]) throw new IllegalArgumentException("Invalid column index for scaling");
    double[][] r = Arrays.copyOf(values,size[0]);
    for (int i = 0; i < size[1]; i++) r[i][col] *= number;
    return new Matrix(r);
  }
  public Matrix rowAddition(int r1, int r2, int number) {
    if (r1 >= size[0] || r2 >= size[0]) throw new IllegalArgumentException("Invalid row indexes for row addition");
    double[][] r = Arrays.copyOf(values,size[0]);
    for (int i = 0; i < size[1]; i++) r[r1][i] += r[r2][i]*number;
    return new Matrix(r);
  }
  public Matrix rowAddition(int r1, int r2) {
    return rowAddition(r1,r2,1);
  }
  public Matrix colAddition(int c1, int c2, int number) {
    if (c1 >= size[0] || c2 >= size[0]) throw new IllegalArgumentException("Invalid column indexes for row addition");
    double[][] r = Arrays.copyOf(values,size[0]);
    for (int i = 0; i < size[1]; i++) r[i][c1] += r[i][c2]*number;
    return new Matrix(r);
  }
  public Matrix colAddition(int c1, int c2) {
    return colAddition(c1,c2,1);
  }
  public Matrix reducedRowEchelon() {
    Matrix r = this.gaussianElimination(); // Start from REF
    int rows = r.size[0];
    int cols = r.size[1];
    for (int i = rows - 1; i >= 0; i--) {
      // Find pivot in row i
      int pivotCol = -1;
      for (int j = 0; j < cols; j++) {
        if (Math.abs(r.values[i][j]) > 1e-10) {
          pivotCol = j;
          break;
        }
      }
      if (pivotCol == -1) continue; // skip zero rows
      // Normalize pivot row so pivot is 1
      double pivotVal = r.values[i][pivotCol];
      for (int j = pivotCol; j < cols; j++) r.values[i][j] /= pivotVal;
      // Eliminate above the pivot
      for (int k = 0; k < i; k++) {
        double factor = r.values[k][pivotCol];
        for (int j = pivotCol; j < cols; j++) {
          r.values[k][j] -= factor*r.values[i][j];
        }
      }
    }
    return r;
  }
  public double mean() {
    return this.sum()/(size[0]*size[1]);
  }
  public double median() {
    double[] allValues = this.flattenRowMajor().values[0];
    Arrays.sort(allValues);
    return (allValues.length % 2 == 1) ? allValues[allValues.length / 2] : (allValues[allValues.length / 2 - 1] + allValues[allValues.length / 2]) / 2.0;
  }
  public Matrix center() {
    return this.subtract(this.mean());
  }
  public double std() {
    return Math.sqrt(this.center().hadamard(this.center()).mean());
  }
  private void calculateDeterminant() {
    if (this.isSquare()) {
      if (this.isDiagonal()) {
        det = 1;
        for (double i : this.getDiagonal()) det *= i;
      } else {
        DeterminantResult lup = this.lupDecomposeWithSwapCount();
        if (lup == null) det = 0;
        else {
          det = lup.swapCount?1:-1;
          int n = size[0];
          for (int i = 0; i < n; i++) det *= lup.U.values[i][i];
        }
      }
    } else det = Double.NaN;
  }
  public boolean isSquare() {
    return size[0]==size[1];
  }
  public boolean isInvertible() {
    return det!=0;
  }
  public boolean isSingular() {
    return det==0;
  }
  public boolean isEqual(Matrix m2) {
    if (!Arrays.equals(size,m2.size)) return false;
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        double diff = Math.abs(values[i][j] - m2.values[i][j]);
        if (diff > 1e-9) return false;
      }
    }
    return true;
  }
  public Matrix transpose() {
    double[][] s = new double[size[1]][size[0]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        s[j][i] = values[i][j];
      }
    }
    return new Matrix(s);
  }
  public boolean isSymmetric() {
    return this.isSquare() && this.isEqual(this.transpose());
  }
  public boolean isOrthogonal() {
    return this.isSquare() && this.multiply(this.transpose()).isIdentity() && this.isOrthonormal();
  }
  public boolean isOrthonormal() {
    return this.transpose().multiply(this).isIdentity();
  }
  public boolean isEmpty() {
    return this.isEqual(Matrix.EMPTY);
  }
  public boolean isIdentity() {
    return this.isSquare() && this.isEqual(Matrix.identity(size[0]));
  }
  public boolean isDiagonal() {
    return this.isSquare() && this.isEqual(this.hadamard(Matrix.identity(size[0])));
  }
  public boolean isUpperTriangular() {
    if (!this.isSquare()) return false;
    for (int i = 1; i < size[0]; i++) {
      for (int j = 0; j < i; j++) {
        if (Math.abs(values[i][j]) > 1e-10) return false;
      }
    }
    return true;
  }
  public boolean isLowerTriangular() {
    if (!this.isSquare()) return false;
    for (int i = 0; i < size[0]; i++) {
      for (int j = i+1; j < size[1]; j++) {
        if (Math.abs(values[i][j]) > 1e-10) return false;
      }
    }
    return true;
  }
  public boolean isPositiveDefinite() {
    try {
      this.choleskyDecompose();
      return true;
    } catch (Exception e) {return false;}
  }
  public Matrix add(double m2) {
    return this.map(a -> a+m2);
  }
  public Matrix add(Matrix m2) {
    if (!Arrays.equals(size,m2.size)) {
      throw new IllegalArgumentException("Matrices must be same size for addition.");
    }
    return this.map(m2,(a,b) -> a+b);
  }
  public Matrix subtract(double m2) {
    return this.map(a -> a-m2);
  }
  public Matrix subtract(Matrix m2) {
    if (!Arrays.equals(size,m2.size)) {
      throw new IllegalArgumentException("Matrices must be same size for subtraction.");
    }
    return this.map(m2,(a,b) -> a-b);
  }
  public Matrix multiply(double number) {
    if (!Double.isNaN(det)) det *= Math.pow(number,size[0]);
    return this.map(a -> a*number);
  }
  public Matrix negate() {
    return this.map(a -> -a);
  }
  private Matrix standardMultiply(Matrix m2) {
    if (size[1] != m2.size[0]) {
      throw new IllegalArgumentException("Width of first matrix must be equal to height of second matrix for matrix multiplication.");
    }
    double[][] s = new double[size[0]][m2.size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < m2.size[1]; j++) {
        s[i][j] = 0;
        for (int k = 0; k < size[1]; k++) {
          s[i][j] += values[i][k]*m2.values[k][j];
        }
      }
    }
    return new Matrix(s);
  }
  public Matrix multiply(Matrix other) {
    if (size[1] != other.size[0]) {
      throw new IllegalArgumentException("Inner matrix dimensions must agree.");
    }
    if (this.isDiagonal()) {
      Matrix result = other.clone();
      for (int i = 0; i < size[0]; i++) result = result.scaleRow(i,values[i][i]);
      return result;
    } else if (other.isDiagonal()) {
      Matrix result = this.clone();
      for (int i = 0; i < size[1]; i++) result = result.scaleColumn(i,other.values[i][i]);
      return result;
    }
    int n = Math.max(Math.max(size[0], size[1]), Math.max(other.size[0], other.size[1]));
    int m = 1;
    while (m < n) m <<= 1;
    Matrix A = this.padToSize(m, m);
    Matrix B = other.padToSize(m, m);
    Matrix C = A.strassenMultiply(B);
    return C.submatrix(size[0], other.size[1]);
  }
  public Matrix padToSize(int newRows, int newCols) {
    if (newRows < size[0] || newCols < size[1]) {
      throw new IllegalArgumentException("New size must be greater than or equal to current size.");
    }
    Matrix padded = new Matrix(newRows,newCols);
    for (int i = 0; i < size[0]; i++) System.arraycopy(values[i], 0, padded.values[i], 0, size[1]);
    return padded;
  }
  private Matrix strassenMultiply(Matrix B) {
    if (size[0] <= 64) { // base case threshold
      return this.standardMultiply(B);
    }
    int newSize = size[0] / 2;
    Matrix A11 = this.submatrix(0, 0, newSize, newSize);
    Matrix A12 = this.submatrix(0, newSize, newSize, newSize*2);
    Matrix A21 = this.submatrix(newSize, 0, newSize*2, newSize);
    Matrix A22 = this.submatrix(newSize, newSize, newSize*2, newSize*2);
    Matrix B11 = B.submatrix(0, 0, newSize, newSize);
    Matrix B12 = B.submatrix(0, newSize, newSize, newSize*2);
    Matrix B21 = B.submatrix(newSize, 0, newSize*2, newSize);
    Matrix B22 = B.submatrix(newSize, newSize, newSize*2, newSize*2);
    Matrix M1 = A11.add(A22).strassenMultiply(B11.add(B22));
    Matrix M2 = A21.add(A22).strassenMultiply(B11);
    Matrix M3 = A11.strassenMultiply(B12.subtract(B22));
    Matrix M4 = A22.strassenMultiply(B21.subtract(B11));
    Matrix M5 = A11.add(A12).strassenMultiply(B22);
    Matrix M6 = A21.subtract(A11).strassenMultiply(B11.add(B12));
    Matrix M7 = A12.subtract(A22).strassenMultiply(B21.add(B22));
    Matrix C11 = M1.add(M4).subtract(M5).add(M7);
    Matrix C12 = M3.add(M5);
    Matrix C21 = M2.add(M4);
    Matrix C22 = M1.subtract(M2).add(M3).add(M6);
    return Matrix.block(new Matrix[][]{{C11, C12}, {C21, C22}});
  }
  public Matrix hadamard(Matrix m2) {
    if (!Arrays.equals(size,m2.size)) {
      throw new IllegalArgumentException("Matrices must be of the same size to find Hadamard Product.");
    }
    return this.map(m2,(a,b) -> a*b);
  }
  public Matrix minor(int row, int col) {
    if (size[0] == 0) throw new IllegalArgumentException("Empty matrix does not have minor");
    double[][] result = new double[size[0]-1][size[1]-1];
    int r = 0;
    for (int i = 0; i < size[0]; i++) {
      if (i == row) continue;
      int c = 0;
      for (int j = 0; j < size[1]; j++) {
        if (j == col) continue;
        result[r][c++] = values[i][j];
      }
      r++;
    }
    return new Matrix(result);
  }
  public Matrix cofactor() {
    if (!this.isSquare()) {
      throw new IllegalArgumentException("Cofactor matrix can only be computed for square matrices.");
    }
    double[][] cofactors = new double[size[0]][size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        cofactors[i][j] = Math.pow(-1, i+j) * this.minor(i,j).det;
      }
    }
    return new Matrix(cofactors);
  }
  public Matrix adjugate() {
    return this.cofactor().transpose();
  }
  public Matrix inverse() {
    if (det == 0) throw new IllegalArgumentException("Matrix with zero determinant does not have inverse");
    if (this.isDiagonal()) return this.diagonalInverse();
    return this.adjugate().divide(det);
  }
  private Matrix diagonalInverse() {
    Matrix inv = new Matrix(size[1],size[0]);
    for (int i = 0; i < Math.min(size[0], size[1]); i++) {
      if (Math.abs(values[i][i]) > 1e-10) inv.values[i][i] = 1.0/values[i][i];
    }
    return inv;
  }
  public Matrix pseudoInverse() {
    Matrix[] svd = this.svd(); // A = U * S * Vᵀ. Compute S⁺: invert non-zero singular values. Compute A⁺ = V * S⁺ * Uᵗ
    return svd[2].multiply(svd[1].diagonalInverse()).multiply(svd[0].transpose());
  }
  public Matrix divide(double number) {
    if (det != Double.NaN) det /= Math.pow(number,size[0]);
    return this.map(a -> a/number);
  }
  public Matrix divide(Matrix m2) {
    if (size[1] != m2.size[0]) {
      throw new IllegalArgumentException("Width of first matrix must be equal to height of second matrix for matrix multiplication.");
    } else if (m2.det == 0) {
      throw new IllegalArgumentException("Divisor matrix cannot have zero determinant.");
    }
    if (det != Double.NaN) det /= m2.det;
    return this.multiply(m2.inverse());
  }
  public double trace() {
    if (size[0] != size[1]) {
      throw new IllegalArgumentException("Matrix must be square to trace.");
    }
    double n = 0;
    for (int i = 0; i < size[0]; i++) n += values[i][i];
    return n;
  }
  public Matrix power(int exponent) {
    if (size[0] != size[1]) {
      throw new IllegalArgumentException("Matrix must be square for exponentiation.");
    }
    Matrix result = Matrix.identity(size[0]);
    if (this.isDiagonal()) {
      for (int i = 0; i < size[0]; i++) result.values[i][i] = Math.pow(values[i][i], exponent);
      return result;
    }
    Matrix base = this;
    while (exponent > 0) {
      if ((exponent & 1) == 1) {
        result = result.multiply(base);
      }
      base = base.multiply(base);
      exponent >>= 1;  // Divide exponent by 2
    }
    return result;
  }
  public Matrix log(int maxIterations, double tolerance) {
    if (!this.isSquare()) {
      throw new IllegalArgumentException("Matrix must be square to compute logarithm.");
    }
    Matrix I = identity(size[0]);
    Matrix X = this.subtract(I); // X = A - I
    Matrix logA = new Matrix(size[0]);
    Matrix term = X.clone();
    for (int k = 1; k <= maxIterations; k++) {
      Matrix currentTerm = term.multiply((k%2==0?-1.0:1.0)/k);
      logA = logA.add(currentTerm);
      term = term.multiply(X); // X^k
      if (currentTerm.norm() < tolerance) break;
    }
    return logA;
  }
  public Matrix gaussianElimination() {
    Matrix data = this.clone();
    int pivotRow = 0;
    for (int col = 0; col < size[1] && pivotRow < size[0]; col++) {
      int maxRow = pivotRow;
      for (int i = pivotRow + 1; i < size[0]; i++) {
        if (Math.abs(data.values[i][col]) > Math.abs(data.values[maxRow][col])) maxRow = i;
      }
      // Skip if column is zero
      if (Math.abs(data.values[maxRow][col]) < 1e-9) continue;
      // Swap current row with pivot row
      data.swapRows(pivotRow, maxRow);
      // Eliminate entries below pivot
      for (int i = pivotRow + 1; i < size[0]; i++) {
        double factor = data.values[i][col] / data.values[pivotRow][col];
        for (int j = col; j < size[1]; j++) data.values[i][j] -= factor * data.values[pivotRow][j];
      }
      pivotRow++;
    }
    return data;  // Return new Matrix in REF
  }
  public int rank() {
    double[][] temp = new Matrix(values).values;
    int rank = size[1];  // Start with full column count
    for (int row = 0; row < rank; row++) {
      // If non-zero diagonal element
      if (temp[row][row] != 0) {
        for (int col = 0; col < size[0]; col++) {
          if (col != row) {
            double mult = temp[col][row] / temp[row][row];
            for (int i = 0; i < rank; i++) {
              temp[col][i] -= mult * temp[row][i];
            }
          }
        }
      } else {
        // Find a row to swap
        boolean reduce = true;
        for (int i = row + 1; i < size[0]; i++) {
          if (temp[i][row] != 0) {
            // Swap rows
            double[] tmp = temp[row];
            temp[row] = temp[i];
            temp[i] = tmp;
            reduce = false;
            break;
          }
        }
        if (reduce) {
          // Reduce rank
          rank--;
          // Copy last column into current column
          for (int i = 0; i < size[0]; i++) {
            temp[i][row] = temp[i][rank];
          }
          row--; // Stay in same row
        }
      }
    }
    return rank;
  }
  public Matrix[] qrDecompose() {
    Matrix A = this.clone();
    Matrix Q = Matrix.identity(size[0]);
    for (int k = 0; k < Math.min(size[0], size[1]); k++) {
      // Get the vector x (k-th column below diagonal)
      Vector x = new Vector(size[0]-k);
      for (int i = k; i < size[0]; i++) {
        x.values[i-k] = A.values[i][k];
      }
      // Compute the norm and the Householder vector v
      double normX = x.norm();
      if (normX == 0) continue;
      Vector e1 = new Vector(size[0]-k);
      e1.values[0] = 1;
      Vector v = x.subtract(e1.multiply(Math.copySign(normX, x.values[0])));
      double vNorm = v.norm();
      if (vNorm < 1e-10) continue;
      v = v.multiply(1.0 / vNorm);
      // Apply the reflection to A[k:m][k:n]
      for (int j = k; j < size[1]; j++) {
        double dot = v.dotProduct(A.getColumnVector(j,k));
        for (int i = k; i < size[0]; i++) A.values[i][j] -= 2*v.values[i-k]*dot;
      }
      // Apply the reflection to Q (accumulating Q^T)
      for (int j = 0; j < size[0]; j++) {
        double dot = v.dotProduct(Q.getColumnVector(j,k));
        for (int i = k; i < size[0]; i++) Q.values[i][j] -= 2*v.values[i-k]*dot;
      }
    }
    Q = Q.transpose(); // Q was built as Q^T
    return new Matrix[]{Q, A};
  }
  public Matrix[] luDecompose() {
    if (!this.isSquare()) {
      throw new IllegalArgumentException("Matrix must be square");
    }
    Matrix L = Matrix.identity(size[0]);
    Matrix U = new Matrix(size[0]);
    for (int k = 0; k < size[0]; k++) {
      // Compute U[k][j] for j >= k
      for (int j = k; j < size[0]; j++) {
        double sum = L.getRowVector(k, 0, k).dotProduct(U.getColumnVector(j, 0, k));
        U.values[k][j] = values[k][j] - sum;
      }
      // Compute L[i][k] for i > k
      for (int i = k + 1; i < size[0]; i++) {
        double sum = L.getRowVector(i,0,k).dotProduct(U.getColumnVector(k,0,k));
        for (int s = 0; s < k; s++) {
          sum += L.values[i][s] * U.values[s][k];
        }
        if (U.values[k][k] == 0) {
          throw new ArithmeticException("Zero pivot encountered at index " + k);
        }
        L.values[i][k] = (values[i][k] - sum) / U.values[k][k];
      }
    }
    return new Matrix[] {L, U};
  }
  public Matrix[] lupDecompose() {
    if (!this.isSquare()) throw new IllegalArgumentException("Matrix must be square");
    else if (this.isSingular()) throw new IllegalArgumentException("Singular matrix cannot be decomposed");
    Matrix A = this.clone(); // Work on a copy
    Matrix L = Matrix.identity(size[0]);
    Matrix U = new Matrix(size[0]);
    Matrix P = Matrix.identity(size[0]);
    for (int k = 0; k < size[0]; k++) {
      // Pivot: find the row with the largest element in column k
      int maxRow = k;
      double maxVal = Math.abs(A.values[k][k]);
      for (int i = k + 1; i < size[0]; i++) {
        double val = Math.abs(A.values[i][k]);
        if (val > maxVal) {
          maxVal = val;
          maxRow = i;
        }
      }
      // Swap rows in A and P (and L partially)
      if (maxRow != k) {
        A.swapRows(k, maxRow);
        P.swapRows(k, maxRow);
        for (int j = 0; j < k; j++) {
          double tempL = L.values[k][j];
          L.values[k][j] = L.values[maxRow][j];
          L.values[maxRow][j] = tempL;
        }
      }
      // Compute U[k][j] for j >= k
      for (int j = k; j < size[0]; j++) {
        double sum = L.getRowVector(k, 0, k).dotProduct(U.getColumnVector(j, 0, k));
        U.values[k][j] = A.values[k][j] - sum;
      }
      // Compute L[i][k] for i > k
      for (int i = k + 1; i < size[0]; i++) {
        double sum = L.getRowVector(i,0,k).dotProduct(U.getColumnVector(k,0,k));
        L.values[i][k] = (A.values[i][k] - sum) / U.values[k][k];
      }
    }
    return new Matrix[]{L,U,P};
  }
  private DeterminantResult lupDecomposeWithSwapCount() {
    if (!this.isSquare()) throw new IllegalArgumentException("Matrix must be square");
    Matrix A = this.clone(); // Work on a copy
    Matrix L = Matrix.identity(size[0]);
    Matrix U = new Matrix(size[0]);
    Matrix P = Matrix.identity(size[0]);
    boolean swapCount = true;
    for (int k = 0; k < size[0]; k++) {
      // Pivot: find the row with the largest element in column k
      int maxRow = k;
      double maxVal = Math.abs(A.values[k][k]);
      for (int i = k + 1; i < size[0]; i++) {
        double val = Math.abs(A.values[i][k]);
        if (val > maxVal) {
          maxVal = val;
          maxRow = i;
        }
      }
      // Swap rows in A and P (and L partially)
      if (maxRow != k) {
        A.swapRows(k, maxRow);
        P.swapRows(k, maxRow);
        for (int j = 0; j < k; j++) {
          double tempL = L.values[k][j];
          L.values[k][j] = L.values[maxRow][j];
          L.values[maxRow][j] = tempL;
        }
        swapCount = !swapCount;
      }
      // Compute U[k][j] for j >= k
      for (int j = k; j < size[0]; j++) {
        double sum = L.getRowVector(k, 0, k).dotProduct(U.getColumnVector(j, 0, k));
        U.values[k][j] = A.values[k][j] - sum;
      }
      // Compute L[i][k] for i > k
      for (int i = k + 1; i < size[0]; i++) {
        double sum = L.getRowVector(i,0,k).dotProduct(U.getColumnVector(k,0,k));
        if (Math.abs(U.values[k][k]) < 1e-10) {
          return null;
        }
        L.values[i][k] = (A.values[i][k] - sum) / U.values[k][k];
      }
    }
    return new DeterminantResult(U,swapCount);
  }
  public Matrix choleskyDecompose() {
    if (!this.isSymmetric()) throw new IllegalArgumentException("Matrix must be square and symmetric.");
    Matrix L = new Matrix(size[0]); // Zero matrix
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j <= i; j++) {
        // Compute sum of L[i][k] * L[j][k] for k = 0..j-1
        double sum = L.getRowVector(i,0,j).dotProduct(L.getRowVector(j,0,j));
        if (i == j) {
          double diag = values[i][i]-sum;
          if (diag <= 0.0) {
            throw new ArithmeticException("Matrix is not positive definite.");
          }
          L.values[i][j] = Math.sqrt(diag);
        } else {
          if (Math.abs(L.values[j][j]) < 1e-12) {
            throw new ArithmeticException("Near-zero pivot encountered at L[" + j + "][" + j + "]. Matrix may not be positive definite.");
          }
          L.values[i][j] = (values[i][j]-sum) / L.values[j][j];
        }
      }
    }
    return L;
  }
  public double[] eigenvalues() {
    return this.eigenvalues(100);
  }
  public double[] eigenvalues(int maxIterations) {
    if (!this.isSquare()) throw new IllegalArgumentException("Matrix must be square to compute eigenvalues.");
    else if (this.isEmpty()) return new double[0];
    else if (this.isDiagonal()) return this.getDiagonal();
    Matrix A = this.clone(); // Work on a copy
    for (int iter = 0; iter < maxIterations; iter++) {
      Matrix[] qr = A.qrDecompose(); // A = Q * R
      Matrix Q = qr[0];
      Matrix R = qr[1];
      A = R.multiply(Q); // A_{k+1} = R * Q
      // Check for convergence: off-diagonal elements should become small
      Vector diagonalVector = A.getDiagonalVector();
      double offDiagonalNorm = A.sumOfSquares()-diagonalVector.dotProduct(diagonalVector);
      if (Math.sqrt(offDiagonalNorm) < 1e-9) break;
    }
    // After convergence, eigenvalues are on the diagonal
    return A.getDiagonal();
  }
  public Vector[] eigenvectors() {
    return this.eigenvectors(this.eigenvalues());
  }
  private Vector[] eigenvectors(double[] eigenvalues) {
    if (size[0] != size[1]) {
      throw new IllegalArgumentException("Matrix must be square to compute eigenvectors.");
    }
    Matrix I = Matrix.identity(size[0]);
    if (this.isDiagonal()) return I.toVectorArray();
    Vector[] eigenvectors = new Vector[eigenvalues.length];
    for (int i = 0; i < eigenvalues.length; i++) {
      double lambda = eigenvalues[i];
      Matrix shifted = this.subtract(I.multiply(lambda)); // A - λI

      // Find null space (i.e., nontrivial solution to (A - λI)v = 0)
      Matrix nullSpace = shifted.nullSpace();

      // Store the first basis vector (can be multiple eigenvectors)
      if (nullSpace.size[1] > 0) {
        eigenvectors[i] = nullSpace.getColumnVector(0);
      } else {
        // fallback: return zero vector (or throw if you want strict behavior)
        eigenvectors[i] = new Vector(size[0]); // zero vector
      }
    }
    return eigenvectors;
  }
  public Matrix nullSpace() {
    Matrix rref = this.reducedRowEchelon();
    int m = rref.size[0];
    int n = rref.size[1];
    // Track pivot columns
    boolean[] isPivot = new boolean[n];
    int pivotRow = 0;
    for (int j = 0; j < n && pivotRow < m; j++) {
      if (Math.abs(rref.values[pivotRow][j]) > 1e-10) {
        isPivot[j] = true;
        pivotRow++;
      }
    }
    // Count number of free variables
    int numFree = 0;
    for (int j = 0; j < n; j++) {
      if (!isPivot[j]) numFree++;
    }
    // Initialize result matrix: n rows, numFree columns
    Matrix nullBasis = new Matrix(n, numFree);
    int basisIndex = 0;
    for (int freeCol = 0; freeCol < n; freeCol++) {
      if (isPivot[freeCol]) continue;
      Vector vec = new Vector(n);
      vec.values[freeCol] = 1.0;
      int row = 0;
      for (int j = 0; j < n; j++) {
        if (isPivot[j]) {
          double sum = rref.getRowVector(row).dotProduct(vec);
          vec.values[j] = -sum;
          row++;
        }
      }
      // Fill this vector into the null space matrix as a column
      nullBasis.setColumn(basisIndex++, vec.values);
    }
    return nullBasis;
  }
  public double[] singularValues() {
    double[] eigenvalues = this.eigenvalues();
    for (int i = 0; i < eigenvalues.length; i++) {
      eigenvalues[i] = eigenvalues[i] > 0 ? Math.sqrt(eigenvalues[i]) : 0.0;
    }
    return eigenvalues;
  }
  public Matrix[] svd() {
    Matrix At = this.transpose();
    Matrix AtA = At.multiply(this); // A^T A
    // Step 1: Eigen-decompose A^T A to get V and singular values squared
    double[] eigVals = AtA.eigenvalues();  // maxIter and tolerance
    Vector[] Vlist = AtA.eigenvectors();   // Columns are eigenvectors

    // Step 2: Singular values are sqrt of eigenvalues (clipped at 0)
    double[] singularValues = new double[eigVals.length];
    for (int i = 0; i < eigVals.length; i++) {
      singularValues[i] = eigVals[i] > 0 ? Math.sqrt(eigVals[i]) : 0.0;
    }
    Matrix V = new Matrix(size[1]);
    for (int i = 0; i < Vlist.length; i++) V.setColumn(i, Vlist[i]);
    // Step 3: Build Σ (sigma matrix)
    Matrix Sigma = new Matrix(size[0], size[1]);
    for (int i = 0; i < Math.min(size[0], size[1]); i++) Sigma.values[i][i] = singularValues[i];

    // Step 4: Compute U = A * V * Σ^{-1}
    Matrix U = new Matrix(size[0]);
    for (int i = 0; i < V.size[1]; i++) {
      if (singularValues[i] < 1e-10) continue; // skip tiny singular values
      U.setColumn(i, V.getColumnVector(i).multiply(this).divide(singularValues[i]).values); // normalize
    }

    // Fill remaining columns of U to make it orthogonal (optional)
    // e.g., via Gram-Schmidt if U is not square

    // Return [U, Sigma, V^T]
    Matrix Vt = V.transpose();
    return new Matrix[] {U,Sigma,Vt};
  }
  public double conditionNumber() {
    if (!this.isEmpty()) return Double.NaN;
    double[] s = this.singularValues();
    double min = s[0], max = s[0];
    for (double i : s) {
      min = Math.min(min,i);
      max = Math.max(max,i);
    }
    return max/min;
  }
  public Matrix[] spectralDecompose() {
    if (!this.isSquare() || !this.isSymmetric()) {
      throw new IllegalArgumentException("Matrix must be square and symmetric for spectral decomposition.");
    }

    // Step 1: Get eigenvalues and eigenvectors
    double[] eigenvalues = this.eigenvalues();
    Vector[] eigenvectors = this.eigenvectors(); // assume orthonormal

    // Step 2: Form Q matrix from eigenvectors (as columns)
    Matrix Q = new Matrix(size[0], eigenvectors.length);
    for (int i = 0; i < eigenvectors.length; i++) Q.setColumn(i,eigenvectors[i]);
    // Step 3: Form Λ (diagonal matrix of eigenvalues)
    Matrix Lambda = Matrix.diagonal(eigenvalues);

    // Optional: You could verify that Q is orthonormal here
    // assert Q.transpose().multiply(Q).approxEquals(Matrix.identity(n), 1e-9);

    // Step 4: Return Q and Λ — user can reconstruct A ≈ Q·Λ·Qᵀ
    return new Matrix[] {Q, Lambda};
  }
  public Matrix[] schurDecompose() {
    return this.schurDecompose(100);
  }
  public Matrix[] schurDecompose(int maxIterations) {
    if (!this.isSquare()) {
      throw new IllegalArgumentException("Matrix must be square for Schur decomposition.");
    }
    int n = size[0];
    Matrix A = this.clone();               // Work on a copy
    Matrix Q_total = Matrix.identity(n);   // Accumulate Q from QR steps
    for (int iter = 0; iter < maxIterations; iter++) {
      Matrix[] qr = A.qrDecompose();      // A = Q * R
      Matrix Q = qr[0];
      Matrix R = qr[1];
      A = R.multiply(Q);                  // A_{k+1} = R * Q
      Q_total = Q_total.multiply(Q);      // Accumulate Q: Q_total = Q_total * Q
      // Check convergence: off-diagonal energy
      Vector diag = A.getDiagonalVector();
      double offDiagonal = A.sumOfSquares() - diag.dotProduct(diag);
      if (Math.sqrt(offDiagonal) < 1e-10) break;
    }
    return new Matrix[] {Q_total,A}; // A has converged to upper triangular T
  }
  public Matrix normalize() {
    return this.divide(this.norm());
  }
  public Matrix det_normalize() {
    return this.divide(det);
  }
  public Matrix maxNormalize() {
    return this.divide(this.max());
  }
  public Matrix minMaxNormalize() {
    double min = this.min();
    return this.subtract(min).divide(this.max()-min);
  }
  public Matrix zScoreNormalize() {
    if (this.std() == 0) return new Matrix(size[0],size[1]);
    else return this.subtract(this.mean()).divide(this.std());
  }
  public Matrix L1_normalize() {
    return this.divide(this.L1_norm());
  }
  public Matrix LInf_normalize() {
    return this.divide(this.LInf_norm());
  }
  public double sumOfSquares() {
    double sum = 0;
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) sum += values[i][j]*values[i][j];
    }
    return sum;
  }
  public double norm() {
    return Math.sqrt(this.sumOfSquares());
  }
  public double L1_norm() {
    double sum = 0;
    for (int j = 0; j < size[1]; j++) {
      for (int i = 0; i < size[0]; i++) sum += Math.abs(values[i][j]);
    }
    return sum;
  }
  public double LInf_norm() {
    double maxRowSum = 0;
    for (int i = 0; i < size[0]; i++) {
      double rowSum = 0;
      for (int j = 0; j < size[1]; j++) rowSum += Math.abs(values[i][j]);
      maxRowSum = Math.max(rowSum,maxRowSum);
    }
    return maxRowSum;
  }
  public double max() {
    if (size[0] == 0) throw new IllegalArgumentException("Empty matrix does not contain maximum value");
    double m = values[0][0];
    for (double[] row : values) {
      for (double v : row) m = Math.max(v,m);
    }
    return m;
  }
  public double min() {
    if (size[0] == 0) throw new IllegalArgumentException("Empty matrix does not contain minimum value");
    double m = values[0][0];
    for (double[] row : values) {
      for (double v : row) m = Math.min(v,m);
    }
    return m;
  }
  public Matrix softmax() {
    Matrix result = new Matrix(size);
    for (int i = 0; i < size[0]; i++) {
      result.setRow(i, this.getRowVector(i).softmax().values);
    }
    return result;
  }
  public double[] rowSums() {
    double[] sums = new double[size[0]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        sums[i] += values[i][j];
      }
    }
    return sums;
  }
  public double[] colSums() {
    double[] sums = new double[size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        sums[j] += values[i][j];
      }
    }
    return sums;
  }
  public double[] rowMeans() {
    double[] means = rowSums();
    for (int i = 0; i < means.length; i++) {
      means[i] /= size[1];
    }
    return means;
  }
  public double[] colMeans() {
    double[] means = colSums();
    for (int i = 0; i < means.length; i++) {
      means[i] /= size[0];
    }
    return means;
  }
  public Matrix flattenRowMajor() {
    return this.reshape(1,size[0]*size[1]);
  }
  public Matrix flattenColMajor() {
    double[][] flat = new double[1][size[0]*size[1]];
    int idx = 0;
    for (int i = 0; i < size[1]; i++) {
      for (int j = 0; j < size[0]; j++) {
        flat[0][idx++] = values[j][i];
      }
    }
    return new Matrix(flat);
  }
  public Matrix clip() {
    return this.conditionalMap((x -> x>1), (x -> 1.0), (x -> x<0?0.0:x));
  }
  public Matrix clip(double min) {
    return this.conditionalMap((x -> x<min), (x -> min), (x -> x));
  }
  public Matrix clip(double min, double max) {
    return this.conditionalMap((x -> x>max), (x -> max), (x -> x<min?min:x));
  }
  public Matrix conditionalMap(Predicate<Double> condition, Function<Double, Double> ifTrue, Function<Double, Double> ifFalse) {
    return this.map(a -> condition.test(a) ? ifTrue.apply(a) : ifFalse.apply(a));
  }
  public Matrix map(Function<Double, Double> f) {
    double[][] r = new double[size[0]][size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        r[i][j] = f.apply(values[i][j]);
      }
    }
    return new Matrix(r);
  }
  public Matrix map(Matrix m2, BiFunction<Double, Double, Double> func) {
    if (!Arrays.equals(size, m2.size)) {
      throw new IllegalArgumentException("Matrices must be the same size for elementwise map.");
    }
    double[][] r = new double[size[0]][size[1]];
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        r[i][j] = func.apply(values[i][j], m2.values[i][j]);
      }
    }
    return new Matrix(r);
  }
  public Matrix applyActivation(ActivationFunction f) {
    return this.map(f.get());
  }
  @Override
  public Matrix clone() {
    double[][] clonedValues = new double[size[0]][size[1]];
    for (int i = 0; i < size[0]; i++) clonedValues[i] = Arrays.copyOf(values[i], size[1]);
    Matrix copy = new Matrix(clonedValues,true);
    copy.det = det;
    return copy;
  }
  @Override
  public String toString() {
    return this.toString(0);
  }
  public String toString(int decimals) {
    int maxLength = 0;
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < size[1]; j++) {
        maxLength = Math.max(maxLength, String.format("%."+decimals+"f", values[i][j]).length());
      }
    }
    StringBuilder s = new StringBuilder();
    s.append("[");
    for (int i = 0; i < size[0]; i++) {
      if (size[0] > 1) s.append("\n  ");
      s.append("[ ");
      for (int j = 0; j < size[1]; j++) {
        s.append(String.format("%"+maxLength+"."+decimals+"f", values[i][j]));
        if (j < size[1] - 1) s.append("  ");
      }
      s.append(" ]");
    }
    if (size[0] > 1) s.append("\n");
    s.append("]");
    return s.toString();
  }
  private class DeterminantResult {
    final Matrix U;
    final boolean swapCount;
    DeterminantResult(Matrix U, boolean swapCount) {
      this.U = U;
      this.swapCount = swapCount;
    }
  }
}

enum ActivationFunction {
  SIGMOID(x -> 1/(1+Math.exp(-x))),
  HARDSIGMOID(x -> Math.max(0, Math.min(1, 0.2 * x + 0.5))),
  BIPOLARSIGMOID(x -> 2/(1+Math.exp(-x))-1),
  TANH(x -> Math.tanh(x)),
  RELU(x -> x>0?x:0),
  LEAKYRELU(x -> x>0?x:0.1*x),
  ELU(x -> x>0?x:(Math.exp(x)-1)),
  SOFTPLUS(x -> Math.log(1+Math.exp(x))),
  SWISH(x -> x/(1+Math.exp(-x))),
  HARDSWISH(x -> x*Math.max(0, Math.min(1, 0.2 * x + 0.5))),
  MISH(x -> x * Math.tanh(Math.log1p(Math.exp(x)))),
  SOFTSIGN(x -> x / (1 + Math.abs(x))),
  GELU(x -> x/2*(1+Math.tanh(Math.sqrt(2/Math.PI)*(x+0.044715*Math.pow(x,3))))),
  BINARYSTEP(x -> x>=0?1.0:0),
  LINEAR(x -> x),
  ARCTAN(x -> Math.atan(x)),
  GAUSSIAN(x -> Math.exp(-x*x)),
  SINE(x -> Math.sin(x));

  private final Function<Double,Double> func;
  ActivationFunction(Function<Double,Double> f) {
    func = f;
  }
  public Function<Double,Double> get() {
    return func;
  }
}
